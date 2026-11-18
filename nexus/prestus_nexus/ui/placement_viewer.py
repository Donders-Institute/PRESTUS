"""Placement viewer: orthogonal T1 slices with segmentation overlay, click-to-place
transducer and focus, bowl outline and distance checks (manual placement)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt6.QtCore import QObject, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QColor, QKeySequence, QShortcut
from PyQt6.QtWidgets import (QButtonGroup, QComboBox, QDoubleSpinBox, QFileDialog, QFrame,
                             QHBoxLayout, QInputDialog, QLabel, QMenu, QMessageBox, QPushButton, QSlider,
                             QToolButton, QVBoxLayout, QWidget, QWidgetAction)

from .. import checks, localite as lc, multitransducer as mt, neuronav_import as ni, paths, placement as pl
from ..bridge import BridgeError, as_floats, as_list
from ..heuristic_import import candidates as heuristic_candidates, selected as heuristic_selected
from ..model import TX, ConfigModel, Issue, is_empty
from ..outputs import dir_placement_run, placement_record
from ..undo import History
from . import icons, style
from .readout import muted, readout, status
from .empty_state import EmptyStack
from .ortho_view import load_volume_async


class _LoadRelay(QObject):
    """Hands a finished T1 / segmentation load from the worker thread to the GUI thread."""
    done = pyqtSignal(object)           # the load's key

try:        # matplotlib is optional: the page explains what is missing
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False


def _mm_status(v: float, warn_below: float | None = None) -> str:
    """Signed distance [mm], red when negative (amber below warn_below), 'n/a' when unknown."""
    if not np.isfinite(v):
        return muted("n/a")
    level = "bad" if v < 0 else ("warn" if warn_below is not None and v < warn_below else "ok")
    return status(f"{v:+.1f} mm", level)


QC_PATH = "placement viewer"     # Issue.path of the viewer's checks (MainWindow opens the viewer)
DEPTH_TOL_MM = 2.0               # target vs focal depth setting: red / a warning above this


# SimNIBS charm labels (final_tissues.nii.gz) and display colours
TISSUE = {1: ("white matter", "#f2f2f2"), 2: ("grey matter", "#9e9e9e"), 3: ("CSF", "#4fa3e0"),
          4: ("bone", "#f5deb3"), 5: ("scalp", "#e8a07a"), 6: ("eye", "#7ec8e3"),
          7: ("compact bone", "#fff3c4"), 8: ("spongy bone", "#d9b779"), 9: ("blood", "#c0392b"),
          10: ("muscle", "#b5651d")}
VIEWS = [("Axial", 2, (0, 1)), ("Coronal", 1, (0, 2)), ("Sagittal", 0, (1, 2))]
AXIS_LABELS = {0: ("L", "R"), 1: ("P", "A"), 2: ("I", "S")}
_OVERHEAD_FRAC = 0.12          # approx. fraction of row height spent on title/axis-label chrome
_LEGEND_FRAC = 0.1             # fraction of the figure the marker legend under the slices takes


class PlacementPage(QWidget):
    changed = pyqtSignal()
    transducer_selected = pyqtSignal(int)   # the transducer clicks now edit
    navigate = pyqtSignal(str)        # config key to show (MainWindow.goto), from the empty state

    def __init__(self, model: ConfigModel, parent: QWidget | None = None, history: History | None = None):
        super().__init__(parent)
        self.position_source = ""         # where a non-manual position came from (banner)
        # Import trajectory of a Localite XML asks MATLAB for its series (MainWindow sets these)
        self.bridge_getter = None         # MainWindow.bridge: the shared MATLAB session
        self.saved_config_getter = None   # saves the study config, returns its path (or None)
        # (file, ti) -> (target, position, trans voxel, focus voxel) from the last import
        self._localite_file_pos: dict[tuple[str, int], tuple] = {}
        self.model = model
        self.ti = 0
        self.t1: pl.Volume | None = None
        self.seg: pl.Volume | None = None
        self._loaded_for: tuple | None = None
        self.cursor = np.zeros(3)
        self.direction: np.ndarray | None = None     # focus -> transducer, when no transducer yet
        self._slice_ratio: float | None = None       # ideal slice-row height per pixel of width

        self._t1_error = ""                          # why the T1 could not be read, if so
        v = QVBoxLayout(self)                        # the sub-tab label names the page
        content = QWidget()
        cv = QVBoxLayout(content)
        cv.setContentsMargins(0, 0, 0, 0)
        self._content = content
        self.stack = EmptyStack(content)
        v.addWidget(self.stack, 1)

        # One status banner (info / warning / error) instead of lines of mixed colours;
        # 'Switch to manual placement' sits in it because it only applies to what it says.
        self.banner_box = QFrame()
        self.banner_box.setObjectName("ViewerBanner")
        bb = QHBoxLayout(self.banner_box)
        bb.setContentsMargins(10, 6, 6, 6)
        self.banner = QLabel()
        self.banner.setWordWrap(True)
        self.banner.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        bb.addWidget(self.banner, 1)
        self.use_manual = QPushButton("Switch to m&anual placement")
        self.use_manual.setToolTip("Set placement.mode to manual, so clicks and the Place actions "
                                   "write the positions to the study config.")
        self.use_manual.clicked.connect(self._use_manual)
        bb.addWidget(self.use_manual, 0, Qt.AlignmentFlag.AlignVCenter)
        self.banner_box.hide()
        cv.addWidget(self.banner_box)

        bar = QHBoxLayout()
        bar.setSpacing(6)
        self.tx_combo = QComboBox()
        self.tx_combo.setToolTip("The transducer that clicks and the Place actions edit")
        self.tx_combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.tx_combo.setMinimumContentsLength(len("Transducer 10"))     # never clipped
        self.tx_combo.currentIndexChanged.connect(self._select_transducer)
        bar.addWidget(self.tx_combo)
        self.tx_add = QPushButton("+")
        self.tx_add.setToolTip("Add a transducer to this study")
        self.tx_add.clicked.connect(self._add_transducer)
        self.tx_remove = QPushButton("−")
        self.tx_remove.setToolTip("Remove the selected transducer from this study")
        self.tx_remove.clicked.connect(self._remove_transducer)
        for b in (self.tx_add, self.tx_remove):
            b.setProperty("compact", True)          # the app's 5px 12px padding clips a glyph
            b.setStyleSheet("padding: 2px 0;")
            b.setFixedWidth(30)
            bar.addWidget(b)
        bar.addSpacing(10)

        click = QLabel("Click sets")
        click.setObjectName("CardSub")
        bar.addWidget(click)
        self.click_mode = QFrame()
        self.click_mode.setObjectName("Segmented")
        seg = QHBoxLayout(self.click_mode)
        seg.setContentsMargins(2, 2, 2, 2)
        seg.setSpacing(2)
        self.rb_focus = QToolButton(text="F&ocus")
        self.rb_focus.setToolTip("Left-click on a slice places the focus (the target) there")
        self.rb_trans = QToolButton(text="T&ransducer")
        self.rb_trans.setToolTip("Left-click on a slice places the transducer (bowl apex) there")
        self.rb_browse = QToolButton(text="&Browse")
        self.rb_browse.setToolTip("Left-click only moves the crosshair; nothing is written to the "
                                  "study config")
        self.click_group = QButtonGroup(self)
        self.click_group.setExclusive(True)
        for b in (self.rb_focus, self.rb_trans, self.rb_browse):
            b.setCheckable(True)
            b.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
            self.click_group.addButton(b)
            seg.addWidget(b)
        self.rb_focus.setChecked(True)
        bar.addWidget(self.click_mode)
        bar.addSpacing(10)

        # undo / redo: the window's edit history (Edit ▸ Undo / Redo), so these buttons and
        # Ctrl+Z step through edits made on any tab in the order they were made
        own = history is None
        self.history = History(model, parent=self) if own else history
        self.history.changed.connect(self._update_undo_buttons)
        self.history.restored.connect(self._history_restored)
        self.undo_btn = QPushButton(icons.icon("undo"), "")
        self.undo_btn.clicked.connect(self.undo)
        self.redo_btn = QPushButton(icons.icon("redo"), "")
        self.redo_btn.clicked.connect(self.redo)
        for btn, keys in ((self.undo_btn, QKeySequence.StandardKey.Undo),
                          (self.redo_btn, QKeySequence.StandardKey.Redo)):
            btn.setProperty("compact", True)
            if own:                        # in the main window, Edit ▸ Undo / Redo own the keys
                sc = QShortcut(QKeySequence(keys), self, btn.click)
                sc.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            bar.addWidget(btn)
        self._update_undo_buttons()
        bar.addSpacing(10)

        b_snap = QPushButton("&Snap to scalp")
        b_snap.setToolTip("Keep the beam direction and move the transducer so its exit plane sits "
                          "the 'gap' distance (More menu) outside the scalp.")
        b_snap.clicked.connect(self._snap)
        b_match = QPushButton("Set focal &depth to target")
        b_match.setToolTip("Set the focal depth setting (focal_distance_ep) to this placement's "
                           "target depth, exit plane to target.")
        b_match.clicked.connect(self._match_focal_distance)
        self.more_btn = QPushButton("&More ▾")
        self.more_btn.setToolTip("Less common placement actions: import a trajectory, put the "
                                 "transducer at the focal depth, scalp gap")
        more = QMenu(self.more_btn)
        a_traj = more.addAction("Import trajectory…")
        a_traj.setToolTip("Place focus and transducer from a Brainsight target export (NIfTI "
                          "coordinates), a 3D Slicer transform, or a k-Plan plan; or use a series "
                          "of a Localite export (localite placement mode; needs MATLAB).")
        a_traj.triggered.connect(self._import_trajectory)
        a_axis = more.addAction("Transducer at focal depth")
        a_axis.setToolTip("Keep the beam direction and put the transducer where the focal depth "
                          "setting lands on the target.")
        a_axis.triggered.connect(self._place_on_axis)
        more.addSeparator()
        self.gap = QDoubleSpinBox(minimum=-20, maximum=50, value=5, singleStep=0.5, suffix=" mm")
        self.gap.setToolTip("Distance between scalp and exit plane that 'Snap to scalp' leaves")
        more.addAction(self._menu_row("Scalp gap", self.gap, more))
        more.setToolTipsVisible(True)
        self.more_btn.clicked.connect(lambda: self._popup(more, self.more_btn))
        self._more_menu = more
        for w in (b_snap, b_match, self.more_btn):
            bar.addWidget(w)
        bar.addStretch(1)

        self.view_btn = QPushButton("Vie&w ▾")
        self.view_btn.setToolTip("Overlay opacity, acoustic field-of-view size and centring")
        view = QMenu(self.view_btn)
        b_centre = view.addAction("Centre on focus")
        b_centre.setToolTip("Move the crosshair (and so the three slices) to the focus")
        b_centre.triggered.connect(self._centre)
        view.addSeparator()
        self.alpha = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=100, value=35)
        self.alpha.setMinimumWidth(140)
        self.alpha.setToolTip("Opacity of the tissue segmentation drawn over the T1")
        self.alpha.valueChanged.connect(lambda _v: self._request_draw())
        view.addAction(self._menu_row("Tissue overlay", self.alpha, view))
        self.fov_diam = QDoubleSpinBox(minimum=0, maximum=300, value=0, singleStep=1, suffix=" mm")
        self.fov_diam.setSpecialValueText("auto")
        self.fov_diam.setToolTip("grid.acoustic_fov_diameter_mm — lateral width of the acoustic "
                                 "simulation crop, centred on the focus. 'auto' = outermost "
                                 "element diameter.")
        self.fov_diam.editingFinished.connect(self._fov_edited)
        view.addAction(self._menu_row("Acoustic FOV ⌀", self.fov_diam, view))
        self.fov_len = QDoubleSpinBox(minimum=0, maximum=400, value=0, singleStep=1, suffix=" mm")
        self.fov_len.setSpecialValueText("auto")
        self.fov_len.setToolTip("grid.acoustic_fov_length_mm — axial length of the acoustic "
                                "simulation crop, starting at the transducer. 'auto' = focal "
                                "depth + grid.z_beyond_focus_mm.")
        self.fov_len.editingFinished.connect(self._fov_edited)
        view.addAction(self._menu_row("Acoustic FOV length", self.fov_len, view))
        view.addSeparator()
        self.grid_z_action = view.addAction("Grid z axis (grid.mode)")
        self.grid_z_action.setCheckable(True)
        self.grid_z_action.setToolTip(
            "The direction PRESTUS rotates onto the simulation grid's z axis under grid.mode, "
            "drawn through the mean focus. On by default with more than one transducer. The "
            "cropped grid box needs the head masks and is only known after preprocessing.")
        self.grid_z_action.toggled.connect(self._grid_z_toggled)
        view.setToolTipsVisible(True)
        self.view_btn.clicked.connect(lambda: self._popup(view, self.view_btn))
        self._view_menu = view
        bar.addWidget(self.view_btn)
        self.coords_btn = QPushButton("&Coordinates…")
        self.coords_btn.setCheckable(True)
        self.coords_btn.setToolTip("Type RAS, voxel, or MNI coordinates to go to or to set the "
                                   "focus or transducer")
        bar.addWidget(self.coords_btn)
        cv.addLayout(bar)

        self.coords_row = QWidget()
        coords = QHBoxLayout(self.coords_row)
        coords.setContentsMargins(0, 0, 0, 0)
        self.coord_space = QComboBox()
        self.coord_space.addItem("RAS mm", "ras")
        self.coord_space.addItem("Voxel", "voxel")
        self.coord_space.addItem("MNI mm", "mni")
        self.coord_space.setCurrentIndex({"ras": 0, "voxel": 1, "mni": 2}.get(self._space(), 1))
        self.coord_space.setToolTip("MNI mm is converted to this subject's space with SimNIBS's "
                                    "mni2subject_coords; RAS and voxel need no conversion. Both "
                                    "the nonlinear and affine transform come from SimNIBS (charm) "
                                    "having already been run for this subject -- the m2m folder "
                                    "must exist either way, setting startup.simnibs_bin_path alone "
                                    "is not enough. For the focus this uses the nonlinear warp "
                                    "(most accurate in-brain); for the transducer it uses the "
                                    "affine, which also still works if that m2m folder's nonlinear "
                                    "warp specifically failed or was never produced.")
        coords.addWidget(self.coord_space)
        self.cx = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        self.cy = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        self.cz = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        for label, spin in (("X", self.cx), ("Y", self.cy), ("Z", self.cz)):
            coords.addWidget(QLabel(label))
            coords.addWidget(spin)
        b_goto = QPushButton("&Go to")
        b_goto.setToolTip("Move the crosshair to these coordinates; writes nothing")
        b_goto.clicked.connect(self._goto_coords)
        b_set_focus = QPushButton("Set focus")
        b_set_focus.setToolTip("Write these coordinates as the focus of the selected transducer")
        b_set_focus.clicked.connect(lambda: self._set_coords("focus_pos"))
        b_set_trans = QPushButton("Set transducer")
        b_set_trans.setToolTip("Write these coordinates as the transducer (bowl apex) position")
        b_set_trans.clicked.connect(lambda: self._set_coords("trans_pos"))
        for w in (b_goto, b_set_focus, b_set_trans):
            coords.addWidget(w)
        coords.addStretch(1)
        self.coords_row.setVisible(False)
        self.coords_btn.toggled.connect(self.coords_row.setVisible)
        cv.addWidget(self.coords_row)

        body = QHBoxLayout()
        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.6), tight_layout=True)
            from .export import subject_stem
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="placement",
                                            filename_hint=subject_stem(lambda: self.model, "placement"))
            self.axes = [self.fig.add_subplot(1, 3, i + 1) for i in range(3)]
            self.canvas.mpl_connect("button_press_event", self._clicked)
            from .tooltip import MarkerTips
            self.marker_tips = MarkerTips(self.canvas)     # hover a marker: its transducer's name
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            body.addWidget(self.canvas, 3)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.info.setMinimumWidth(280)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.info.setToolTip(
            "Exit plane → scalp: along the beam axis, from the SimNIBS labels.\n"
            "Rim → scalp: from 36 points on the exit-plane rim along the beam; a negative minimum "
            "means the transducer cuts into the head, a large SD uneven contact.\n"
            "Bone on the beam: 37 rays from the aperture (centre + 3 rings) to the target.")
        body.addWidget(self.info, 1)
        cv.addLayout(body, 1)
        self._body = body
        # the how-to stays visible; the pointer to the Beam axis tab moves to the tooltip
        hint = QLabel("Click a slice to set the chosen point; scroll to step through slices.<br>"
                      "Positions are saved to the study config at once.")
        hint.setToolTip("The Beam axis tab shows the medium and the free-water profile along the beam.")
        hint.setObjectName("CardSub")
        hint.setWordWrap(True)
        cv.addWidget(hint)
        self._restyle()

    @staticmethod
    def _menu_row(label: str, widget: QWidget, menu: QMenu) -> QWidgetAction:
        """A labelled control inside a drop-down menu."""
        row = QWidget()
        h = QHBoxLayout(row)
        h.setContentsMargins(12, 4, 12, 4)
        text = QLabel(label)
        text.setMinimumWidth(130)
        h.addWidget(text)
        h.addWidget(widget, 1)
        act = QWidgetAction(menu)
        act.setDefaultWidget(row)
        return act

    @staticmethod
    def _popup(menu: QMenu, button: QWidget) -> None:
        menu.popup(button.mapToGlobal(button.rect().bottomLeft()))

    def _restyle(self) -> None:
        """Palette-dependent looks (segmented click mode), re-applied on refresh so a theme
        switch is followed."""
        c = style.active()
        self.click_mode.setStyleSheet(
            f"QFrame#Segmented {{ background: {c['surface']}; border: 1px solid {c['line2']}; "
            f"border-radius: 6px; }}"
            f"QFrame#Segmented QToolButton {{ border: none; border-radius: 4px; padding: 3px 10px; "
            f"color: {c['ink2']}; }}"
            f"QFrame#Segmented QToolButton:hover {{ background: {c['surface3']}; }}"
            f"QFrame#Segmented QToolButton:checked {{ background: {c['brand']}; color: white; "
            f"font-weight: 600; }}")

    def _set_banner(self, level: str, lines: list[str]) -> None:
        """One compact banner: level 'info' | 'warning' | 'error' picks the palette colour."""
        if not lines:
            self.banner_box.hide()
            return
        c = style.active()
        col = QColor(c[{"error": "red", "warning": "amber"}.get(level, "info")])
        rgba = lambda a: f"rgba({col.red()}, {col.green()}, {col.blue()}, {a})"   # noqa: E731
        self.banner_box.setStyleSheet(
            f"QFrame#ViewerBanner {{ background: {rgba(0.10)}; border: 1px solid {rgba(0.45)}; "
            f"border-radius: 6px; }}"
            f"QFrame#ViewerBanner QLabel {{ color: {c['ink']}; background: transparent; }}")
        self.banner.setText("<br>".join(lines))
        self.banner_box.show()

    def _show_empty(self) -> None:
        """No T1 to draw on: one panel that says what is missing and opens that field."""
        m = self.model
        sid, anat = m.get("subject_id"), m.get("path.anat")
        if not HAVE_MPL:
            self.stack.show_empty("The viewer needs matplotlib and nibabel",
                                  "Install them in the Nexus's Python environment and restart.")
        elif not isinstance(sid, (int, float)) or isinstance(sid, bool):
            self.stack.show_empty("No subject selected",
                                  "The viewer shows one subject's T1 image; set subject_id first.",
                                  "Set subject…", lambda: self.navigate.emit("subject_id"))
        elif is_empty(anat):
            self.stack.show_empty("No T1 image yet",
                                  "The viewer draws the transducer and focus on the subject's T1 "
                                  "image. Set the anatomy folder (path.anat) on Setup ▸ Paths.",
                                  "Set T1 image…", lambda: self.navigate.emit("path.anat"))
        elif self._t1_error:
            self.stack.show_empty("Cannot read the T1 image", self._t1_error,
                                  "Check T1 image…", lambda: self.navigate.emit("path.anat"))
        else:
            pattern = m.get("path.t1_pattern") or ""
            self.stack.show_empty(f"No T1 image found for sub-{int(sid):03d}",
                                  f"Nothing matches path.t1_pattern '{pattern}' in {anat}. "
                                  "Check the anatomy folder and the file pattern.",
                                  "Check T1 image…", lambda: self.navigate.emit("path.anat"))

    def resizeEvent(self, event) -> None:      # noqa: N802
        super().resizeEvent(event)
        self._apply_row_cap()

    # ── state from the model ────────────────────────────────────────────
    def _space(self) -> str:
        return str(self.model.get("placement.manual.coordinate_space") or "voxel")

    def _editable(self) -> bool:
        return self.model.get("placement.mode") in ("manual", None, "") and self._space() in ("voxel", "ras")

    # The T1 and segmentation are read on ortho_view's loader thread (shared cache, one read
    # per file however many pages ask). refresh() does not wait: it shows 'Loading' and draws
    # when the read is done. Callers that need the images at once (ensure_images, the Beam
    # axis and 3D pages, the checks) wait for the same read.
    _pending_load: tuple | None = None      # (key, {'t1': (path, Future), 'seg': ...})
    _load_hooked: object = None
    _t1: pl.Volume | None = None
    _seg: pl.Volume | None = None

    # Reading t1 / seg waits for a read refresh() started, so every user of the images
    # (Beam axis, 3D, PlanTUS pages, the checks) sees them as before.
    @property
    def t1(self) -> pl.Volume | None:
        self._finish_load()
        return self._t1

    @t1.setter
    def t1(self, vol: pl.Volume | None) -> None:
        self._t1 = vol

    @property
    def seg(self) -> pl.Volume | None:
        self._finish_load()
        return self._seg

    @seg.setter
    def seg(self, vol: pl.Volume | None) -> None:
        self._seg = vol

    def _load_images(self, wait: bool = True) -> None:
        m = self.model
        sid = m.get("subject_id")
        anat, pattern = m.get("path.anat"), m.get("path.t1_pattern")
        seg_root = m.get("path.seg") or anat
        key = (sid, anat, pattern, seg_root)
        if key != self._loaded_for:
            self._loaded_for = key
            self.t1 = self.seg = None
            self._overlay_cache = {}            # per-image caches (keyed by object id) start over
            self._axes_for = None
            self._t1_error = ""
            self._pending_load = None
            if not isinstance(sid, (int, float)) or is_empty(anat):
                return
            futures = {}
            t1 = pl.find_t1(anat, pattern or "", int(sid))
            seg = pl.find_segmentation(seg_root, int(sid)) if not is_empty(seg_root) else None
            for role, path in (("t1", t1), ("seg", seg)):
                if path is None:
                    continue
                try:
                    futures[role] = (path, load_volume_async(path))
                except OSError as exc:
                    if role == "t1":
                        self._t1_error = f"{path}: {exc}"
            self._pending_load = (key, futures) if futures else None
        if self._pending_load is None:
            return
        key, futures = self._pending_load
        if wait or all(f.done() for _p, f in futures.values()):
            self._apply_load()
        elif self._load_hooked is not key:
            self._load_hooked = key
            relay = self.__dict__.get("_relay")
            if relay is None:
                relay = self._relay = _LoadRelay(self)
                relay.done.connect(self._load_done, Qt.ConnectionType.QueuedConnection)
            for _p, f in futures.values():
                def done(_f, k=key, r=relay):
                    try:
                        r.done.emit(k)
                    except RuntimeError:          # page closed meanwhile
                        pass
                f.add_done_callback(done)

    def _apply_load(self) -> None:
        """Take the finished (or waited-for) reads of the pending load."""
        _key, futures = self._pending_load
        self._pending_load = None
        for role, (path, fut) in futures.items():
            try:
                vol = fut.result()
            except Exception as exc:
                if role == "t1":
                    self._t1_error = f"{path}: {exc}"
                vol = None
            if role == "t1":
                self.t1 = vol
                if vol is not None:
                    self.cursor = vol.to_ras(np.array(vol.shape) / 2)[0]
            else:
                self.seg = vol

    def _load_done(self, key) -> None:
        if self._pending_load is None or self._pending_load[0] != key:
            return                                # applied already, or superseded
        if all(f.done() for _p, f in self._pending_load[1].values()):
            self._apply_load()
            self.refresh()

    def _finish_load(self) -> None:
        """Wait for a T1 / segmentation read refresh() started (positions need the T1)."""
        if self._pending_load is not None:
            self._load_images()

    def images_loading(self) -> bool:
        """A T1 / segmentation read is still running (refresh shows 'Loading')."""
        return self._pending_load is not None

    def positions(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        """(transducer, focus) in RAS mm for the current transducer, if known."""
        return self._positions_for(self.ti)

    # One draw asks for every transducer's positions and geometry several times, and in the
    # Localite / heuristic modes each answer reads files on the project volume; so answers
    # are kept while the model's values, the T1 and the files' refresh() stay the same.
    _memo_key: tuple | None = None

    def _memo(self, name: str, ti: int | None, compute):
        key = (repr((self.model.overrides, self.model.tx)), id(self.t1), id(self.seg))
        if self._memo_key != key:
            self._memo_key, self._memo_store = key, {}
        if (name, ti) not in self._memo_store:
            self._memo_store[(name, ti)] = compute()
        return self._memo_store[(name, ti)]

    def _forget(self) -> None:
        """Read positions and geometry again on next use (files may have changed)."""
        self._memo_key = None

    def _positions_for(self, ti: int) -> tuple[np.ndarray | None, np.ndarray | None]:
        def compute():
            before = self.position_source
            trans, focus = self._read_positions(ti)
            source = self.position_source if self.position_source != before else None
            self.position_source = before
            return trans, focus, source
        trans, focus, source = self._memo("positions", ti, compute)
        if source is not None:
            self.position_source = source
        return (None if trans is None else trans.copy()), (None if focus is None else focus.copy())

    def _read_positions(self, ti: int) -> tuple[np.ndarray | None, np.ndarray | None]:
        """(transducer, focus) in RAS mm for transducer ti, if known.

        Manual placement reads transducer.trans_pos/focus_pos directly. Localite and
        heuristic placement leave those empty (PRESTUS resolves them at run time from
        ingested markers / the heuristic search) -- for those two modes, fall back to
        what neuronav_gui / transducer_positioning already wrote to disk, so the viewer
        isn't blank for every non-manual mode."""
        self._finish_load()
        out = []
        for key in ("trans_pos", "focus_pos"):
            v = self.model.get(f"{TX}.{key}", ti)
            if not (isinstance(v, list) and len(v) == 3) or self.t1 is None:
                out.append(None)
            elif self._space() == "ras":
                out.append(np.asarray(v, float))
            elif self._space() == "voxel":
                out.append(pl.prestus_voxel_to_ras(self.t1, v))
            else:
                out.append(None)          # MNI: needs SimNIBS to map to subject space
        if out[0] is None and out[1] is None:
            mode = self.model.get("placement.mode")
            found, self.position_source = (None, None), ""
            if mode == "localite":
                found = self._localite_positions_for(ti)
                self.position_source = ("the Localite file import"
                                        if not is_empty(self.model.get("placement.localite.file"))
                                        else "neuronav_gui")
            elif mode == "heuristic":
                found, self.position_source = self._heuristic_positions_for(ti), "transducer_positioning"
            if all(v is None for v in found) and mode not in ("manual", None, ""):
                found = self._recorded_positions_for(ti)
                rec = placement_record(self.model)
                self.position_source = f"the last run's placement stage ({rec.get('created', '')})" \
                    if rec else ""
            return found
        return out[0], out[1]

    def target_depth_ep(self, ti: int) -> float | None:
        """Planned target depth of transducer ti: exit plane -> focus_pos along the beam [mm],
        from the positions this viewer shows (loads the T1 if the tab was never opened)."""
        self._load_images()
        source = self.position_source
        trans, focus = self._positions_for(ti)
        self.position_source = source           # keep the banner of the shown transducer
        bowl = self._bowl_for(ti)
        if trans is None or focus is None or bowl is None:
            return None
        return float(np.linalg.norm(focus - trans)) - bowl.offset

    def _recorded_positions_for(self, ti: int) -> tuple[np.ndarray | None, np.ndarray | None]:
        """Positions a previous run's placement stage resolved (sub-XXX_placement<affix>.json,
        write_placement_positions.m) -- covers every placement mode."""
        rec = placement_record(self.model)
        if rec is None or self.t1 is None or rec.get("mode") != self.model.get("placement.mode"):
            return None, None
        txs = rec["transducers"]
        if ti >= len(txs) or not isinstance(txs[ti], dict):
            return None, None
        out = []
        for key in ("trans_pos", "focus_pos"):
            v = txs[ti].get(key)
            out.append(pl.prestus_voxel_to_ras(self.t1, v) if isinstance(v, list) and len(v) == 3 else None)
        return out[0], out[1]

    def _localite_positions_for(self, ti: int) -> tuple[np.ndarray | None, np.ndarray | None]:
        file = self.model.get("placement.localite.file")
        if not is_empty(file):
            # placement.localite.file: the positions the import got from MATLAB, while the
            # series chosen then is still the one configured
            hit = self._localite_file_pos.get((str(file), ti))
            sel = (repr(self.model.get("placement.localite.target")),
                   repr(self.model.get("placement.localite.position")))
            if hit is None or hit[:2] != sel or self.t1 is None:
                return None, None
            return tuple(None if v is None else pl.prestus_voxel_to_ras(self.t1, v) for v in hit[2:])
        post = self.model.get("path.localite_post") or self.model.get("path.localite")
        sid = self.model.get("subject_id")
        session = self.model.get("placement.localite.session")
        target = self.model.get("placement.localite.target")
        if is_empty(post) or is_empty(session) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None, None
        if isinstance(target, list):
            if ti >= len(target):
                return None, None
            target = target[ti]
        elif ti > 0:
            return None, None            # single target string: only transducer 1 is known
        if is_empty(target):
            return None, None
        try:
            entry = ni.coordinates_for(str(post), session, str(target), int(sid))
        except OSError:
            return None, None
        if entry is None:
            return None, None
        return entry.trans_ras, entry.focus_ras

    def _heuristic_positions_for(self, ti: int) -> tuple[np.ndarray | None, np.ndarray | None]:
        if ti > 0:
            return None, None            # transducer_positioning writes one target per run
        base = dir_placement_run(self.model, "heuristic")
        sid = self.model.get("subject_id")
        target = self.model.get("placement.heuristic.target_name")
        if base is None or is_empty(target) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None, None
        sel = heuristic_selected(base, int(sid), str(target))
        if sel is None:
            return None, None
        return sel.trans_ras, sel.focus_ras

    def _write(self, key: str, ras: np.ndarray) -> None:
        self._finish_load()
        self._record("Move " + ("focus" if key == "focus_pos" else "transducer"))
        if self._space() == "ras":
            value = [round(float(c), 2) for c in ras]
        else:
            value = pl.ras_to_prestus_voxel(self.t1, ras)
        self.model.set(f"{TX}.{key}", value, self.ti)
        self.changed.emit()

    # ── undo / redo ─────────────────────────────────────────────────────
    def _record(self, label: str = "Placement change") -> None:
        """Make the values this user action writes (however many) one step of the window's
        edit history."""
        self.history.group_pass(label)

    def undo(self) -> None:
        self.history.undo()

    def redo(self) -> None:
        self.history.redo()

    def _history_restored(self, _label: str) -> None:
        self._cursor_set_by_user = False
        if self.history.parent() is self:  # standalone; in the main window every page refreshes
            self.changed.emit()
            self.refresh()

    def _update_undo_buttons(self) -> None:
        h = self.history
        self.undo_btn.setEnabled(h.can_undo())
        self.redo_btn.setEnabled(h.can_redo())
        self.undo_btn.setToolTip(f"{h.undo_text()} (Ctrl+Z); one history for all tabs")
        self.redo_btn.setToolTip(f"{h.redo_text()} (Ctrl+Shift+Z)")

    def _bowl(self) -> pl.Bowl | None:
        return self._bowl_for(self.ti)

    def _transducers(self) -> list[dict]:
        """model.effective_transducers(), merged once per model state."""
        return self._memo("transducers", None, self.model.effective_transducers)

    def _bowl_for(self, ti: int) -> pl.Bowl | None:
        txs = self._transducers()
        return pl.Bowl.from_transducer(txs[ti]) if ti < len(txs) else None

    def _heuristic_candidate_points(self) -> tuple[np.ndarray, np.ndarray]:
        """RAS mm positions and prop_intersect of every heuristic candidate (tpars_*.csv),
        for the current subject/target; empty arrays if not heuristic mode / not found."""
        if self.model.get("placement.mode") != "heuristic" or self.t1 is None:
            return np.empty((0, 3)), np.empty(0)
        base = dir_placement_run(self.model, "heuristic")
        sid = self.model.get("subject_id")
        target = self.model.get("placement.heuristic.target_name")
        if base is None or is_empty(target) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return np.empty((0, 3)), np.empty(0)
        cands = heuristic_candidates(base, int(sid), str(target))
        if not cands:
            return np.empty((0, 3)), np.empty(0)
        pts = self.t1.to_ras(np.array([c.trans_vox for c in cands]) - 1)
        return pts, np.array([c.prop_intersect for c in cands])

    def _resolved_fov(self, trans, focus, bowl: pl.Bowl | None) -> tuple[float | None, float | None, bool, bool]:
        """(diameter_mm, length_mm, diameter_is_auto, length_is_auto) for the selected
        transducer, resolving grid.acoustic_fov_diameter_mm/length_mm 'auto' the way
        acoustic_grid_fov.m does (doc_partial-fov.md)."""
        d, l = self.model.get("grid.acoustic_fov_diameter_mm"), self.model.get("grid.acoustic_fov_length_mm")
        diam_auto = not (isinstance(d, (int, float)) and not isinstance(d, bool) and d > 0)
        len_auto = not (isinstance(l, (int, float)) and not isinstance(l, bool) and l > 0)
        diameter = (bowl.aperture if bowl is not None else None) if diam_auto else float(d)
        if not len_auto:
            length = float(l)
        elif trans is not None and focus is not None:
            zbf = self.model.get("grid.z_beyond_focus_mm")
            zbf = float(zbf) if isinstance(zbf, (int, float)) and not isinstance(zbf, bool) else 40.0
            length = float(np.linalg.norm(focus - trans)) + zbf
        else:
            length = None
        return diameter, length, diam_auto, len_auto

    def _expected_bowl_distance(self, ti: int | None = None) -> float | None:
        ti = self.ti if ti is None else ti
        tx = self._transducers()[ti]
        bowl = self._bowl_for(ti)
        ep, fb = tx.get("focal_distance_ep"), tx.get("focal_distance_bowl")
        if isinstance(fb, (int, float)) and not isinstance(fb, bool):
            return float(fb)
        if isinstance(ep, (int, float)) and not isinstance(ep, bool) and bowl is not None:
            return float(ep) + bowl.offset
        return None

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        self._forget()
        self.tx_combo.blockSignals(True)
        self.tx_combo.clear()
        for i in range(self.model.n_transducers):
            self.tx_combo.addItem(f"Transducer {i + 1}")
            self.tx_combo.setItemData(i, self.model.transducer_tooltip(i), Qt.ItemDataRole.ToolTipRole)
        self.ti = min(self.ti, self.model.n_transducers - 1)
        self.tx_combo.setCurrentIndex(self.ti)
        self.tx_combo.setToolTip(f"{self.model.transducer_tooltip(self.ti)}\n"
                                 "The transducer that clicks and the Place actions edit")
        self.tx_combo.blockSignals(False)
        self.tx_remove.setEnabled(self.model.n_transducers > 1)
        self._load_images(wait=False)
        if self._grid_z_user is None:
            self.grid_z_action.blockSignals(True)
            self.grid_z_action.setChecked(self.model.n_transducers > 1)
            self.grid_z_action.blockSignals(False)

        self._restyle()
        if self.images_loading():
            path = next(iter(self._pending_load[1].values()))[0]
            self.stack.show_empty("Loading the T1 image…", str(path))
            self._slice_ratio = None
            return
        if self.t1 is None:
            self._show_empty()
            self._slice_ratio = None
            return
        self.stack.show_content()
        mode = self.model.get("placement.mode")
        msgs = []
        level = "info"
        self.use_manual.setVisible(mode not in ("manual", None, ""))
        if mode not in ("manual", None, ""):
            r_trans, r_focus = self.positions()
            if r_trans is not None or r_focus is not None:
                source = self.position_source or "PRESTUS"
                msgs.append(f"placement.mode is <b>{mode}</b>: showing the position written by {source} "
                            "(read-only here; PRESTUS resolves it again at run time).")
            else:
                level = "warning"
                where = {"localite": "placement.localite.file (Import trajectory), or path.localite_post "
                                     "with placement.localite.session/.target",
                         "heuristic": "the heuristic placement run folder (tpars_*.csv)"
                         }.get(mode, "")
                where = (where + " and " if where else "") + "sub-XXX_placement*.json in io.dir_output"
                msgs.append(f"placement.mode is <b>{mode}</b>: no resolved position found yet "
                            f"(checked {where}). Run the placement stage once "
                            "(e.g. with only modules.run_transducer_placement), or switch to manual placement.")
        if self._space() == "mni":
            level = "warning"
            msgs.append("coordinate_space is <b>mni</b>: positions need SimNIBS to map to this "
                        "subject and are not shown. Use 'ras' or 'voxel' to place here.")
        if self.seg is None:
            msgs.append("No segmentation (m2m_sub-XXX/final_tissues.nii.gz): no tissue overlay "
                        "and no scalp snapping.")
        self._set_banner(level, msgs)
        self.fov_diam.blockSignals(True)
        self.fov_len.blockSignals(True)
        self.fov_diam.setValue(self._mm_or(self.model.get("grid.acoustic_fov_diameter_mm")))
        self.fov_len.setValue(self._mm_or(self.model.get("grid.acoustic_fov_length_mm")))
        self.fov_diam.blockSignals(False)
        self.fov_len.blockSignals(False)
        trans, focus = self.positions()
        if focus is not None and not self._cursor_set_by_user:
            self.cursor = focus.copy()
        if self.window().isVisible() and not self.isVisible():
            self._stale = True              # another sub-tab is showing: draw when shown
            return
        self._draw()

    _stale = False

    def showEvent(self, event) -> None:        # noqa: N802
        super().showEvent(event)
        if self._stale:
            self._stale = False
            self._draw()

    @staticmethod
    def _mm_or(v, default: float = 0.0) -> float:
        return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else default

    def _fov_edited(self) -> None:
        self._record("Change acoustic FOV")
        self.model.set("grid.acoustic_fov_diameter_mm", self.fov_diam.value() or [])
        self.model.set("grid.acoustic_fov_length_mm", self.fov_len.value() or [])
        self.changed.emit()
        self._draw()

    _cursor_set_by_user = False

    # ── drawing ─────────────────────────────────────────────────────────
    def _request_draw(self) -> None:
        """Redraw the slices once control returns to the event loop: a burst of scroll
        steps or slider moves costs one draw. The readout does not depend on the crosshair,
        so it is left as it is."""
        if not hasattr(self, "_draw_timer"):
            self._draw_timer = QTimer(self, singleShot=True, interval=0)
            self._draw_timer.timeout.connect(lambda: self._draw(info=False))
        self._draw_timer.start()

    def _window(self) -> tuple[float, float]:
        """Display range of the T1 (1st to 99.5th percentile), once per image."""
        if getattr(self, "_window_for", None) is not self.t1:
            self._window_for = self.t1
            self._window_lohi = tuple(np.percentile(self.t1.canon[::4, ::4, ::4], [1, 99.5]))
        return self._window_lohi

    def _draw(self, info: bool = True) -> None:
        if not HAVE_MPL:
            return
        if hasattr(self, "_draw_timer"):
            self._draw_timer.stop()                 # this draw covers a pending request
        if self.t1 is None:
            self._setup_axes()
            self._slice_ratio = None
            self._apply_row_cap()
            self.canvas.draw_idle()
            self._update_info()
            return
        t1 = self.t1
        c = t1.ras_to_canon(self.cursor)[0]
        data = t1.canon
        vs = np.sqrt((t1.canon_affine[:3, :3] ** 2).sum(axis=0))
        lo, hi = self._window()
        trans, focus = self.positions()
        # One entry per transducer that has a drawable bowl: (is_selected, shapes, trans, focus).
        # Non-selected transducers are drawn dimmed, so a multi-transducer arrangement is visible
        # as a whole rather than only one bowl at a time (doc_multitransducer.md).
        tx_shapes = []
        for i in range(self.model.n_transducers):
            ti_trans, ti_focus = self._positions_for(i)
            ti_bowl = self._bowl_for(i)
            if ti_trans is not None and ti_focus is not None and ti_bowl is not None and \
                    np.linalg.norm(ti_focus - ti_trans) > 0:
                tx_shapes.append((i == self.ti, ti_bowl.outline(ti_trans, ti_focus - ti_trans),
                                  ti_trans, ti_focus, self.model.transducer_tooltip(i)))
        # Acoustic FOV box of every transducer, each resolved from its own geometry as
        # resolve_acoustic_fov_params.m does; coherent coupling crops to their union,
        # async coupling uses one per transducer (doc_partial-fov.md).
        fovs = []
        for i in range(self.model.n_transducers):
            ti_trans, ti_focus = (trans, focus) if i == self.ti else self._positions_for(i)
            ti_diam, ti_len, _, _ = self._resolved_fov(ti_trans, ti_focus, self._bowl_for(i))
            if ti_trans is not None and ti_focus is not None and ti_diam and ti_len \
                    and np.linalg.norm(ti_focus - ti_trans) > 0:
                fovs.append((i == self.ti, pl.fov_box(ti_trans, ti_focus, ti_diam, ti_len)))
        heur_pts, heur_prop = self._memo("heuristic", None, self._heuristic_candidate_points)
        grid_z = self._grid_z_line() if self.grid_z_action.isChecked() else None
        row_ratio = 0.0     # tallest per-column (height/width) need across the three views
        # The slice images, titles and legend are made once per T1/segmentation and then only
        # given new data; the markers below are redrawn each time (clearing and rebuilding the
        # axes cost about as much as rendering them).
        images = self._setup_axes(lo, hi, vs)
        alpha = self.alpha.value()
        for ax, (title, plane, (h, w)), (im, over) in zip(self.axes, VIEWS, images):
            k = int(np.clip(round(c[plane]), 0, data.shape[plane] - 1))
            sl = data[(slice(None),) * plane + (k,)]     # a view; remaining axes in order (h, w)
            row_ratio = max(row_ratio, (vs[w] / vs[h]) * (sl.shape[1] / sl.shape[0]))
            im.set_data(sl.T)
            if over is not None:
                over.set_visible(alpha > 0)
                if alpha > 0:
                    over.set_data(self._overlay(plane, k, sl.shape).T)
                    over.set_alpha(alpha / 100)
            for artist in list(ax.lines) + list(ax.collections):
                artist.remove()
            ax.axhline(c[w], color="#ffd54f", lw=0.5, alpha=0.6)
            ax.axvline(c[h], color="#ffd54f", lw=0.5, alpha=0.6)

            def proj(pts):
                p = t1.ras_to_canon(pts)
                return p[:, h], p[:, w]
            self.marker_tips.clear(ax)
            for selected, shapes, s_trans, s_focus, tip in tx_shapes:
                a = 1.0 if selected else 0.35
                rim_c, ax_c = ("#e53935", "#e53935") if selected else ("#9e9e9e", "#9e9e9e")
                x, y = proj(shapes["rim"])
                ax.plot(x, y, color=rim_c, lw=1, alpha=a)
                x, y = proj(np.vstack([s_trans, s_focus]))
                ax.plot(x, y, color=ax_c, lw=0.8, ls="--", alpha=a)
                x, y = proj(shapes["exit"][None])
                ax.plot(x, y, "s", color=rim_c, ms=3, alpha=a)
                x, y = proj(shapes["geom_focus"][None])
                ax.plot(x, y, "x", color="#bdbdbd", ms=5, alpha=a)
                for point, colour, marker in ((s_trans, rim_c, "o"), (s_focus, "#1e88e5" if selected else "#9e9e9e", "+")):
                    x, y = proj(point[None])
                    ax.plot(x, y, marker, color=colour, ms=9, mew=2, mfc="none" if marker == "o" else colour, alpha=a)
                    self.marker_tips.add(ax, x[0], y[0], tip if marker == "o" else f"{tip} · focus")
            for selected, fov in fovs:
                for edge in fov["edges"]:
                    x, y = proj(edge)
                    ax.plot(x, y, color="#43a047", lw=0.8, ls=(0, (3, 2)), alpha=1.0 if selected else 0.35)
            if grid_z is not None:
                x, y = proj(grid_z)
                ax.plot(x, y, color=style.active()["brand"], lw=1.0, ls=(0, (6, 2, 1, 2)))
            if len(heur_pts):
                x, y = proj(heur_pts)
                ax.scatter(x, y, c=heur_prop, cmap="cividis", vmin=0, vmax=0.3, s=8,
                          alpha=0.6, linewidths=0)
            ax.set_xlim(-0.5, sl.shape[0] - 0.5)
            ax.set_ylim(-0.5, sl.shape[1] - 0.5)
        # Row height needed per pixel of width, for _apply_row_cap: the slices keep a fixed
        # real-world (mm) aspect ratio, so with a plain Expanding size policy the canvas just
        # gets wider without the slices getting any bigger — the row's height is set once by
        # the surrounding layout and never revisited, and the aspect-locked axes just centre
        # themselves in the extra width as blank margin. _OVERHEAD_FRAC approximates the title
        # and axis-label chrome tight_layout reserves, as a fraction of the row height.
        self._slice_ratio = row_ratio / len(self.axes) / (1 - _OVERHEAD_FRAC - _LEGEND_FRAC)
        self.canvas.draw_idle()
        self._apply_row_cap()
        if info:
            self._update_info()

    def _setup_axes(self, lo: float = 0.0, hi: float = 1.0, vs=None) -> list:
        """(T1 image, tissue overlay image or None) per view, made when the T1 or the
        segmentation changed (or cleared when there is no T1)."""
        key = (id(self.t1), id(self.seg))
        if getattr(self, "_axes_for", None) == key:
            return self._images
        self._axes_for = key
        self._images = []
        for ax in self.axes:
            ax.clear()
            ax.set_xticks([])
            ax.set_yticks([])
        if self.t1 is None:
            return self._images
        shape = self.t1.canon.shape
        for ax, (title, plane, (h, w)) in zip(self.axes, VIEWS):
            blank = np.zeros((shape[w], shape[h]))
            im = ax.imshow(blank, cmap="gray", origin="lower", vmin=lo, vmax=hi,
                           aspect=vs[w] / vs[h], interpolation="nearest")
            over = None
            if self.seg is not None:
                over = ax.imshow(np.full((shape[w], shape[h]), np.nan), origin="lower",
                                 cmap=_TISSUE_CMAP, vmin=0, vmax=10, aspect=vs[w] / vs[h],
                                 interpolation="nearest")
            ax.set_title(f"{title}", fontsize=style.PLOT_PT["title"])
            ax.set_xlabel(f"{AXIS_LABELS[h][0]}  ←→  {AXIS_LABELS[h][1]}", fontsize=style.PLOT_PT["normal"])
            ax.set_ylabel(f"{AXIS_LABELS[w][0]}  ←→  {AXIS_LABELS[w][1]}", fontsize=style.PLOT_PT["normal"])
            ax.set_autoscale_on(False)
            self._images.append((im, over))
        self._draw_legend()
        return self._images

    def _apply_row_cap(self) -> None:
        if not HAVE_MPL:
            return
        if not self._slice_ratio:
            self.canvas.setMaximumHeight(16_777_215)           # Qt's QWIDGETSIZE_MAX: no cap
            return
        self.layout().activate()                # make self._body.geometry() reflect this resize
        self._content.layout().activate()
        avail = self._body.geometry().height() or self.canvas.height()
        ideal = self.canvas.width() * self._slice_ratio
        self.canvas.setMaximumHeight(max(50, round(min(ideal, avail))))

    def _overlay(self, plane: int, k: int, shape2d: tuple) -> np.ndarray:
        """Tissue labels on T1 slice k of the given plane (NaN outside the head), cached
        per slice so scrolling back and redrawing the same slices does not resample."""
        key = (id(self.t1), id(self.seg), plane, k)
        cache = self.__dict__.setdefault("_overlay_cache", {})
        if key not in cache:
            if len(cache) > 96:
                cache.clear()
            cache[key] = self._overlay_slice(plane, k, shape2d)
        return cache[key]

    def _overlay_slice(self, plane: int, k: int, shape2d: tuple) -> np.ndarray:
        if self.seg.shape == self.t1.shape and np.allclose(self.seg.affine, self.t1.affine, atol=1e-4):
            # same grid as the T1 (SimNIBS writes final_tissues on it): the slice itself
            labels = self.seg.canon[(slice(None),) * plane + (k,)].astype(float)
            labels[labels == 0] = np.nan
            return labels
        ii, jj = np.meshgrid(np.arange(shape2d[0]), np.arange(shape2d[1]), indexing="ij")
        vox = np.zeros((ii.size, 3))
        others = [a for a in range(3) if a != plane]
        vox[:, others[0]] = ii.ravel()
        vox[:, others[1]] = jj.ravel()
        vox[:, plane] = k
        labels = self.seg.sample(self.t1.canon_to_ras(vox)).astype(float)
        labels[labels == 0] = np.nan
        return labels.reshape(shape2d)

    def _draw_legend(self) -> None:
        """One line under the slices naming the markers (the Beam axis tab has its own key)."""
        for leg in list(self.fig.legends):
            leg.remove()
        marks = [
            Line2D([], [], ls="none", marker="+", color="#1e88e5", ms=8, mew=2, label="target"),
            Line2D([], [], ls="none", marker="o", mfc="none", color="#e53935", ms=7, mew=2,
                   label="bowl apex"),
            Line2D([], [], color="#e53935", lw=1, label="bowl rim"),
            Line2D([], [], ls="none", marker="s", color="#e53935", ms=3, label="exit plane"),
            Line2D([], [], color="#e53935", lw=0.8, ls="--", label="beam axis"),
            Line2D([], [], ls="none", marker="x", color="#9e9e9e", ms=5, label="geometric focus"),
            Line2D([], [], color="#43a047", lw=0.8, ls=(0, (3, 2)), label="acoustic FOV")]
        if self.grid_z_action.isChecked():
            marks.append(Line2D([], [], color=style.active()["brand"], lw=1.0, ls=(0, (6, 2, 1, 2)),
                                label="grid z"))
        self.fig.legend(handles=marks, loc="lower center", ncol=len(marks), fontsize=style.PLOT_PT["small"],
                        frameon=False, handlelength=1.6, columnspacing=1.2)
        self.fig.get_layout_engine().set(rect=(0, _LEGEND_FRAC * 0.6, 1, 1 - _LEGEND_FRAC * 0.6))

    # ── info panel ──────────────────────────────────────────────────────
    def _update_info(self) -> None:
        """The readout: target and depth, contact, skull on the beam, beam plane,
        positions, grid; the coordinate space and T1 as a footer."""
        trans, focus = self.positions()
        bowl = self._bowl()
        target, contact, skull = [], [], []
        if trans is not None and focus is not None and bowl is not None:
            axis = focus - trans
            dist = float(np.linalg.norm(axis))
            ep = dist - bowl.offset
            target.append(("Target depth", f"{ep:.1f} mm" + muted(" from the exit plane")))
            exp = self._expected_bowl_distance()
            if exp is not None:
                diff = dist - exp
                where = ("on target" if abs(diff) < 0.05 else
                         f"target {abs(diff):.1f} mm {'beyond' if diff > 0 else 'short of'} the focus")
                target.append(("Focal depth setting", f"{exp - bowl.offset:.1f} mm · "
                               + status(where, "ok" if abs(diff) <= DEPTH_TOL_MM else "bad")))
            dev = self.model.device(self.ti)
            if dev is not None and dev.min_foc is not None and dev.max_foc is not None:
                ok = dev.min_foc <= ep <= dev.max_foc
                target.append(("Steering range", f"{dev.min_foc:g}–{dev.max_foc:g} mm · "
                               + status("inside" if ok else "outside", "ok" if ok else "bad")))
            if self.seg is not None and dist > 0:
                scalp = pl.scalp_exit(self.seg, focus, trans - focus)
                if scalp is not None:
                    gap = float(np.dot(trans + axis / dist * bowl.offset - scalp, (trans - focus) / dist))
                    contact.append(("Exit plane → scalp", _mm_status(gap)
                                    + muted(" on the axis (outside)" if gap >= 0 else " (inside the head)")))
                rows = self._beam_qc_rows(trans, focus, bowl)
                contact += rows["contact"]
                skull = rows["skull"]
        positions = []
        for label, p in (("Target", focus), ("Bowl apex", trans)):
            if p is None:
                positions.append((label, muted("not set")))
                continue
            vox = pl.ras_to_prestus_voxel(self.t1, p) if self.t1 is not None else None
            tissue = self._tissue_at(p)
            positions.append((label, f"RAS ({p[0]:.1f}, {p[1]:.1f}, {p[2]:.1f}) mm<br>"
                              + muted(f"voxel {vox}" + (f" · {tissue}" if tissue else ""))))
        grid = []
        fov_diam, fov_len, fov_diam_auto, fov_len_auto = self._resolved_fov(trans, focus, bowl)
        if fov_diam and fov_len:
            auto = " (auto)" if fov_diam_auto and fov_len_auto else ""
            grid.append(("Acoustic FOV", f"⌀ {fov_diam:.0f} × {fov_len:.0f} mm" + muted(auto)))
        if self.model.n_transducers > 1:
            freq_issue = any(i.where() == "transducer.freq_hz" for i in checks.rule_issues(self.model))
            grid.append(("Transducers", f"{self.model.n_transducers} · grid.mode "
                         f"{self.model.get('grid.mode') or 'focal_axis'}"
                         + (" · " + status("frequencies differ", "bad") if freq_issue else "")))
        title = "Target" + (f" · transducer {self.ti + 1}" if self.model.n_transducers > 1 else "")
        footer = f"Coordinates: {self._space()}"
        if self.t1 is not None:
            footer += (f"<br>T1 {self.t1.path.name} "      # the image on its own line
                       f"({'×'.join(f'{s:.2g}' for s in self.t1.voxel_size)} mm)")
        self.info.setText(readout([(title, target), ("Contact", contact),
                                   ("Skull along the beam", skull),
                                   ("Positions", positions), ("Simulation grid", grid)], footer))

    def _beam_qc(self, trans, focus, bowl: pl.Bowl) -> dict[str, float]:
        """placement.beam_qc, cached per placement (and shared with qc_issues)."""
        key = (tuple(np.round(trans, 2)), tuple(np.round(focus, 2)),
               (bowl.curv_radius, bowl.aperture, bowl.dist_geom_ep), str(self.seg.path))
        cache = self.__dict__.setdefault("_qc_cache", {})
        if key not in cache:
            if len(cache) > 16:
                cache.clear()
            cache[key] = pl.beam_qc(self.seg, trans, focus, bowl)
        return cache[key]

    def qc_issues(self) -> list[Issue]:
        """The viewer's red flags as Issues-panel warnings, for every transducer: the target
        away from the focal depth setting, the exit-plane rim inside the scalp, and beam rays
        through an eye. Manual placement only (other modes are resolved again at run time),
        and only once the images are loaded (the viewer or the Run tab loaded them), so
        editing never waits for a NIfTI."""
        if self._loaded_for is None or self.model.get("placement.mode") not in ("manual", None, ""):
            return []
        self._load_images()                    # follows a changed subject or path
        if self.t1 is None:
            return []
        issues: list[Issue] = []
        source = self.position_source
        n = self.model.n_transducers
        try:
            for ti in range(n):
                trans, focus = self._positions_for(ti)
                bowl = self._bowl_for(ti)
                if trans is None or focus is None or bowl is None:
                    continue
                dist = float(np.linalg.norm(focus - trans))
                if dist <= 0:
                    continue
                who = f"Transducer {ti + 1}: " if n > 1 else ""
                exp = self._expected_bowl_distance(ti)
                if exp is not None and abs(dist - exp) > DEPTH_TOL_MM:
                    issues.append(Issue(
                        "warning", QC_PATH,
                        f"{who}the target is {abs(dist - exp):.1f} mm "
                        f"{'beyond' if dist > exp else 'short of'} the focus (target depth "
                        f"{dist - bowl.offset:.1f} mm, focal depth setting {exp - bowl.offset:.1f} mm). "
                        "'Set focal depth to target' or move the transducer.", ti))
                if self.seg is None:
                    continue
                q = self._beam_qc(trans, focus, bowl)
                # the exit plane's centre on the axis ('Exit plane vs scalp') and its rim
                # (beam QC): on a curved head either can be the one inside
                gaps = [q["rim_gap_min_mm"]]
                scalp = pl.scalp_exit(self.seg, focus, trans - focus)
                if scalp is not None:
                    gaps.append(float(np.dot(trans + (focus - trans) / dist * bowl.offset - scalp,
                                             (trans - focus) / dist)))
                worst = min((g for g in gaps if np.isfinite(g)), default=np.nan)
                if np.isfinite(worst) and worst < 0:
                    issues.append(Issue(
                        "warning", QC_PATH,
                        f"{who}the exit plane is up to {-worst:.1f} mm inside the scalp "
                        "(the transducer cuts into the head). 'Snap to scalp' or move it.", ti))
                if q["eye_in_beam"] > 0:
                    issues.append(Issue(
                        "warning", QC_PATH,
                        f"{who}{100 * q['eye_in_beam']:.0f}% of the beam rays pass an eye.", ti))
        finally:
            self.position_source = source
        return issues

    _grid_z_user: bool | None = None       # the View menu choice; None: on with 2+ transducers

    def _grid_z_toggled(self, on: bool) -> None:
        self._grid_z_user = on
        if self.t1 is not None and HAVE_MPL:
            self._draw_legend()
            self._request_draw()

    def _grid_z_line(self) -> np.ndarray | None:
        """Two RAS points 80 mm either side of the mean focus along the grid z axis that
        grid.mode gives (multitransducer.grid_z_axis, in the T1's RAS+ voxel grid as
        preproc_head.m), or None when it cannot be determined."""
        t1 = self.t1
        if t1 is None:
            return None
        def compute():
            to_canon = lambda r: t1.ras_to_canon(np.atleast_2d(r))[0]      # noqa: E731
            txs = mt.model_transducers(self.model, self._positions_for, to_canon)
            o = mt.orientation(self.model, txs)
            foci = [np.asarray(t["focus_pos"], float) for t in txs if t.get("focus_pos") is not None]
            if o["z"] is None or not foci:
                return None
            centre = np.mean(foci, axis=0)
            vs = np.sqrt((t1.canon_affine[:3, :3] ** 2).sum(axis=0))
            step = o["z"] * (80.0 / float(np.linalg.norm(o["z"] * vs)))   # 80 mm in voxels
            return t1.canon_to_ras(np.vstack([centre - step, centre + step]))
        return self._memo("grid_z", None, compute)

    def ensure_images(self) -> None:
        """Load the T1 and segmentation (if not yet), so qc_issues can check the placement."""
        self._load_images()

    def show_transducer(self, ti: int) -> None:
        if not 0 <= ti < self.model.n_transducers:
            return
        if ti < self.tx_combo.count():
            self.tx_combo.setCurrentIndex(ti)
        elif ti != self.ti:            # list not filled yet (Viewer not shown since the change)
            self.ti = ti
            self._cursor_set_by_user = False
            self.refresh()

    def _beam_qc_rows(self, trans, focus, bowl: pl.Bowl) -> dict[str, list]:
        """Contact and skull rows of placement.beam_qc (37 rays from the aperture to the
        focus; rim gap: exit-plane rim to scalp along the beam), cached per placement."""
        q = self._beam_qc(trans, focus, bowl)
        rim_min = q["rim_gap_min_mm"]
        contact = []
        if np.isfinite(rim_min):
            even = q["rim_gap_sd_mm"]
            contact.append(("Rim → scalp", f"{q['rim_gap_mean_mm']:.1f} ± {even:.1f} mm · min "
                            + _mm_status(rim_min)))
        inc = q["beam_skull_incidence_mean_deg"]
        skull = [("Bone on the beam", f"{q['beam_skull_path_mean_mm']:.1f} ± "
                  f"{q['beam_skull_path_sd_mm']:.1f} mm" + muted(f" · central ray {q['skull_path_mm']:.1f} mm")),
                 ("Incidence", (status(f"{inc:.0f}°", "warn" if inc > 20 else "ok") if np.isfinite(inc)
                                else muted("n/a")) + muted(" beam mean"))]
        if q["eye_in_beam"] > 0:
            skull.append(("Eye", status(f"{100 * q['eye_in_beam']:.0f}% of beam rays", "bad")))
        return {"contact": contact, "skull": skull}

    def _tissue_at(self, ras) -> str:
        if self.seg is None:
            return ""
        label = int(self.seg.sample(np.asarray(ras)[None])[0])
        return TISSUE.get(label, ("outside head" if label == 0 else f"label {label}", ""))[0]

    # ── interaction ─────────────────────────────────────────────────────
    def _select_transducer(self, index: int) -> None:
        if index >= 0:
            self.ti = index
            self._cursor_set_by_user = False
            self.refresh()
            self.transducer_selected.emit(self.ti)

    def _add_transducer(self) -> None:
        self._record("Add transducer")
        self.ti = self.model.add_transducer()
        self._cursor_set_by_user = False
        self.changed.emit()
        self.refresh()
        self.transducer_selected.emit(self.ti)

    def _remove_transducer(self) -> None:
        self._record("Remove transducer")
        self.model.remove_transducer(self.ti)
        self.ti = max(0, self.ti - 1)
        self._cursor_set_by_user = False
        self.changed.emit()
        self.refresh()
        self.transducer_selected.emit(self.ti)

    def _event_ras(self, event) -> np.ndarray | None:
        if self.t1 is None or event.inaxes not in self.axes or event.xdata is None:
            return None
        _, plane, (h, w) = VIEWS[self.axes.index(event.inaxes)]
        c = self.t1.ras_to_canon(self.cursor)[0]
        c[h], c[w] = event.xdata, event.ydata
        return self.t1.canon_to_ras(c[None])[0]

    def _clicked(self, event) -> None:
        ras = self._event_ras(event)
        if ras is None:
            return
        self.cursor = ras
        self._cursor_set_by_user = True
        if event.button == 1 and self._editable() and not self.rb_browse.isChecked():
            self._write("focus_pos" if self.rb_focus.isChecked() else "trans_pos", ras)
            self.refresh()
        else:
            self._request_draw()

    def _scrolled(self, event) -> None:
        if self.t1 is None or event.inaxes not in self.axes:
            return
        _, plane, _ = VIEWS[self.axes.index(event.inaxes)]
        c = self.t1.ras_to_canon(self.cursor)[0]
        c[plane] += 1 if event.button == "up" else -1
        self.cursor = self.t1.canon_to_ras(c[None])[0]
        self._cursor_set_by_user = True
        self._request_draw()

    def _centre(self) -> None:
        _, focus = self.positions()
        if focus is not None:
            self.cursor = focus.copy()
            self._request_draw()

    def _entry_ras(self, transformation_type: str = "nonl") -> np.ndarray | None:
        """RAS mm for the typed X/Y/Z in the chosen entry space, or None (after a message) if
        it can't be resolved yet. For MNI, `transformation_type` picks the SimNIBS transform:
        'nonl' (nonlinear warp, most accurate in-brain) or '12dof'/'6dof' (affine, well-defined
        everywhere including outside the brain) -- see doc_placement_mni.md."""
        xyz = np.array([self.cx.value(), self.cy.value(), self.cz.value()])
        space = self.coord_space.currentData()
        if space == "ras":
            return xyz
        if space == "voxel":
            if self.t1 is None:
                QMessageBox.information(self, "Coordinates", "No T1 image loaded.")
                return None
            return pl.prestus_voxel_to_ras(self.t1, xyz)
        sid = self.model.get("subject_id")
        seg_root = self.model.get("path.seg") or self.model.get("path.anat")
        simnibs = self.model.get("startup.simnibs_bin_path")
        if not isinstance(sid, (int, float)) or is_empty(seg_root):
            QMessageBox.information(self, "MNI coordinates", "Set subject_id and path.anat "
                                    "(or path.seg) first.")
            return None
        m2m = paths.m2m_dir(seg_root, int(sid))
        if not m2m.is_dir():
            QMessageBox.information(self, "MNI coordinates", f"{m2m.name} not found under "
                                    f"{seg_root}. SimNIBS (charm) must be run for this subject "
                                    "before MNI coordinates can be converted here -- both the "
                                    "nonlinear and affine transform come from that output; "
                                    "setting startup.simnibs_bin_path alone is not enough.")
            return None
        if is_empty(simnibs):
            QMessageBox.information(self, "MNI coordinates", "Set startup.simnibs_bin_path (the "
                                    "SimNIBS bin folder) to run mni2subject_coords.")
            return None
        ld_lib = self.model.get("hpc.ld_library_path")
        try:
            return pl.mni_to_subject(xyz, m2m, simnibs, ld_lib, transformation_type)
        except pl.SimnibsCoordError as exc:
            if transformation_type != "nonl":
                QMessageBox.warning(self, "MNI coordinates", str(exc))
                return None
            try:                                            # nonlinear warp missing/failed: fall
                ras = pl.mni_to_subject(xyz, m2m, simnibs, ld_lib, "12dof")          # back to the
            except pl.SimnibsCoordError:                                            # affine, which
                QMessageBox.warning(self, "MNI coordinates", str(exc))              # SimNIBS always
                return None                                                         # produces
            QMessageBox.warning(self, "MNI coordinates", f"The nonlinear transform failed ({exc}); "
                                "used the affine (12dof) transform instead, which is less accurate "
                                "for an in-brain point.")
            return ras

    def _transform_for(self, key: str) -> str:
        """PRESTUS's own convention for MNI conversion (doc_placement_mni.md): the nonlinear
        warp for an in-brain focus, the affine for a transducer/scalp point (where the nonlinear
        field is unreliable/extrapolated, and which the affine covers without needing SimNIBS's
        nonlinear registration to have succeeded)."""
        return "12dof" if key == "trans_pos" else "nonl"

    def _goto_coords(self) -> None:
        key = "trans_pos" if self.rb_trans.isChecked() else "focus_pos"
        ras = self._entry_ras(self._transform_for(key))
        if ras is None:
            return
        self.cursor = ras
        self._cursor_set_by_user = True
        self._request_draw()

    def _set_coords(self, key: str) -> None:
        if not self._editable():
            QMessageBox.information(self, "Placement", "Switch to manual placement with "
                                    "coordinate space 'voxel' or 'ras' to place the transducer here.")
            return
        ras = self._entry_ras(self._transform_for(key))
        if ras is None:
            return
        self.cursor = ras
        self._cursor_set_by_user = True
        self._write(key, ras)
        self.refresh()

    def _direction(self) -> np.ndarray | None:
        trans, focus = self.positions()
        if focus is None:
            return None
        if trans is not None and np.linalg.norm(trans - focus) > 0:
            return trans - focus
        if self.direction is not None:
            return self.direction
        if self.seg is not None:           # radial from the head's centre
            idx = np.argwhere(self.seg.data[::4, ::4, ::4] > 0) * 4
            if len(idx):
                return focus - self.seg.to_ras(idx.mean(axis=0))[0]
        return None

    def _need(self, what: str) -> bool:
        self._finish_load()
        if not self._editable():
            QMessageBox.information(self, "Placement", "Switch to manual placement with coordinate "
                                    "space 'voxel' or 'ras' to place the transducer here.")
            return False
        if what == "focus" and self.positions()[1] is None:
            QMessageBox.information(self, "Placement", "Set the focus first.")
            return False
        return True

    def _place_on_axis(self) -> None:
        if not self._need("focus"):
            return
        d, dist = self._direction(), self._expected_bowl_distance()
        if d is None or dist is None:
            QMessageBox.information(self, "Placement", "Needs a beam direction (transducer, trajectory, "
                                    "or segmentation) and a focal depth setting (focal_distance_ep or _bowl).")
            return
        self._write("trans_pos", pl.place_on_axis(self.positions()[1], d, dist))
        self.refresh()

    def _snap(self) -> None:
        if not self._need("focus"):
            return
        d, bowl = self._direction(), self._bowl()
        if self.seg is None or d is None or bowl is None:
            QMessageBox.information(self, "Placement", "Needs the segmentation, a beam direction, "
                                    "and the transducer geometry.")
            return
        apex = pl.snap_to_scalp(self.seg, self.positions()[1], d, bowl, self.gap.value())
        if apex is None:
            QMessageBox.warning(self, "Placement", "No scalp found along the beam axis.")
            return
        self._write("trans_pos", apex)
        self.refresh()

    def _match_focal_distance(self) -> None:
        trans, focus = self.positions()
        bowl = self._bowl()
        if trans is None or focus is None or bowl is None:
            QMessageBox.information(self, "Placement", "Needs transducer, focus, and transducer geometry.")
            return
        ep = round(float(np.linalg.norm(focus - trans)) - bowl.offset, 2)
        self._record("Set focal depth to target")
        self.model.set(f"{TX}.focal_distance_ep", ep, self.ti)
        self.model.reset(f"{TX}.focal_distance_bowl", self.ti)      # derived from ep by PRESTUS
        self.changed.emit()
        self.refresh()

    def _import_trajectory(self) -> None:
        self._finish_load()
        path, _ = QFileDialog.getOpenFileName(
            self, "Import trajectory", "",
            "Trajectories (*.txt *.tfm *.kdp *.kps *.h5 *.xml);;Brainsight / 3D Slicer (*.txt *.tfm);;"
            "k-Plan plan (*.kdp *.kps *.h5);;Localite export (*.xml);;All files (*)")
        if not path:
            return
        if Path(path).suffix.lower() == ".xml":
            self._import_localite(path)
            return
        if not self._need("any"):
            return
        try:
            traj = pl.read_trajectory(path, t1=self.t1)
        except Exception as exc:
            QMessageBox.warning(self, "Import trajectory", str(exc))
            return
        if traj.transducer is not None:
            self._import_pose(traj)
            return
        self._write("focus_pos", traj.target)
        self.direction = traj.towards_transducer
        dist = self._expected_bowl_distance()
        if dist is not None:
            self._write("trans_pos", pl.place_on_axis(traj.target, traj.towards_transducer, dist))
        self.cursor = traj.target.copy()
        self._cursor_set_by_user = True
        self.refresh()
        QMessageBox.information(self, "Import trajectory",
                                f"{traj.source} target '{traj.name}' imported"
                                + ("" if dist is not None else
                                   ". Set the focal depth setting (focal_distance_ep) to place the transducer on its axis") + ".")

    def _import_localite(self, path: str) -> None:
        """A Localite XML export: transducer ti uses one of its series at run time
        (placement.mode localite, placement.localite.file), where PRESTUS averages the
        series' triggers and converts the pose. MATLAB lists the series and their positions
        with that same code (nexus_neuronav 'file_series')."""
        bridge = self.bridge_getter() if self.bridge_getter is not None else None
        cfg = self.saved_config_getter() if self.saved_config_getter is not None else None
        if bridge is None or cfg is None:
            if bridge is None:
                QMessageBox.information(self, "Import trajectory",
                                        "Reading a Localite export needs MATLAB (Settings ▸ MATLAB).")
            return
        ti = self.ti
        self.setCursor(Qt.CursorShape.BusyCursor)

        def failed(err: BridgeError) -> None:
            self.unsetCursor()
            QMessageBox.warning(self, "Import trajectory", err.message)

        bridge.call("nexus_neuronav", {"action": "file_series", "config": cfg, "file": path},
                    on_done=lambda out: self._localite_series_read(path, ti, out), on_error=failed)

    def _localite_series_read(self, path: str, ti: int, out: dict) -> None:
        self.unsetCursor()
        series = [dict(s) for s in as_list(out.get("series")) if isinstance(s, dict)]
        usable = lc.usable_series(series)
        if not usable:
            QMessageBox.warning(self, "Import trajectory",
                                f"{Path(path).name} has no recorded, tracked transducer pose "
                                f"({len(series)} series: planned points or untracked).")
            return
        chosen = usable[0]
        if len(usable) > 1:
            labels = [lc.series_label(s) for s in usable]
            what = f"transducer {ti + 1}" if self.model.n_transducers > 1 else "the transducer"
            label, ok = QInputDialog.getItem(
                self, "Import trajectory",
                f"{Path(path).name} has {len(usable)} series. Place {what} with:", labels, 0, False)
            if not ok:
                return
            chosen = usable[labels.index(label)]
        self._record("Import Localite pose")
        changes = lc.file_import_settings(self.model, path, chosen, series, ti)
        for key, value in changes:
            if value is None:
                self.model.reset(key)
            else:
                self.model.set(key, value)
        for key in ("trans_pos", "focus_pos"):         # the run resolves them from the file
            self.model.reset(f"{TX}.{key}", ti)
        trans, focus = (chosen.get(k) for k in ("trans_pos", "focus_pos"))
        trans = as_floats(trans) if trans is not None and len(as_floats(trans)) == 3 else None
        focus = as_floats(focus) if focus is not None and len(as_floats(focus)) == 3 else None
        self._localite_file_pos[(path, ti)] = (
            repr(self.model.get("placement.localite.target")),
            repr(self.model.get("placement.localite.position")), trans, focus)
        self._forget()
        self.changed.emit()
        self.refresh()
        tgt = dict(changes).get("placement.localite.target")
        descr = chosen.get("description") or ""
        by_name = bool(descr) and (tgt[ti] if isinstance(tgt, list) else tgt) == descr
        pick = f"target '{descr}'" if by_name else f"series {int(chosen['index'])}"
        lines = [f"Placement mode set to Localite: {pick} of {Path(path).name}, "
                 f"the mean of {int(chosen.get('n_markers') or 0)} marker(s), resolved again at run time."]
        if trans is None:
            lines.append(chosen.get("error") or out.get("note") or "The position shows once the run resolves it.")
        QMessageBox.information(self, "Import trajectory", "\n".join(lines))

    def _import_pose(self, traj: "pl.Trajectory") -> None:
        """A trajectory that gives the bowl apex (k-Plan): the transducer goes exactly
        there, and the focal depth setting is matched to the target with this
        transducer's own geometry (k-Plan's model may put its exit plane elsewhere)."""
        self._write("trans_pos", traj.transducer)
        self.direction = traj.towards_transducer
        bowl = self._bowl()
        target = traj.target
        if target is None:
            dist = self._expected_bowl_distance()
            if dist is not None:
                target = pl.place_on_axis(traj.transducer, -traj.towards_transducer, dist)
        if target is not None:
            self._write("focus_pos", target)
            self.cursor = np.asarray(target, float).copy()
            self._cursor_set_by_user = True
        ep = None
        if traj.target is not None and bowl is not None:
            ep = round(float(np.linalg.norm(traj.target - traj.transducer)) - bowl.offset, 2)
            self._record("Set focal depth to target")
            self.model.set(f"{TX}.focal_distance_ep", ep, self.ti)
            self.model.reset(f"{TX}.focal_distance_bowl", self.ti)
            self.changed.emit()
        self.refresh()
        lines = [f"{traj.source} plan '{traj.name}' imported: transducer at the plan's bowl apex"
                 + (", focus at its target" if traj.target is not None else
                    "; the file has no target, so the focus is on the beam axis at the current focal depth"
                    if target is not None else "; the file has no target: place the focus on the beam axis") + "."]
        if ep is not None:
            lines.append(f"Focal depth set to {ep:.2f} mm from this transducer's exit plane.")
        elif traj.target is not None:
            lines.append("Choose the transducer to set the focal depth from its geometry.")
        lines += traj.notes
        QMessageBox.information(self, "Import trajectory", "\n".join(lines))

    def _use_manual(self) -> None:
        self._record("Use manual placement")
        self.model.set("placement.mode", "manual")
        self.changed.emit()
        self.refresh()


if HAVE_MPL:
    _TISSUE_CMAP = ListedColormap(["#000000"] + [TISSUE[i][1] for i in range(1, 11)]
                                  ).with_extremes(bad=(0, 0, 0, 0))
