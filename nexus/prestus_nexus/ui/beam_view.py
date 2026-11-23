"""Beam axis tab: a section through the beam axis of the selected transducer with the simulation
medium, the transducer drawn like the transducer schematic, and the free-water profile on the
same depth axis.

The tab reads the images, positions and geometry from the placement viewer (PlacementPage), so
the transducer selection, image loading and position sources stay in one place."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtWidgets import (QButtonGroup, QFrame, QHBoxLayout, QLabel, QPushButton, QToolButton,
                             QVBoxLayout, QWidget)

from .. import placement as pl
from ..bridge import as_columns, as_floats
from ..model import TX, is_empty
from ..outputs import (axial_profile_paths, dir_output, placement_record_path, read_axial_profile,
                       read_placement_record)
from . import style
from .empty_state import EmptyStack
from .errors import log_current
from .export import subject_stem
from .readout import muted, readout, status

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

# Simulation media (medium_properties names)
MEDIUM_COLOURS = {"water": "#cfe6f5", "brain": "#a6a6a6", "skin": "#e8a07a", "skull": "#f5deb3",
                  "skull_cortical": "#fff3c4", "skull_trabecular": "#d9b779"}
_EXTRA_COLOURS = ["#8e7cc3", "#76a5af", "#c27ba0", "#6aa84f", "#e69138"]
BEAM_STEP_MM = 0.5             # sampling of the beam-plane section
BEAM_TURN_DEG = 15.0           # one turn of the plane about the axis (scroll or buttons)
GAP_MODES = ("heuristic", "mni", "plantus", "tus_entry")    # modes that apply placement.skin_gap_mm
PROFILE_SESSION = "placement_profile"     # nexus_calibration session of this tab
_PROFILE_FRAC = 0.42           # profile strip height relative to the beam plane
_PROFILE_GAP = 0.1             # gap between beam plane and profile strip (beam-plane heights)
FREE_WATER = style.SERIES[0]   # categorical (style.SERIES); the simulated profile is drawn in
                               # the theme's ink (style.active()["ink"]), so it shows on dark too


def profile_from_payload(out) -> dict | None:
    """{'z_ep', 'I' (normalised to its peak), 'focus_ep', 'model'} from a nexus_calibration
    'view' payload at one depth: the library model, else the geometric one."""
    res = (out or {}).get("results") or {} if isinstance(out, dict) else {}
    for kind in ("library", "geometric"):
        r = res.get(kind)
        if not isinstance(r, dict):
            continue
        z, cols = as_floats(r.get("z")), as_columns(r.get("Imod"))
        if cols and len(cols) == len(z) and all(len(c) == 1 for c in cols):
            cols = [[c[0] for c in cols]]        # one depth sent as a flat column
        if not z or not cols or len(cols[0]) != len(z):
            continue
        i = np.asarray(cols[0], float)
        peak = np.nanmax(i) if np.isfinite(i).any() else 0
        if not peak > 0:
            continue
        focus = as_floats(r.get("focus"))
        return {"z_ep": z, "I": i / peak, "focus_ep": focus[0] if focus else float("nan"),
                "model": kind}
    return None


def _mm_status(v: float) -> str:
    """Signed distance [mm], ✗ when negative, 'n/a' when unknown."""
    if not np.isfinite(v):
        return muted("n/a")
    return status(f"{v:+.1f} mm", "bad" if v < 0 else "ok")


def _finite(a) -> bool:
    return a is not None and bool(np.isfinite(np.asarray(a, float)).any())


# The run's outputs live on the project volume, where even a stat can take seconds: the tab
# reads them here, off the GUI thread, and draws with what it has meanwhile.
_IO = ThreadPoolExecutor(max_workers=1, thread_name_prefix="nexus-beam-io")


def _read_run_outputs(medium_path: Path | None, profile_paths, record_path: Path | None,
                      medium_key, medium):
    """(medium key, medium volume with its voxels read, axial profile, placement record) of
    the last run. The volume is reused while its file is unchanged (same key)."""
    key = vol = None
    if medium_path is not None:
        try:
            key = (str(medium_path), medium_path.stat().st_mtime)
        except OSError:
            key = None
    if key is not None:
        if key == medium_key and medium is not None:
            vol = medium
        else:
            try:
                vol = pl.Volume(medium_path)
                vol.data                                # read now, not on the GUI thread
            except Exception:  # best-effort: the beam is drawn without the medium (nibabel raises many types)
                log_current(f"beam view: cannot read {medium_path}")
                vol = None
    profile = read_axial_profile(*profile_paths) if profile_paths is not None else None
    record = read_placement_record(record_path) if record_path is not None else None
    return key, vol, profile, record


class BeamPage(QWidget):
    navigate = pyqtSignal(str)        # config key to show (MainWindow.goto), from the empty state

    def __init__(self, viewer, parent: QWidget | None = None):
        super().__init__(parent)
        self.v = viewer                             # PlacementPage: images, positions, geometry
        self.model = viewer.model
        self.beam_angle = 0.0
        self.scale = "scaled"                       # profile scale: 'scaled' | 'absolute'
        self.medium: pl.Volume | None = None        # last run's medium_masks, if written
        self._medium_key: tuple | None = None
        self._beam_cache: dict = {}
        self.bridge_getter = None                   # MainWindow.bridge: the shared MATLAB session
        self._profiles: dict = {}                   # MATLAB prediction per (device, depth)
        self._profile_pending = None
        self._profile_errors: dict = {}
        self._run_profile: dict | None = None       # last run's axial_profile CSV (on refresh)
        self._outputs_read = False                  # the run's outputs were looked for once
        self._placement_record: dict | None = None  # the run's resolved positions (JSON)
        self._io = None                             # pending _read_run_outputs
        self._io_again = False
        self._io_timer = QTimer(self, interval=50)
        self._io_timer.timeout.connect(self._poll_outputs)
        self._present: list = []                    # media in the section (key swatches)
        self._gap_line = False
        self._rows: dict[str, list] = {}

        v = QVBoxLayout(self)                        # the sub-tab label names the page
        content = QWidget()
        cv = QVBoxLayout(content)
        cv.setContentsMargins(0, 0, 0, 0)
        self.stack = EmptyStack(content)
        v.addWidget(self.stack, 1)

        bar = QHBoxLayout()
        bar.setSpacing(6)
        self.tx_label = QLabel()
        self.tx_label.setObjectName("CardSub")
        bar.addWidget(self.tx_label)
        bar.addSpacing(10)
        bar.addWidget(QLabel("Plane"))
        self.turn_left = QPushButton("⟲")
        self.turn_left.setToolTip(f"Turn the section {BEAM_TURN_DEG:g}° about the beam axis")
        self.turn_left.clicked.connect(lambda: self._turn(-BEAM_TURN_DEG))
        self.angle_label = QLabel()
        self.turn_right = QPushButton("⟳")
        self.turn_right.setToolTip(f"Turn the section {BEAM_TURN_DEG:g}° the other way")
        self.turn_right.clicked.connect(lambda: self._turn(BEAM_TURN_DEG))
        for w in (self.turn_left, self.angle_label, self.turn_right):
            if isinstance(w, QPushButton):
                w.setProperty("compact", True)
                w.setFixedWidth(30)
            bar.addWidget(w)
        bar.addStretch(1)
        bar.addWidget(QLabel("Profile"))
        self.scale_box = QFrame()
        self.scale_box.setObjectName("Segmented")
        seg = QHBoxLayout(self.scale_box)
        seg.setContentsMargins(2, 2, 2, 2)
        seg.setSpacing(2)
        self.rb_scaled = QToolButton(text="&Scaled")
        self.rb_scaled.setToolTip("Each profile scaled to its own maximum: compare where they peak")
        self.rb_absolute = QToolButton(text="&Absolute")
        self.rb_absolute.setToolTip("Intensity in W/cm²: compare how much reaches the target")
        self.scale_group = QButtonGroup(self)
        for b in (self.rb_scaled, self.rb_absolute):
            b.setCheckable(True)
            self.scale_group.addButton(b)
            seg.addWidget(b)
        self.rb_scaled.setChecked(True)
        self.scale_group.buttonClicked.connect(self._scale_changed)
        bar.addWidget(self.scale_box)
        self.compute = QPushButton("Compute free-water profile")
        self.compute.setToolTip("Predict the device's free-water profile at the focal depth "
                                "setting in MATLAB, as the Calibration page does (starts MATLAB)")
        self.compute.clicked.connect(lambda: self._request_profile(start=True))
        bar.addWidget(self.compute)
        cv.addLayout(bar)

        body = QHBoxLayout()
        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 6))
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="beam",
                                            filename_hint=subject_stem(self.model, "beam"))
            gs = self.fig.add_gridspec(2, 2, width_ratios=[3.2, 1],
                                       height_ratios=[1, _PROFILE_GAP + _PROFILE_FRAC + 0.08],
                                       left=0.07, right=0.99, top=0.95, bottom=0.08, wspace=0.08,
                                       hspace=0)
            self._beam_cell = gs[0, 0]
            self.beam_ax = self.fig.add_subplot(self._beam_cell)
            self.key_ax = self.fig.add_subplot(gs[:, 1])
            self._new_profile_axes()
            self.canvas.mpl_connect("button_press_event", self._clicked)
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            body.addWidget(self.canvas, 3)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.info.setMinimumWidth(280)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.info.setToolTip(
            "Exit plane → scalp: along the beam axis, from the SimNIBS labels.\n"
            "Exit plane → medium: to the first tissue of the simulation medium in this section.\n"
            "Free water: the last run's on-axis profile, or the library model's prediction.")
        body.addWidget(self.info, 1)
        cv.addLayout(body, 1)
        hint = QLabel("A section through the beam axis of the transducer selected on the Viewer. "
                      "Scroll on it (or use ⟲ ⟳) to turn it about the axis; a click moves the "
                      "Viewer's crosshair.")
        hint.setObjectName("CardSub")
        hint.setWordWrap(True)
        cv.addWidget(hint)

    # ── state ───────────────────────────────────────────────────────────
    def _bridge(self):
        return self.bridge_getter() if self.bridge_getter is not None else None

    def _media(self) -> tuple[list[str], dict]:
        """(medium_properties names in config order, layers): a medium_masks value v is
        medium_names[v - 1] (preproc_medium_mask.m)."""
        eff = self.model.effective()
        props, layers = eff.get("medium_properties"), eff.get("layers")
        return (list(props) if isinstance(props, dict) else [],
                layers if isinstance(layers, dict) else {})

    def _layered(self) -> bool:
        return (self.model.get("simulation.medium") or "layered") == "layered"

    def _medium_path(self) -> Path | None:
        """Where the last layered run writes its medium_masks in T1 space
        (nifti_medium_masks.m); not checked on disk."""
        from ..study import affix_path
        base, sid = dir_output(self.model), self.model.get("subject_id")
        if base is None or not self._layered():
            return None
        return affix_path(lambda a: base / "nii" / f"sub-{int(sid):03d}_layered_T1w{a}_medium_masks.nii.gz",
                          self.model.get("io.output_affix") or "")

    def _read_outputs(self) -> None:
        """Look for the run's medium and profile: at once while the window is not shown
        (offscreen), else in the background and redraw when they change."""
        args = (self._medium_path(),
                axial_profile_paths(self.model) if self.model.n_transducers == 1 else None,
                placement_record_path(self.model))
        if not self.window().isVisible():
            self._take_outputs(_read_run_outputs(*args, self._medium_key, self.medium))
            return
        if self._io is not None:                    # read again once this one is in
            self._io_again = True
            return
        self._io = _IO.submit(_read_run_outputs, *args, self._medium_key, self.medium)
        self._io_timer.start()

    def _poll_outputs(self) -> None:
        if self._io is None or not self._io.done():
            return
        self._io_timer.stop()
        future, self._io = self._io, None
        try:
            changed = self._take_outputs(future.result())
        except Exception:  # best-effort: the view keeps what it shows; logged
            log_current("beam view: cannot read the run outputs")
            changed = not self._outputs_read
            self._outputs_read = True
        if self._io_again:
            self._io_again = False
            self._read_outputs()
        if changed and self.isVisible():
            self._draw()

    def _take_outputs(self, result) -> bool:
        """Keep what _read_run_outputs found; True when it changes what is drawn."""
        key, vol, profile, record = result
        stamp = lambda p: None if p is None else (p.get("path"), p.get("mtime"))   # noqa: E731
        changed = (not self._outputs_read or key != self._medium_key
                   or stamp(profile) != stamp(self._run_profile))
        if key != self._medium_key:
            self._beam_cache.clear()
        self._medium_key, self.medium, self._run_profile = key, vol, profile
        changed |= (record or {}).get("created") != (self._placement_record or {}).get("created")
        self._placement_record = record
        self._outputs_read = True
        return changed

    def _profile_key(self) -> tuple | None:
        """(serial, driving system, library folder, focal depth [mm from the exit plane]) of the
        selected transducer, or None without a library device or a focal depth setting."""
        ti = self.v.ti
        g = lambda k: self.model.get(f"{TX}.{k}", ti)        # noqa: E731
        serial, bowl, exp = g("serial"), self.v._bowl(), self.v._expected_bowl_distance()
        if is_empty(serial) or bowl is None or exp is None:
            return None
        ds, lib = g("combo.ds_serial"), g("library_path")
        return (str(serial), "" if is_empty(ds) else str(ds), "" if is_empty(lib) else str(lib),
                round(exp - bowl.offset, 1))

    def _run_usable(self) -> bool:
        run = self._run_profile
        return run is not None and (_finite(run.get("oneil")) or _finite(run.get("intensity")))

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        self.v._load_images()
        n = self.model.n_transducers
        self.tx_label.setText(f"Transducer {self.v.ti + 1} of {n}" if n > 1 else "")
        self._read_outputs()
        self._apply_style()
        self._draw()

    def _apply_style(self) -> None:
        c = style.active()
        self.scale_box.setStyleSheet(
            f"QFrame#Segmented {{ background: {c['surface']}; border: 1px solid {c['line2']}; "
            f"border-radius: 6px; }}"
            f"QFrame#Segmented QToolButton {{ border: none; border-radius: 4px; padding: 3px 10px; "
            f"color: {c['ink2']}; }}"
            f"QFrame#Segmented QToolButton:hover {{ background: {c['surface3']}; }}"
            f"QFrame#Segmented QToolButton:checked {{ background: {c['brand']}; color: white; "
            f"font-weight: 600; }}"
            f"QFrame#Segmented QToolButton:disabled {{ color: {c['ink4']}; }}")

    def _show_empty(self) -> bool:
        """An empty state for what is missing before a section can be drawn; False if none."""
        v, m = self.v, self.model
        if not HAVE_MPL:
            self.stack.show_empty("The beam view needs matplotlib and nibabel",
                                  "Install them in the Nexus's Python environment and restart.")
            return True
        if v.t1 is None:
            self.stack.show_empty("No T1 image yet", "The beam section is cut from the subject's "
                                  "T1 image and segmentation. Set the anatomy folder first.",
                                  "Set T1 image…", lambda: self.navigate.emit("path.anat"))
            return True
        trans, focus = v.positions()
        if trans is None or focus is None or np.linalg.norm(focus - trans) == 0:
            if v._space() == "mni":
                text = ("The positions are in MNI coordinates, which this view cannot map to the "
                        "subject. Use 'ras' or 'voxel' coordinates to see the beam.")
            elif m.get("placement.mode") not in ("manual", None, ""):
                text = ("No resolved position for this placement mode yet. Run the placement "
                        "stage once (only modules.run_transducer_placement), or place manually.")
            else:
                text = "Set the focus and the transducer (bowl apex) on the Viewer tab first."
            self.stack.show_empty("No placement to show", text)
            return True
        if v._bowl() is None:
            self.stack.show_empty("No transducer geometry", "Choose a device, or enter the "
                                  "radius of curvature and element diameters.",
                                  "Choose transducer…", lambda: self.navigate.emit("transducer.serial"))
            return True
        return False

    # ── drawing ─────────────────────────────────────────────────────────
    def _new_profile_axes(self) -> None:
        """The profile strip: an inset of the beam plane, so it follows the beam plane's 1:1
        box (whatever width the aspect leaves it) and shares its depth axis."""
        self.beam_ax.set_anchor("N")
        self.prof_ax = self.beam_ax.inset_axes(
            [0, -_PROFILE_GAP - _PROFILE_FRAC, 1, _PROFILE_FRAC], sharex=self.beam_ax)

    def _draw(self) -> None:
        if self._show_empty():
            self.info.setText("")
            return
        self.stack.show_content()
        self.beam_ax.clear()                    # also drops the profile inset: make it again
        self._new_profile_axes()
        self.key_ax.clear()
        self.key_ax.set_axis_off()
        self._rows = {"contact": [], "plane": [], "water": [], "missing": []}
        self.angle_label.setText(f"{self.beam_angle:.0f}°")
        trans, focus = self.v.positions()
        self._draw_beam(trans, focus)
        self._draw_profile(trans, focus)
        self._draw_key()
        self.canvas.draw_idle()
        self._update_info()

    def _beam_section(self, trans, focus, depth_max: float, half: float):
        """(depth, offsets, T1 section, medium section, source) for the current plane; medium
        values index _media()[0] (1-based; 0 = outside the simulation grid)."""
        names, layers = self._media()
        seg = self.v.seg
        source = "run" if self.medium is not None else \
            ("labels" if seg is not None and self._layered() else "")
        key = (tuple(np.round(trans, 2)), tuple(np.round(focus, 2)), round(self.beam_angle, 3),
               round(depth_max, 1), round(half, 1), source, self._medium_key,
               str(seg.path) if seg is not None else "", str(self.v.t1.path), tuple(names),
               repr(sorted(layers.items())))
        if key not in self._beam_cache:
            if len(self._beam_cache) > 8:
                self._beam_cache.clear()
            axis, lat = pl.beam_plane_basis(trans, focus, self.beam_angle)
            depth = np.arange(-10.0, depth_max + BEAM_STEP_MM / 2, BEAM_STEP_MM)
            offsets = np.arange(-half, half + BEAM_STEP_MM / 2, BEAM_STEP_MM)
            pts = pl.beam_plane_points(trans, axis, lat, depth, offsets)
            flat = pts.reshape(-1, 3)
            t1 = self.v.t1.sample(flat).astype(float).reshape(pts.shape[:2])
            med = None
            if source == "run":
                med = self.medium.sample(flat).astype(int).reshape(pts.shape[:2])
            elif source == "labels":
                med = pl.segmentation_to_medium(seg.sample(flat).reshape(pts.shape[:2]), layers, names)
            self._beam_cache[key] = (depth, offsets, t1, med, source)
        return self._beam_cache[key]

    def _draw_beam(self, trans, focus) -> None:
        """The section, laid out like the transducer schematic (bowl on the left, depth to the
        right): simulation medium, bowl, exit plane, scalp on the axis, target, focal depth
        setting, device steering range and acoustic FOV."""
        v, ax = self.v, self.beam_ax
        bowl = v._bowl()
        names, _ = self._media()
        dist = float(np.linalg.norm(focus - trans))
        exp = v._expected_bowl_distance()
        fov_diam, fov_len, _, _ = v._resolved_fov(trans, focus, bowl)
        a = bowl.aperture / 2
        depth_max = min(250.0, max(dist, exp or 0.0, fov_len or 0.0) + 15)
        half = float(np.ceil(max(a * 1.4, (fov_diam or 0) / 2 + 5, 30.0)))
        depth, offsets, t1, med, source = self._beam_section(trans, focus, depth_max, half)
        ext = (depth[0] - BEAM_STEP_MM / 2, depth[-1] + BEAM_STEP_MM / 2,
               offsets[0] - BEAM_STEP_MM / 2, offsets[-1] + BEAM_STEP_MM / 2)
        data = v.t1.canon
        lo, hi = np.percentile(data[::4, ::4, ::4], [1, 99.5])
        ax.imshow(t1, cmap="gray", origin="lower", vmin=lo, vmax=hi, extent=ext,
                  interpolation="nearest")
        water = names.index("water") + 1 if "water" in names else -1
        self._present = []
        if med is not None and names:
            colours = [MEDIUM_COLOURS.get(n, _EXTRA_COLOURS[i % len(_EXTRA_COLOURS)])
                       for i, n in enumerate(names)]
            cmap = ListedColormap(["#000000"] + colours).with_extremes(bad=(0, 0, 0, 0))
            shown = np.where((med > 0) & (med <= len(names)), med, np.nan).astype(float)
            ax.imshow(shown, cmap=cmap, vmin=0, vmax=len(names), origin="lower", extent=ext,
                      alpha=0.8, interpolation="nearest")
            self._present = [(names[i - 1], colours[i - 1]) for i in range(1, len(names) + 1)
                             if (med == i).any()]
        if fov_diam and fov_len:          # the FOV's true section in this plane
            ax.plot([0, fov_len, fov_len, 0, 0], [-fov_diam / 2, -fov_diam / 2, fov_diam / 2,
                                                 fov_diam / 2, -fov_diam / 2],
                    color="#43a047", lw=0.8, ls=(0, (3, 2)))
        dev = self.model.device(v.ti)
        if dev is not None and dev.min_foc is not None and dev.max_foc is not None:
            y = -a - 3                     # below the aperture, clear of the axis markers
            ax.plot([bowl.offset + dev.min_foc, bowl.offset + dev.max_foc], [y, y], color="#1e88e5",
                    lw=4, alpha=0.4, solid_capstyle="butt")
        ax.axhline(0, color="#e53935", lw=0.6, ls="--")
        R = bowl.curv_radius
        r = np.linspace(-a, a, 81)
        ax.plot(R - np.sqrt(np.maximum(R * R - r * r, 0)), r, color="#e53935", lw=2.5,
                solid_capstyle="butt")
        ax.plot([bowl.offset, bowl.offset], [-a * 1.1, a * 1.1], color="#e53935", lw=1, ls=":")
        # at the top edge, at the exit plane's depth: clear of the bowl's rim
        ax.text(bowl.offset, 0.99, " exit plane", fontsize=style.PLOT_PT["small"], color="#e53935", va="top",
                transform=ax.get_xaxis_transform())
        ax.plot([R], [0], "x", color="#757575", ms=6)
        ax.plot([dist], [0], "+", color="#1e88e5", ms=11, mew=2)
        if exp is not None:
            ax.plot([exp], [0], "x", color="#ff7043", ms=8, mew=2)

        scalp_depth = None
        if v.seg is not None:
            scalp = pl.scalp_exit(v.seg, focus, trans - focus)
            if scalp is not None:
                scalp_depth = dist - float(np.linalg.norm(scalp - focus))
                ax.plot([scalp_depth] * 2, [-a * 0.6, a * 0.6], color="#ff9800", lw=1.4)
                ax.annotate("", xy=(scalp_depth, -a * 0.75), xytext=(bowl.offset, -a * 0.75),
                            arrowprops=dict(arrowstyle="<->", color="#ff9800", lw=0.9))
                # beyond the scalp line, on a light box: clear of the bowl and the exit plane
                ax.text(max(scalp_depth, bowl.offset) + 1.5, -a * 0.75,
                        f"{scalp_depth - bowl.offset:+.1f} mm", fontsize=style.PLOT_PT["small"], color="#b26a00",
                        ha="left", va="center", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.2", fc=style.active()["surface"], ec="none", alpha=0.85))
                self._rows["contact"].append(("Exit plane → scalp", _mm_status(scalp_depth - bowl.offset)
                                              + muted(" on the axis")))
        mode, gap = self.model.get("placement.mode"), self.model.get("placement.skin_gap_mm")
        gap_ok = mode in GAP_MODES and isinstance(gap, (int, float)) and not isinstance(gap, bool)
        self._gap_line = gap_ok and scalp_depth is not None
        if self._gap_line:
            ax.plot([scalp_depth - gap] * 2, [-a, a], color="#ff9800", lw=0.9, ls="--")

        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_aspect("equal", adjustable="box")
        ax.tick_params(labelsize=style.PLOT_PT["small"], labelbottom=False)
        ax.set_ylabel("lateral [mm]", fontsize=style.PLOT_PT["normal"])

        if med is not None:
            on_axis, worst = pl.medium_clearance(med, depth, offsets, bowl, water)
            self._rows["contact"] += [
                ("Exit plane → medium", _mm_status(on_axis) + muted(" on the axis")),
                ("Closest in this plane", _mm_status(worst) + muted(" across the aperture"))]
        if gap_ok:
            self._rows["plane"].append(("Skin gap setting", f"{gap:g} mm" + muted(" (dashed line)")))
        if source == "run":
            when = datetime.fromtimestamp(self._medium_key[1]).strftime("%d %b %Y %H:%M")
            self._rows["plane"].append(("Medium", "last run" + muted(f" · {when}")))
        elif source == "labels":
            self._rows["plane"].append(("Medium", "SimNIBS labels" + muted(
                " via layers, unsmoothed " + ("(no run yet)" if self._outputs_read
                                              else "(looking for the last run's medium)"))))
        elif not self._layered():
            self._rows["missing"].append(("Medium", muted(f"simulation.medium is "
                                          f"{self.model.get('simulation.medium')}: no tissue layers")))
        if v.seg is None:
            self._rows["missing"].append(("Segmentation", muted("not found: no medium and no scalp "
                                                                "before a run")))
        self._rows["plane"].append(("Plane", f"{self.beam_angle:.0f}°"
                                    + muted(" · 0°: lateral axis towards superior")))

    # ── free-water profile ──────────────────────────────────────────────
    def _request_profile(self, start: bool = False) -> None:
        """Predict the free-water profile in MATLAB (nexus_calibration open + view). Without
        start only when MATLAB is already running, so a redraw never starts it."""
        key, bridge = self._profile_key(), self._bridge()
        if key is None or bridge is None or key in self._profiles or key == self._profile_pending:
            return
        if not start and (bridge.state not in ("ready", "busy") or key in self._profile_errors):
            return
        self._profile_errors.pop(key, None)
        self._profile_pending = key
        serial, ds, lib, ep = key
        isppa = self.model.get(f"{TX}.target_isppa_wcm2", self.v.ti)
        isppa = float(isppa) if isinstance(isppa, (int, float)) and not isinstance(isppa, bool) else 30.0

        def failed(err) -> None:
            self._profile_pending = None
            text = str(err).strip()
            self._profile_errors[key] = text.splitlines()[0] if text else "failed"
            self._draw()

        def viewed(out) -> None:
            self._profile_pending = None
            prof = profile_from_payload(out)
            if prof is None:
                self._profile_errors[key] = "MATLAB returned no profile"
            else:
                self._profiles[key] = prof
            self._draw()

        def opened(_out) -> None:
            bridge.call("nexus_calibration",
                        {"action": "view", "session": PROFILE_SESSION, "predicted": True,
                         "depths": [ep], "settings": {"isppa_wcm2": isppa}},
                        on_done=viewed, on_error=failed)

        bridge.call("nexus_calibration",
                    {"action": "open", "session": PROFILE_SESSION, "serial": serial,
                     "driving_system": ds, "library_path": lib},
                    on_done=opened, on_error=failed)
        if start:                                # a redraw calls this without start
            self._draw()

    def _curves(self) -> tuple[list[dict], str]:
        """The profiles to draw: dicts with z (exit-plane depth), I, label, style, absolute
        (W/cm²) and peak_kind; and their source ('run' | 'prediction' | '')."""
        if self._run_usable():
            run = self._run_profile
            z = np.asarray(run["depth_ep"], float)
            out = []
            if _finite(run.get("oneil")):
                out.append({"z": z, "I": np.asarray(run["oneil"], float), "label": "free water (O'Neil)",
                            "colour": FREE_WATER, "ls": "-", "lw": 1.6, "absolute": True})
            if _finite(run.get("kwave")):
                out.append({"z": z, "I": np.asarray(run["kwave"], float), "label": "free water (k-Wave)",
                            "colour": FREE_WATER, "ls": "--", "lw": 1.3, "absolute": True})
            front = z >= 0                       # behind the exit plane: the source
            if _finite(run.get("intensity")):
                out.append({"z": z[front], "I": np.asarray(run["intensity"], float)[front],
                            "label": "simulated, in the head", "colour": style.active()["ink"], "ls": "-",
                            "lw": 1.1, "absolute": True})
            return out, "run"
        key = self._profile_key()
        prof = self._profiles.get(key) if key is not None else None
        if prof is not None:
            return [{"z": np.asarray(prof["z_ep"], float), "I": np.asarray(prof["I"], float),
                     "label": f"free water ({prof['model']} model)", "colour": FREE_WATER,
                     "ls": "-", "lw": 1.6, "absolute": False}], "prediction"
        return [], ""

    def _draw_profile(self, trans, focus) -> None:
        """Profiles along the beam on the beam plane's depth axis, scaled to their own maxima or
        absolute [W/cm²], with the exit plane, target and focal depth setting."""
        ax, v = self.prof_ax, self.v
        bowl = v._bowl()
        ax.set_autoscalex_on(False)
        ax.tick_params(labelsize=style.PLOT_PT["small"])
        ax.set_xlabel("depth from bowl apex [mm]", fontsize=style.PLOT_PT["normal"])
        if not self._run_usable():
            self._request_profile()
        curves, source = self._curves()
        has_abs = any(c["absolute"] for c in curves)
        self.rb_absolute.setEnabled(has_abs)
        self.rb_absolute.setToolTip("Intensity in W/cm²: compare how much reaches the target"
                                    if has_abs else "No absolute profile yet: the prediction is "
                                    "relative. Run a simulation for profiles in W/cm².")
        key = self._profile_key()
        self.compute.setVisible(source != "run" and key is not None and key not in self._profiles)
        self.compute.setEnabled(key != self._profile_pending)
        absolute = self.scale == "absolute" and has_abs
        shown = [c for c in curves if c["absolute"] or not absolute]

        def note(text: str) -> None:
            # on a box of the plot's background: the depth lines behind it do not run through it;
            # broken into lines here (matplotlib's wrap follows the figure, not this narrow plot)
            import textwrap
            text = "\n".join(textwrap.fill(part, 55) for part in text.split("\n"))
            ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=style.PLOT_PT["normal"], color=style.active()["ink3"],
                    transform=ax.transAxes, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.5", fc=style.active()["surface"], ec="none"))

        ax.axvline(bowl.offset, color="#e53935", lw=1, ls=":")
        ax.axvline(float(np.linalg.norm(focus - trans)), color="#1e88e5", lw=1.3)
        exp = v._expected_bowl_distance()
        if exp is not None:
            ax.axvline(exp, color="#ff7043", lw=1.1, ls="--")
        if not shown:
            ax.set_yticks([])
            ax.set_ylabel("free water", fontsize=style.PLOT_PT["normal"])
            if not self._outputs_read:
                note("Looking for the last run's profile…")
            elif key is None and source == "":
                note("No free-water profile: it needs a run of this study, or a library device "
                     "with a focal depth setting.")
            elif key == self._profile_pending:
                note("Computing the free-water profile in MATLAB…")
            elif key in self._profile_errors:
                note(f"Free-water profile not available: {self._profile_errors[key]}\n"
                     "'Compute free-water profile' tries again.")
            elif self._bridge() is None:
                note("No free-water profile: run the study, or open the Nexus with MATLAB.")
            else:
                note("No free-water profile yet: run the study, or click "
                     "'Compute free-water profile' (starts MATLAB).")
            self._rows["water"].append(("Free water", muted("no profile yet")))
            return
        handles, top = [], 0.0
        for c in shown:
            i = c["I"]
            peak = np.nanmax(i) if np.isfinite(i).any() else np.nan
            y = i / peak if not absolute and peak > 0 else i
            handles.append(ax.plot(c["z"] + bowl.offset, y, color=c["colour"], ls=c["ls"],
                                   lw=c["lw"], label=c["label"])[0])
            top = max(top, float(np.nanmax(y)) if np.isfinite(y).any() else 0.0)
        if absolute:
            ax.set_ylabel("I [W/cm²]", fontsize=style.PLOT_PT["normal"])
            ax.set_ylim(0, (top or 1) * 1.15)
        else:
            ax.set_ylabel("I / own max", fontsize=style.PLOT_PT["normal"])
            ax.set_ylim(0, 1.15)
            ax.set_yticks([0, 0.5, 1])
            ax.axhline(0.5, color=style.active()["ink4"], lw=0.6, ls=":")
        ax.legend(handles=handles, loc="upper right", fontsize=style.PLOT_PT["small"], frameon=False)

        # free-water focus marker and readout (none when the run wrote no free-water column)
        fw = next((c for c in curves if c["label"].startswith("free water")), None)
        if fw is None:
            self._rows["water"].append(("Free water", muted("not in this run's profile")))
            return
        i = int(np.nanargmax(fw["I"]))
        if source == "run":
            depth, how = float(fw["z"][i]), "peak, " + fw["label"].split("(")[1].rstrip(")")
            when = datetime.fromtimestamp(self._run_profile["mtime"]).strftime("%d %b %Y %H:%M")
            caption = f"last run · {when}"
        else:
            prof = self._profiles[key]
            depth = prof["focus_ep"] if np.isfinite(prof["focus_ep"]) else float(fw["z"][i])
            how = f"FLHM centre, {prof['model']} model"
            caption = f"prediction at the focal depth setting ({key[3]:g} mm) · relative only"
        if any(c is fw for c in shown):
            y_top = ax.get_ylim()[1]
            ax.plot([depth + bowl.offset], [y_top * 0.93], "v", color=FREE_WATER, ms=6)
        differ = source == "run" and self._run_settings_differ()
        if differ:
            caption += " · ▲ settings differ from now"
        if self.scale == "absolute" and not has_abs:
            caption += " · Absolute needs a run"
        ax.text(0.01, 0.95, caption, ha="left", va="top", fontsize=style.PLOT_PT["small"],
                color=style.active()["amber" if differ else "ink3"], transform=ax.transAxes)
        value = f"{depth:.1f} mm" + muted(f" from the exit plane ({how})")
        if differ:
            value += "<br>" + status("the run used another placement or focal depth", "warn")
        self._rows["water"].append(("Free-water focus", value))
        if source == "run":                  # the run's own peak in brain, not the near field
            run = self._run_profile
            z, sim = np.asarray(run["depth_ep"], float), np.asarray(run["intensity"], float)
            brain = np.array([t == "brain" for t in run.get("tissue", [])], bool)
            if brain.size == sim.size and brain.any() and np.isfinite(sim[brain]).any():
                j = int(np.nanargmax(np.where(brain, sim, np.nan)))
                self._rows["water"].append(("Simulated peak in brain", f"{sim[j]:.1f} W/cm²"
                                            + muted(f" at {z[j]:.1f} mm from the exit plane, on the axis")))

    def _run_settings_differ(self) -> bool:
        """The last run used another focal depth setting or placement than the config now."""
        run, v = self._run_profile, self.v
        bowl, exp = v._bowl(), v._expected_bowl_distance()
        if bowl is not None and exp is not None and np.isfinite(run["planned_ep"]) \
                and abs(run["planned_ep"] - (exp - bowl.offset)) > 0.05:
            return True
        rec_trans, rec_focus = self._recorded_positions()
        trans, focus = v.positions()
        for a, b in ((rec_trans, trans), (rec_focus, focus)):
            if a is not None and b is not None and np.linalg.norm(a - b) > 1.0:
                return True
        return False

    def _recorded_positions(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        """The run's resolved positions of the shown transducer in RAS mm, from the record the
        background read found (as PlacementPage._recorded_positions_for, without the disk)."""
        rec, v = self._placement_record, self.v
        if rec is None or v.t1 is None or rec.get("mode") != self.model.get("placement.mode"):
            return None, None
        txs = rec["transducers"]
        if v.ti >= len(txs) or not isinstance(txs[v.ti], dict):
            return None, None
        out = []
        for key in ("trans_pos", "focus_pos"):
            p = txs[v.ti].get(key)
            out.append(pl.prestus_voxel_to_ras(v.t1, p) if isinstance(p, list) and len(p) == 3 else None)
        return out[0], out[1]

    def _draw_key(self) -> None:
        ax = self.key_ax
        marks = [
            Line2D([], [], ls="none", marker="+", color="#1e88e5", ms=10, mew=2, label="target (focus_pos)"),
            Line2D([], [], ls="none", marker="x", color="#ff7043", ms=8, mew=2, label="focal depth setting"),
            Line2D([], [], ls="none", marker="x", color="#757575", ms=6, label="geometric focus"),
            Line2D([], [], color="#e53935", lw=2.5, label="bowl"),
            Line2D([], [], color="#e53935", lw=1, ls=":", label="exit plane"),
            Line2D([], [], color="#e53935", lw=0.8, ls="--", label="beam axis"),
            Line2D([], [], color="#ff9800", lw=1.4, label="scalp on the axis"),
        ]
        if self._gap_line:
            marks.append(Line2D([], [], color="#ff9800", lw=0.9, ls="--", label="skin gap setting"))
        marks += [Line2D([], [], color="#1e88e5", lw=4, alpha=0.4, label="device steering range"),
                  Line2D([], [], color="#43a047", lw=0.8, ls=(0, (3, 2)), label="acoustic FOV"),
                  Line2D([], [], ls="none", marker="v", color=FREE_WATER, ms=6,
                         label="free-water focus")]
        first = ax.legend(handles=marks, loc="upper left", fontsize=style.PLOT_PT["normal"], frameon=False,
                          borderaxespad=0, handlelength=2)
        if self._present:
            ax.add_artist(first)
            ax.legend(handles=[Patch(color=c, label=n.replace("_", " ")) for n, c in self._present],
                      loc="lower left", fontsize=style.PLOT_PT["normal"], frameon=False, borderaxespad=0,
                      title="medium", title_fontsize=style.PLOT_PT["normal"], alignment="left")

    def _update_info(self) -> None:
        r = self._rows
        self.info.setText(readout([("Contact", r["contact"]), ("Section", r["plane"]),
                                   ("Free water", r["water"]), ("Not available", r["missing"])]))

    # ── interaction ─────────────────────────────────────────────────────
    def _scale_changed(self, _button=None) -> None:
        self.scale = "absolute" if self.rb_absolute.isChecked() else "scaled"
        self._draw()

    def _turn(self, deg: float) -> None:
        self.beam_angle = (self.beam_angle + deg) % 360
        self._draw()

    def _scrolled(self, event) -> None:
        if event.inaxes in (self.beam_ax, self.prof_ax):
            self._turn(BEAM_TURN_DEG if event.button == "up" else -BEAM_TURN_DEG)

    def _clicked(self, event) -> None:
        if event.inaxes is self.prof_ax:
            if not self._run_usable():
                self._request_profile(start=True)
            return
        if event.inaxes is not self.beam_ax or event.xdata is None:
            return
        trans, focus = self.v.positions()
        if trans is None or focus is None or np.linalg.norm(focus - trans) == 0:
            return
        axis, lat = pl.beam_plane_basis(trans, focus, self.beam_angle)
        self.v.cursor = trans + axis * event.xdata + lat * event.ydata
        self.v._cursor_set_by_user = True       # the Viewer's slices follow on its next draw
