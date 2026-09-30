"""Placement viewer: orthogonal T1 slices with segmentation overlay, click-to-place
transducer and focus, bowl outline and distance checks (manual placement)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QButtonGroup, QComboBox, QDoubleSpinBox, QFileDialog, QHBoxLayout,
                             QLabel, QMessageBox, QPushButton, QRadioButton, QSlider, QVBoxLayout,
                             QWidget)

from .. import paths, placement as pl
from ..model import TX, ConfigModel, is_empty

try:        # matplotlib is optional: the page explains what is missing
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False


# SimNIBS charm labels (final_tissues.nii.gz) and display colours
TISSUE = {1: ("white matter", "#f2f2f2"), 2: ("grey matter", "#9e9e9e"), 3: ("CSF", "#4fa3e0"),
          4: ("bone", "#f5deb3"), 5: ("scalp", "#e8a07a"), 6: ("eye", "#7ec8e3"),
          7: ("compact bone", "#fff3c4"), 8: ("spongy bone", "#d9b779"), 9: ("blood", "#c0392b"),
          10: ("muscle", "#b5651d")}
VIEWS = [("Axial", 2, (0, 1)), ("Coronal", 1, (0, 2)), ("Sagittal", 0, (1, 2))]
AXIS_LABELS = {0: ("L", "R"), 1: ("P", "A"), 2: ("I", "S")}
_OVERHEAD_FRAC = 0.12          # approx. fraction of row height spent on title/axis-label chrome


class PlacementPage(QWidget):
    changed = pyqtSignal()

    def __init__(self, model: ConfigModel, parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        self.ti = 0
        self.t1: pl.Volume | None = None
        self.seg: pl.Volume | None = None
        self._loaded_for: tuple | None = None
        self.cursor = np.zeros(3)
        self.direction: np.ndarray | None = None     # focus -> transducer, when no transducer yet
        self._slice_ratio: float | None = None       # ideal slice-row height per pixel of width

        v = QVBoxLayout(self)
        title = QLabel("Placement viewer")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        self.banner = QLabel()
        self.banner.setWordWrap(True)
        v.addWidget(self.banner)

        bar = QHBoxLayout()
        self.tx_combo = QComboBox()
        self.tx_combo.currentIndexChanged.connect(self._select_transducer)
        bar.addWidget(QLabel("Transducer"))
        bar.addWidget(self.tx_combo)
        self.tx_add = QPushButton("+")
        self.tx_add.setToolTip("Add a transducer")
        self.tx_add.setMaximumWidth(28)
        self.tx_add.clicked.connect(self._add_transducer)
        bar.addWidget(self.tx_add)
        self.tx_remove = QPushButton("−")
        self.tx_remove.setToolTip("Remove this transducer")
        self.tx_remove.setMaximumWidth(28)
        self.tx_remove.clicked.connect(self._remove_transducer)
        bar.addWidget(self.tx_remove)
        self.rb_focus = QRadioButton("Click sets focus")
        self.rb_trans = QRadioButton("Click sets transducer")
        self.rb_browse = QRadioButton("Click to browse")
        self.rb_browse.setToolTip("Left-click just moves the crosshair, like scrolling; "
                                  "nothing is written to the study config")
        self.rb_focus.setChecked(True)
        grp = QButtonGroup(self)
        for rb in (self.rb_focus, self.rb_trans, self.rb_browse):
            grp.addButton(rb)
            bar.addWidget(rb)
        bar.addStretch(1)
        self.use_manual = QPushButton("Switch to manual placement")
        self.use_manual.clicked.connect(self._use_manual)
        bar.addWidget(self.use_manual)
        v.addLayout(bar)

        tools = QHBoxLayout()
        b_axis = QPushButton("Transducer at focal distance")
        b_axis.setToolTip("Keep the beam direction; put the bowl apex focal_distance_bowl from the focus")
        b_axis.clicked.connect(self._place_on_axis)
        self.gap = QDoubleSpinBox(minimum=-20, maximum=50, value=5, singleStep=0.5, suffix=" mm")
        self.gap.setToolTip("Distance between scalp and exit plane")
        b_snap = QPushButton("Snap to scalp")
        b_snap.setToolTip("Keep the beam direction; put the exit plane this far outside the scalp")
        b_snap.clicked.connect(self._snap)
        b_match = QPushButton("Match focal distance")
        b_match.setToolTip("Set focal_distance_ep to the exit plane → focus distance of this placement")
        b_match.clicked.connect(self._match_focal_distance)
        b_traj = QPushButton("Import trajectory…")
        b_traj.setToolTip("Brainsight target export (NIfTI coordinates) or 3D Slicer transform")
        b_traj.clicked.connect(self._import_trajectory)
        b_centre = QPushButton("Centre on focus")
        b_centre.clicked.connect(self._centre)
        self.alpha = QSlider(Qt.Orientation.Horizontal, minimum=0, maximum=100, value=35)
        self.alpha.setMaximumWidth(120)
        self.alpha.valueChanged.connect(self._draw)
        for w in (b_axis, QLabel("gap"), self.gap, b_snap, b_match, b_traj, b_centre,
                  QLabel("overlay"), self.alpha):
            tools.addWidget(w)
        tools.addStretch(1)
        v.addLayout(tools)

        coords = QHBoxLayout()
        coords.addWidget(QLabel("Enter coordinates:"))
        self.coord_space = QComboBox()
        self.coord_space.addItem("RAS mm", "ras")
        self.coord_space.addItem("Voxel", "voxel")
        self.coord_space.addItem("MNI mm", "mni")
        self.coord_space.setCurrentIndex({"ras": 0, "voxel": 1, "mni": 2}.get(self._space(), 1))
        self.coord_space.setToolTip("MNI mm is converted to this subject's space with SimNIBS's "
                                    "mni2subject_coords (needs startup.simnibs_bin_path and the "
                                    "m2m folder); RAS and voxel need no conversion. For the focus "
                                    "this uses the nonlinear warp (most accurate in-brain); for "
                                    "the transducer it uses the affine, which also still works if "
                                    "the nonlinear transform is missing.")
        coords.addWidget(self.coord_space)
        self.cx = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        self.cy = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        self.cz = QDoubleSpinBox(minimum=-500, maximum=500, decimals=2)
        for label, spin in (("X", self.cx), ("Y", self.cy), ("Z", self.cz)):
            coords.addWidget(QLabel(label))
            coords.addWidget(spin)
        b_goto = QPushButton("Go to")
        b_goto.setToolTip("Move the crosshair here; writes nothing")
        b_goto.clicked.connect(self._goto_coords)
        b_set_focus = QPushButton("Set focus")
        b_set_focus.clicked.connect(lambda: self._set_coords("focus_pos"))
        b_set_trans = QPushButton("Set transducer")
        b_set_trans.clicked.connect(lambda: self._set_coords("trans_pos"))
        for w in (b_goto, b_set_focus, b_set_trans):
            coords.addWidget(w)
        coords.addStretch(1)
        v.addLayout(coords)

        body = QHBoxLayout()
        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.4), tight_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.axes = [self.fig.add_subplot(1, 3, i + 1) for i in range(3)]
            self.canvas.mpl_connect("button_press_event", self._clicked)
            self.canvas.mpl_connect("scroll_event", self._scrolled)
            body.addWidget(self.canvas, 3)
        else:
            body.addWidget(QLabel("Install matplotlib and nibabel to use the viewer."), 3)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.info.setMinimumWidth(280)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        body.addWidget(self.info, 1)
        v.addLayout(body, 1)
        self._body = body
        hint = QLabel("Click a slice to place the focus or transducer, or switch to 'Click to "
                      "browse' to just move around; scroll to move through slices either way. "
                      "Or type RAS, voxel or MNI coordinates above and use Go to / Set focus / "
                      "Set transducer. Positions are written to the study config immediately.")
        hint.setStyleSheet("color: gray;")
        v.addWidget(hint)

    def resizeEvent(self, event) -> None:      # noqa: N802
        super().resizeEvent(event)
        self._apply_row_cap()

    # ── state from the model ────────────────────────────────────────────
    def _space(self) -> str:
        return str(self.model.get("placement.manual.coordinate_space") or "voxel")

    def _editable(self) -> bool:
        return self.model.get("placement.mode") in ("manual", None, "") and self._space() in ("voxel", "ras")

    def _load_images(self) -> None:
        m = self.model
        sid = m.get("subject_id")
        anat, pattern = m.get("path.anat"), m.get("path.t1_pattern")
        seg_root = m.get("path.seg") or anat
        key = (sid, anat, pattern, seg_root)
        if key == self._loaded_for:
            return
        self._loaded_for = key
        self.t1 = self.seg = None
        if not isinstance(sid, (int, float)) or is_empty(anat):
            return
        t1 = pl.find_t1(anat, pattern or "", int(sid))
        if t1 is not None:
            try:
                self.t1 = pl.Volume(t1)
                self.cursor = self.t1.to_ras(np.array(self.t1.shape) / 2)[0]
            except Exception as exc:
                self.banner.setText(f"<span style='color:#c9302c'>Cannot read {t1}: {exc}</span>")
        seg = pl.find_segmentation(seg_root, int(sid)) if not is_empty(seg_root) else None
        if seg is not None:
            try:
                self.seg = pl.Volume(seg)
            except Exception:
                self.seg = None

    def positions(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        """(transducer, focus) in RAS mm for the current transducer, if known."""
        out = []
        for key in ("trans_pos", "focus_pos"):
            v = self.model.get(f"{TX}.{key}", self.ti)
            if not (isinstance(v, list) and len(v) == 3) or self.t1 is None:
                out.append(None)
            elif self._space() == "ras":
                out.append(np.asarray(v, float))
            elif self._space() == "voxel":
                out.append(pl.prestus_voxel_to_ras(self.t1, v))
            else:
                out.append(None)          # MNI: needs SimNIBS to map to subject space
        return out[0], out[1]

    def _write(self, key: str, ras: np.ndarray) -> None:
        if self._space() == "ras":
            value = [round(float(c), 2) for c in ras]
        else:
            value = pl.ras_to_prestus_voxel(self.t1, ras)
        self.model.set(f"{TX}.{key}", value, self.ti)
        self.changed.emit()

    def _bowl(self) -> pl.Bowl | None:
        txs = self.model.effective_transducers()
        return pl.Bowl.from_transducer(txs[self.ti]) if self.ti < len(txs) else None

    def _expected_bowl_distance(self) -> float | None:
        tx = self.model.effective_transducers()[self.ti]
        bowl = self._bowl()
        ep, fb = tx.get("focal_distance_ep"), tx.get("focal_distance_bowl")
        if isinstance(fb, (int, float)) and not isinstance(fb, bool):
            return float(fb)
        if isinstance(ep, (int, float)) and not isinstance(ep, bool) and bowl is not None:
            return float(ep) + bowl.offset
        return None

    # ── refresh ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        self.tx_combo.blockSignals(True)
        self.tx_combo.clear()
        self.tx_combo.addItems([f"Transducer {i + 1}" for i in range(self.model.n_transducers)])
        self.ti = min(self.ti, self.model.n_transducers - 1)
        self.tx_combo.setCurrentIndex(self.ti)
        self.tx_combo.blockSignals(False)
        self.tx_remove.setEnabled(self.model.n_transducers > 1)
        self._load_images()

        mode = self.model.get("placement.mode")
        msgs = []
        self.use_manual.setVisible(mode != "manual")
        if mode != "manual":
            msgs.append(f"placement.mode is <b>{mode}</b>: PRESTUS computes the positions. "
                        "The viewer shows configured positions only.")
        if self._space() == "mni":
            msgs.append("coordinate_space is <b>mni</b>: positions need SimNIBS to map to this "
                        "subject and are not shown. Use 'ras' or 'voxel' to place here.")
        if self.t1 is None:
            msgs.append("<span style='color:#c9302c'>No T1 image found. Set subject_id, "
                        "path.anat and path.t1_pattern.</span>")
        elif self.seg is None:
            msgs.append("No segmentation (m2m_sub-XXX/final_tissues.nii.gz): no overlay, no scalp snapping.")
        self.banner.setText("<br>".join(msgs))
        trans, focus = self.positions()
        if focus is not None and not self._cursor_set_by_user:
            self.cursor = focus.copy()
        self._draw()

    _cursor_set_by_user = False

    # ── drawing ─────────────────────────────────────────────────────────
    def _draw(self) -> None:
        if not HAVE_MPL:
            return
        for ax in self.axes:
            ax.clear()
            ax.set_xticks([])
            ax.set_yticks([])
        if self.t1 is None:
            self._slice_ratio = None
            self._apply_row_cap()
            self.canvas.draw_idle()
            self._update_info()
            return
        t1 = self.t1
        c = t1.ras_to_canon(self.cursor)[0]
        data = t1.canon
        vs = np.sqrt((t1.canon_affine[:3, :3] ** 2).sum(axis=0))
        lo, hi = np.percentile(data[::4, ::4, ::4], [1, 99.5])
        trans, focus = self.positions()
        bowl = self._bowl()
        shapes = []
        if trans is not None and focus is not None and bowl is not None and \
                np.linalg.norm(focus - trans) > 0:
            shapes = bowl.outline(trans, focus - trans)
        row_ratio = 0.0     # tallest per-column (height/width) need across the three views
        for ax, (title, plane, (h, w)) in zip(self.axes, VIEWS):
            k = int(np.clip(round(c[plane]), 0, data.shape[plane] - 1))
            sl = np.take(data, k, axis=plane)            # remaining axes in order (h, w)
            row_ratio = max(row_ratio, (vs[w] / vs[h]) * (sl.shape[1] / sl.shape[0]))
            ax.imshow(sl.T, cmap="gray", origin="lower", vmin=lo, vmax=hi,
                      aspect=vs[w] / vs[h], interpolation="nearest")
            if self.seg is not None and self.alpha.value() > 0:
                ax.imshow(self._overlay(plane, k, sl.shape).T, origin="lower",
                          cmap=_TISSUE_CMAP, vmin=0, vmax=10, alpha=self.alpha.value() / 100,
                          aspect=vs[w] / vs[h], interpolation="nearest")
            ax.set_title(f"{title}", fontsize=9)
            ax.set_xlabel(f"{AXIS_LABELS[h][0]}  ←→  {AXIS_LABELS[h][1]}", fontsize=8)
            ax.set_ylabel(f"{AXIS_LABELS[w][0]}  ←→  {AXIS_LABELS[w][1]}", fontsize=8)
            ax.axhline(c[w], color="#ffd54f", lw=0.5, alpha=0.6)
            ax.axvline(c[h], color="#ffd54f", lw=0.5, alpha=0.6)

            def proj(pts):
                p = t1.ras_to_canon(pts)
                return p[:, h], p[:, w]
            if shapes:
                x, y = proj(shapes["rim"])
                ax.plot(x, y, color="#e53935", lw=1)
                x, y = proj(np.vstack([trans, focus]))
                ax.plot(x, y, color="#e53935", lw=0.8, ls="--")
                x, y = proj(shapes["exit"][None])
                ax.plot(x, y, "s", color="#e53935", ms=3)
                x, y = proj(shapes["geom_focus"][None])
                ax.plot(x, y, "x", color="#bdbdbd", ms=5)
            for point, colour, marker in ((trans, "#e53935", "o"), (focus, "#1e88e5", "+")):
                if point is None:
                    continue
                x, y = proj(point[None])
                ax.plot(x, y, marker, color=colour, ms=9, mew=2, mfc="none" if marker == "o" else colour)
            ax.set_xlim(-0.5, sl.shape[0] - 0.5)
            ax.set_ylim(-0.5, sl.shape[1] - 0.5)
        # Row height needed per pixel of width, for _apply_row_cap: the slices keep a fixed
        # real-world (mm) aspect ratio, so with a plain Expanding size policy the canvas just
        # gets wider without the slices getting any bigger — the row's height is set once by
        # the surrounding layout and never revisited, and the aspect-locked axes just centre
        # themselves in the extra width as blank margin. _OVERHEAD_FRAC approximates the title
        # and axis-label chrome tight_layout reserves, as a fraction of the row height.
        self._slice_ratio = row_ratio / len(self.axes) / (1 - _OVERHEAD_FRAC)
        self.canvas.draw_idle()
        self._apply_row_cap()
        self._update_info()

    def _apply_row_cap(self) -> None:
        if not HAVE_MPL:
            return
        if not self._slice_ratio:
            self.canvas.setMaximumHeight(16_777_215)           # Qt's QWIDGETSIZE_MAX: no cap
            return
        self.layout().activate()                # make self._body.geometry() reflect this resize
        avail = self._body.geometry().height() or self.canvas.height()
        ideal = self.canvas.width() * self._slice_ratio
        self.canvas.setMaximumHeight(max(50, round(min(ideal, avail))))

    def _overlay(self, plane: int, k: int, shape2d: tuple) -> np.ndarray:
        ii, jj = np.meshgrid(np.arange(shape2d[0]), np.arange(shape2d[1]), indexing="ij")
        vox = np.zeros((ii.size, 3))
        others = [a for a in range(3) if a != plane]
        vox[:, others[0]] = ii.ravel()
        vox[:, others[1]] = jj.ravel()
        vox[:, plane] = k
        labels = self.seg.sample(self.t1.canon_to_ras(vox)).astype(float)
        labels[labels == 0] = np.nan
        return labels.reshape(shape2d)

    # ── info panel ──────────────────────────────────────────────────────
    def _update_info(self) -> None:
        trans, focus = self.positions()
        lines = [f"<b>Coordinate space:</b> {self._space()}"]
        if self.t1 is not None:
            lines.append(f"<b>T1:</b> {self.t1.path.name} "
                         f"({'×'.join(f'{s:.2g}' for s in self.t1.voxel_size)} mm)")
        for label, p in (("Focus", focus), ("Transducer (bowl apex)", trans)):
            if p is None:
                lines.append(f"<b>{label}:</b> not set")
                continue
            vox = pl.ras_to_prestus_voxel(self.t1, p) if self.t1 is not None else None
            tissue = self._tissue_at(p)
            lines.append(f"<b>{label}:</b> RAS [{p[0]:.2f}, {p[1]:.2f}, {p[2]:.2f}] mm · "
                         f"voxel {vox}{' · ' + tissue if tissue else ''}")
        bowl = self._bowl()
        if trans is not None and focus is not None and bowl is not None:
            axis = focus - trans
            dist = float(np.linalg.norm(axis))
            exp = self._expected_bowl_distance()
            ep = dist - bowl.offset
            line = f"<b>Bowl → focus:</b> {dist:.2f} mm (exit plane → focus {ep:.2f} mm)"
            if exp is not None:
                diff = dist - exp
                colour = "#3c763d" if abs(diff) <= 2 else "#c9302c"
                line += (f"<br>Configured focal distance (bowl): {exp:.2f} mm — "
                         f"<span style='color:{colour}'>Δ {diff:+.2f} mm</span>")
            dev = self.model.device(self.ti)
            if dev is not None and dev.min_foc is not None and dev.max_foc is not None:
                ok = dev.min_foc <= ep <= dev.max_foc
                line += (f"<br>Device focal range {dev.min_foc:g}–{dev.max_foc:g} mm: "
                         + ("<span style='color:#3c763d'>OK</span>" if ok
                            else "<span style='color:#c9302c'>outside</span>"))
            lines.append(line)
            if self.seg is not None and dist > 0:
                exit_pt = trans + axis / dist * bowl.offset
                scalp = pl.scalp_exit(self.seg, focus, trans - focus)
                if scalp is not None:
                    gap = float(np.dot(exit_pt - scalp, (trans - focus) / dist))
                    colour = "#3c763d" if gap >= 0 else "#c9302c"
                    lines.append(f"<b>Exit plane vs scalp:</b> <span style='color:{colour}'>"
                                 f"{gap:+.1f} mm</span> ({'outside' if gap >= 0 else 'inside the head'})")
        self.info.setText("<br><br>".join(lines))

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

    def _add_transducer(self) -> None:
        self.ti = self.model.add_transducer()
        self._cursor_set_by_user = False
        self.changed.emit()
        self.refresh()

    def _remove_transducer(self) -> None:
        self.model.remove_transducer(self.ti)
        self.ti = max(0, self.ti - 1)
        self._cursor_set_by_user = False
        self.changed.emit()
        self.refresh()

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
            self._draw()

    def _scrolled(self, event) -> None:
        if self.t1 is None or event.inaxes not in self.axes:
            return
        _, plane, _ = VIEWS[self.axes.index(event.inaxes)]
        c = self.t1.ras_to_canon(self.cursor)[0]
        c[plane] += 1 if event.button == "up" else -1
        self.cursor = self.t1.canon_to_ras(c[None])[0]
        self._cursor_set_by_user = True
        self._draw()

    def _centre(self) -> None:
        _, focus = self.positions()
        if focus is not None:
            self.cursor = focus.copy()
            self._draw()

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
        if not isinstance(sid, (int, float)) or is_empty(seg_root) or is_empty(simnibs):
            QMessageBox.information(self, "MNI coordinates", "Needs subject_id, path.anat (or "
                                    "path.seg) and startup.simnibs_bin_path, to run SimNIBS's "
                                    "mni2subject_coords.")
            return None
        m2m = paths.m2m_dir(seg_root, int(sid))
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
        self._draw()

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
            QMessageBox.information(self, "Placement", "Needs a beam direction (transducer, trajectory "
                                    "or segmentation) and focal_distance_ep or _bowl.")
            return
        self._write("trans_pos", pl.place_on_axis(self.positions()[1], d, dist))
        self.refresh()

    def _snap(self) -> None:
        if not self._need("focus"):
            return
        d, bowl = self._direction(), self._bowl()
        if self.seg is None or d is None or bowl is None:
            QMessageBox.information(self, "Placement", "Needs the segmentation, a beam direction "
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
            QMessageBox.information(self, "Placement", "Needs transducer, focus and transducer geometry.")
            return
        ep = round(float(np.linalg.norm(focus - trans)) - bowl.offset, 2)
        self.model.set(f"{TX}.focal_distance_ep", ep, self.ti)
        self.model.reset(f"{TX}.focal_distance_bowl", self.ti)      # derived from ep by PRESTUS
        self.changed.emit()
        self.refresh()

    def _import_trajectory(self) -> None:
        if not self._need("any"):
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import trajectory", "",
                                              "Trajectories (*.txt *.tfm);;All files (*)")
        if not path:
            return
        try:
            traj = pl.read_trajectory(path)
        except Exception as exc:
            QMessageBox.warning(self, "Import trajectory", str(exc))
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
                                   ". Set focal_distance_ep to place the transducer on its axis") + ".")

    def _use_manual(self) -> None:
        self.model.set("placement.mode", "manual")
        self.changed.emit()
        self.refresh()


if HAVE_MPL:
    _TISSUE_CMAP = ListedColormap(["#000000"] + [TISSUE[i][1] for i in range(1, 11)]
                                  ).with_extremes(bad=(0, 0, 0, 0))
