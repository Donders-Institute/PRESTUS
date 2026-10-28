"""Pseudo-CT tab: shows the generated pseudoCT.nii.gz (if any) over the subject's T1,
and, when the pipeline has been run with pct.enabled=1, the empirical in-skull
property distributions cached by medium_setup.m (sub-XXX_<medium>_medium_stats<affix>.json
in dir_cache — see functions/medium/medium_setup.m) next to the Default/Conservative/
Liberal bounds (simulation.uncertainty_<variant>, else the files shipped in config/uncertainty/)."""

from __future__ import annotations

import json
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QLabel, QVBoxLayout, QWidget

from .. import placement as pl
from ..model import is_empty
from ..outputs import dir_output, dir_pct, pseudoct_file
from ..uncertainty import variant_file
from ..yamlio import load_file
from . import style
from .empty_state import EmptyStack
from .ortho_view import OrthoViewer

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

PROPERTIES = [
    ("pseudoCT_HU", "pseudo-HU", None),
    ("density", "Density [kg/m3]", "density"),
    ("sound_speed", "Sound speed [m/s]", "sound_speed"),
    ("alpha_coeff", "Attenuation [dB/(cm.MHz)]", "alpha_coeff"),
]
BOUND_STYLE = {"run": ("This run", "#333333", "-"),
              "conservative": ("Conservative", "#1565c0", "--"),
              "liberal": ("Liberal", "#c62828", "--")}


class PctPage(QWidget):
    """T1 + pseudoCT overlay, plus the cached in-skull property histograms (if any)."""

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self._loaded_for: tuple | None = None

        outer = QVBoxLayout(self)
        title = QLabel("Pseudo-CT")
        title.setObjectName("CardTitle")
        outer.addWidget(title)
        content = QWidget()
        self.stack = EmptyStack(content)
        outer.addWidget(self.stack, 1)
        v = QVBoxLayout(content)
        v.setContentsMargins(0, 0, 0, 0)
        self.banner = QLabel()
        self.banner.setObjectName("CardSub")
        self.banner.setWordWrap(True)
        self.banner.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.banner)

        self.viewer = OrthoViewer()
        v.addWidget(self.viewer, 2)

        # The histograms appear only with a cached distribution; else one line says how
        # to get one (no four empty axes).
        self.stats_title = QLabel("In-skull property distributions (from the last run with "
                                  "pct.enabled=1)")
        self.stats_title.setObjectName("SubGroupTitle")
        v.addWidget(self.stats_title)
        self.stats_note = QLabel()
        self.stats_note.setObjectName("CardSub")
        self.stats_note.setWordWrap(True)
        v.addWidget(self.stats_note)
        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.2), tight_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.axes = [self.fig.add_subplot(1, 4, i + 1) for i in range(4)]
            v.addWidget(self.canvas, 1)
            self.canvas.hide()
        else:
            self.canvas = None
            self.stats_note.setText("Install matplotlib to see property distributions.")

    @property
    def model(self):
        return self.window.model

    # ── loading ─────────────────────────────────────────────────────────
    def refresh(self) -> None:
        m = self.model
        sid, anat = m.get("subject_id"), m.get("path.anat")
        pattern, pct = m.get("path.t1_pattern"), m.get("path.pct")
        seg = m.get("path.seg")
        key = (sid, anat, pattern, pct, seg, m.get("path.sim"), m.get("simulation.medium"))
        if key == self._loaded_for:
            return
        self._loaded_for = key
        self._load_images()
        self._load_stats()

    def _t1(self) -> Path | None:
        m = self.model
        sid, anat, pattern = m.get("subject_id"), m.get("path.anat"), m.get("path.t1_pattern")
        if is_empty(anat) or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None
        return pl.find_t1(anat, pattern or "", int(sid))

    def _goto_run(self) -> None:
        from .main_window import RUN
        self.window.show_special(RUN)

    def _load_images(self) -> None:
        msgs = []
        enabled = bool(self.model.get("pct.enabled"))
        t1 = self._t1()
        pct_dir = dir_pct(self.model)
        pct_file, legacy = pseudoct_file(self.model)
        sid = self.model.get("subject_id")
        sub = f"sub-{int(sid):03d}" if isinstance(sid, (int, float)) and not isinstance(sid, bool) else "this subject"
        # One panel for what is missing, most basic first; the page shows content otherwise.
        if not enabled and pct_file is None:
            self.stack.show_empty("Pseudo-CT is off for this study",
                                  "Skull properties come from the tissue segmentation. Turn on "
                                  "pct.enabled to derive them from a UTE/PETRA image instead.",
                                  "Turn on pseudo-CT…", lambda: self.window.goto("pct.enabled"))
        elif t1 is None:
            anat = self.model.get("path.anat")
            self.stack.show_empty("No T1 image yet" if is_empty(anat) else f"No T1 image found for {sub}",
                                  "The pseudo-CT is shown over the subject's T1 image. Set the "
                                  "anatomy folder (path.anat) and path.t1_pattern on Setup ▸ Paths.",
                                  "Set T1 image…", lambda: self.window.goto("path.anat"))
        elif pct_file is None:
            self.stack.show_empty(f"No pseudo-CT for {sub} yet",
                                  "Running the pipeline with pct.enabled on writes pseudoCT.nii.gz"
                                  + (f" to {pct_dir}" if pct_dir is not None else "") + ".",
                                  "Go to Run", self._goto_run)
        else:
            self.stack.show_content()
        self.viewer.set_anatomy(t1)
        if not enabled:
            msgs.append("pct.enabled is off: this study does not use the pseudo-CT shown here.")
        if pct_file is not None:
            self.viewer.set_overlay(pct_file)
            msgs.append(f"pseudo-CT: {pct_file}")
            if legacy:
                msgs.append(f"<span style='color:{style.active()['red']}'>This pseudo-CT is in the SimNIBS folder, "
                            "where PRESTUS no longer writes. Move it (with tissues_mask.nii.gz) to "
                            f"{pct_dir}.</span>")
        else:
            self.viewer.set_overlay(None)
        self.banner.setText("<br>".join(msgs))

    def _stats_file(self) -> Path | None:
        base = dir_output(self.model)
        if base is None or not (base / "cache").is_dir():
            return None
        sid = self.model.get("subject_id")
        if not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None
        matches = sorted((base / "cache").glob(f"sub-{int(sid):03d}_*_medium_stats*.json"),
                         key=lambda p: p.stat().st_mtime)
        return matches[-1] if matches else None

    def _bounds(self) -> dict[str, dict]:
        """{'conservative': medium_properties dict, 'liberal': medium_properties dict}."""
        root = self.window.prestus_root
        out = {}
        for variant in ("conservative", "liberal"):
            path = variant_file(self.window.model, root, variant)
            if path.is_file():
                try:
                    out[variant] = load_file(path).get("medium_properties", {})
                except Exception:
                    pass
        return out

    def _load_stats(self) -> None:
        if self.canvas is None:
            return
        for ax in self.axes:
            ax.clear()
        for legend in list(self.fig.legends):
            legend.remove()
        stats_file = self._stats_file()
        self.stats_title.setText("In-skull property distributions")
        if stats_file is None:
            self.canvas.hide()
            self.stats_note.setText("None cached for this subject yet: a run with pct.enabled on "
                                    "writes them (sub-XXX_&lt;medium&gt;_medium_stats.json in the "
                                    "output's cache folder).")
            self.stats_note.show()
            return
        try:
            stats = json.loads(stats_file.read_text())
        except Exception as exc:
            self.canvas.hide()
            self.stats_note.setText(f"Could not read {stats_file}: {exc}")
            self.stats_note.show()
            return
        self.stats_note.hide()
        self.canvas.show()
        bounds = self._bounds()
        reference = stats.get("reference", {}).get("skull", {})
        for ax, (key, label, bound_key) in zip(self.axes, PROPERTIES):
            prop = stats.get("properties", {}).get(key)
            ax.set_title(label, fontsize=9)
            if prop is None or not prop.get("histogram", {}).get("counts"):
                ax.set_xticks([]); ax.set_yticks([])
                ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes,
                        fontsize=8, color="gray")
                continue
            edges = prop["histogram"]["edges"]
            counts = prop["histogram"]["counts"]
            centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(counts))]
            widths = [edges[i + 1] - edges[i] for i in range(len(counts))]
            ax.bar(centers, counts, width=widths, color="#8d9db6", edgecolor="none")
            if bound_key is not None:
                lines = {"run": reference.get(bound_key)}
                lines.update({v: bounds[v]["skull"].get(bound_key) for v in ("conservative", "liberal")
                             if v in bounds and "skull" in bounds[v]})
                for kind, value in lines.items():
                    if value is None:
                        continue
                    label_, colour, ls = BOUND_STYLE[kind]
                    ax.axvline(value, color=colour, linestyle=ls, linewidth=1.2, label=label_)
            ax.tick_params(labelsize=7)
        handles, labels = self.axes[1].get_legend_handles_labels()
        if handles:
            self.fig.legend(handles, labels, loc="upper center", ncol=len(handles), fontsize=8,
                            bbox_to_anchor=(0.5, 1.08))
        self.canvas.draw_idle()
        self.stats_title.setText("In-skull property distributions · "
                                 f"{stats_file.name} (n={stats.get('n_voxels', '?')} skull voxels)")
