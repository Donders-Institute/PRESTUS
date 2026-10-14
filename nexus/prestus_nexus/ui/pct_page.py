"""Pseudo-CT tab: shows the generated pseudoCT.nii.gz (if any) over the subject's T1,
and, when the pipeline has been run with pct.enabled=1, the empirical in-skull
property distributions cached by medium_setup.m (sub-XXX_<medium>_medium_stats<affix>.json
in dir_cache — see functions/medium/medium_setup.m) next to the Default/Conservative/
Liberal bounds shipped in config/uncertainty/."""

from __future__ import annotations

import json
from pathlib import Path

from PyQt6.QtWidgets import QLabel, QVBoxLayout, QWidget

from .. import placement as pl
from ..model import is_empty
from ..outputs import dir_output, dir_pct, pseudoct_file
from ..yamlio import load_file
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

        v = QVBoxLayout(self)
        title = QLabel("Pseudo-CT")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        self.banner = QLabel()
        self.banner.setWordWrap(True)
        v.addWidget(self.banner)

        self.viewer = OrthoViewer()
        v.addWidget(self.viewer, 2)

        v.addWidget(QLabel("In-skull property distributions (from the last run with "
                           "pct.enabled=1, if cached)"))
        if HAVE_MPL:
            self.fig = Figure(figsize=(9, 3.2), tight_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.axes = [self.fig.add_subplot(1, 4, i + 1) for i in range(4)]
            v.addWidget(self.canvas, 1)
        else:
            self.canvas = None
            v.addWidget(QLabel("Install matplotlib to see property distributions."), 1)

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

    def _load_images(self) -> None:
        msgs = []
        if not bool(self.model.get("pct.enabled")):
            msgs.append("pct.enabled is off: this study does not use a pseudo-CT.")
        t1 = self._t1()
        self.viewer.set_anatomy(t1)
        if t1 is None:
            msgs.append("<span style='color:#c9302c'>No T1 image found. Set subject_id, "
                        "path.anat and path.t1_pattern.</span>")
        pct_dir = dir_pct(self.model)
        pct_file, legacy = pseudoct_file(self.model)
        if pct_file is not None:
            self.viewer.set_overlay(pct_file)
            msgs.append(f"pseudo-CT: {pct_file}")
            if legacy:
                msgs.append("<span style='color:#c9302c'>This pseudo-CT is in the SimNIBS folder, "
                            "where PRESTUS no longer writes. Move it (with tissues_mask.nii.gz) to "
                            f"{pct_dir}.</span>")
        else:
            self.viewer.set_overlay(None)
            if bool(self.model.get("pct.enabled")):
                msgs.append("No pseudoCT.nii.gz found yet at the expected location"
                            + (f" ({pct_dir})" if pct_dir is not None else "")
                            + " — generate it by running the pipeline, or set path.pct/path.seg.")
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
            path = Path(root) / "config" / "uncertainty" / f"config_medium_{variant}.yaml"
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
        if stats_file is None:
            for ax, (_, label, _) in zip(self.axes, PROPERTIES):
                ax.set_title(label, fontsize=9)
                ax.text(0.5, 0.5, "no cached data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=8, color="gray")
                ax.set_xticks([]); ax.set_yticks([])
            self.canvas.draw_idle()
            note = ("<br>No cached in-skull property distribution for this subject yet — "
                    "run the pipeline with pct.enabled=1 (adds a sub-XXX_&lt;medium&gt;_"
                    "medium_stats.json cache automatically).")
            self.banner.setText(self.banner.text() + note)
            return
        try:
            stats = json.loads(stats_file.read_text())
        except Exception as exc:
            self.banner.setText(self.banner.text() + f"<br>Could not read {stats_file}: {exc}")
            return
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
                    label_, colour, style = BOUND_STYLE[kind]
                    ax.axvline(value, color=colour, linestyle=style, linewidth=1.2, label=label_)
            ax.tick_params(labelsize=7)
        handles, labels = self.axes[1].get_legend_handles_labels()
        if handles:
            self.fig.legend(handles, labels, loc="upper center", ncol=len(handles), fontsize=8,
                            bbox_to_anchor=(0.5, 1.08))
        self.canvas.draw_idle()
        self.banner.setText(self.banner.text() + f"<br>Distribution cache: {stats_file.name} "
                            f"(n={stats.get('n_voxels', '?')} skull voxels)")
