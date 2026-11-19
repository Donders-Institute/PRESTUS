"""Skull CT tab: shows the CT the skull mapping reads over the subject's T1 (the measured CT in
path.ct_pattern, else the generated pseudoCT.nii.gz, if any), and the in-skull property
distributions the default (CT-mapped) run models against the range the uncertainty bound
variants model.

With pct.enabled, only the default variant maps the skull from the CT; the four bound
variants keep its head geometry but give the skull layer the homogeneous values of their medium
file (simulation.uncertainty_<variant>, else the files shipped in config/uncertainty/). The
histograms come from, in order of preference:
  - the default run's own cache, sub-XXX_<medium>_medium_stats<affix>.json in dir_cache
    (functions/medium/medium_setup.m, written for every model with a skull layer: skull_model
    'pct', 'ct', or 'layered'; beam_path holds the same distributions on the beam path), else
    the newest sub-XXX_*_medium_stats*.json there;
  - the prediction sub-XXX_<medium>_pct_prediction<affix>.json that
    pct_predict_skull_properties.m writes without simulating (the "Predict skull properties"
    button, nexus_pct through the bridge).
The band spans the minimum to maximum skull value over the four variant files, and each axis
states the fraction of skull voxels outside it."""

from __future__ import annotations

import json
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

from .. import placement as pl
from ..bridge import BridgeError, as_list
from ..model import is_empty
from ..outputs import dir_output, dir_pct, measured_ct_file, pseudoct_file
from ..study import sanitize_affix
from ..uncertainty import property_sets
from . import style
from .empty_state import EmptyStack
from .export import subject_stem
from .matlab_status import MatlabStatus
from .ortho_view import OrthoViewer

try:
    from .mpl_theme import ThemedCanvas as FigureCanvasQTAgg   # drawn in the Nexus theme
    from matplotlib.figure import Figure
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

PROPERTIES = [
    ("pseudoCT_HU", "(pseudo-)HU", None),
    ("density", "Density [kg/m³]", "density"),
    ("sound_speed", "Sound speed [m/s]", "sound_speed"),
    ("alpha_coeff", "Attenuation [dB/(cm.MHz)]", "alpha_coeff"),
]
# reference.skull is the config's homogeneous skull value, which a mapped pCT run does not use
RUN_STYLE = ("Skull layer value (unused when mapped)", "ink", "--")      # colour: palette role
BAND_STYLE = ("Bound variants (homogeneous skull)", "#d98a1f", 0.22)
SOURCE_TEXT = {"grid": "grid voxels: the exact preprocessed grid",
               "native": "native voxels: the CT image before cropping/resampling/smoothing"}


def skull_band(bounds: dict[str, dict], key: str) -> tuple[float, float] | None:
    """(min, max) of the skull value over the variant files that set it."""
    values = []
    for props in bounds.values():
        v = (props.get("skull") or {}).get(key)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            values.append(float(v))
    return (min(values), max(values)) if values else None


def fraction_outside(edges: list, counts: list, lo: float, hi: float) -> float | None:
    """Share of voxels in bins whose centre lies outside [lo, hi]."""
    total = sum(counts)
    if not total:
        return None
    out = sum(n for i, n in enumerate(counts) if not lo <= (edges[i] + edges[i + 1]) / 2 <= hi)
    return out / total


def constant_value(prop: dict) -> float | None:
    """The single value of a distribution with min == max (e.g. attenuation 'k-plan'), else None."""
    lo, hi = prop.get("min"), prop.get("max")
    if isinstance(lo, (int, float)) and isinstance(hi, (int, float)) and lo == hi:
        return float(lo)
    return None


def histogram_bars(prop: dict) -> tuple[list, list, list] | None:
    """(centres, counts, widths) of a distribution_stats histogram; None without voxels.

    jsonencode writes a one-bin histogram as scalars, so both arrays are read as lists. A
    constant distribution becomes one bar at its value; its width is set when drawn."""
    hist = prop.get("histogram") or {}
    counts, edges = as_list(hist.get("counts")), as_list(hist.get("edges"))
    total = sum(n for n in counts if isinstance(n, (int, float)))
    if not total:
        return None
    value = constant_value(prop)
    if value is not None or len(edges) != len(counts) + 1:
        if value is None:
            return None
        return [value], [total], [0.0]
    return ([(edges[i] + edges[i + 1]) / 2 for i in range(len(counts))], counts,
            [edges[i + 1] - edges[i] for i in range(len(counts))])


def outside_share(prop: dict, lo: float, hi: float) -> float | None:
    """Share of skull voxels outside [lo, hi]: exact for a constant distribution, else by bin centre."""
    hist = prop.get("histogram") or {}
    counts, edges = as_list(hist.get("counts")), as_list(hist.get("edges"))
    if not sum(counts):
        return None
    value = constant_value(prop)
    if value is not None:
        return 0.0 if lo <= value <= hi else 1.0
    if len(edges) != len(counts) + 1:
        return None
    return fraction_outside(edges, counts, lo, hi)


def stats_candidates(cache: Path, subject_id: int, medium: str, affix: str) -> tuple[Path | None, Path | None]:
    """(default run's medium_stats, pct prediction) in a cache folder, either None if absent.

    affix: io.output_affix as written. The default variant's own file
    (sub-XXX_<medium>_medium_stats<affix>.json, medium_setup.m; affix as load_parameters writes
    it, else as written: study.affix_path) is preferred; else the newest
    sub-XXX_*_medium_stats*.json, another run's (stats_title says so). The prediction is matched
    exactly (pct_predict_skull_properties.m)."""
    from ..study import affix_path
    if not cache.is_dir():
        return None, None
    sub = f"sub-{subject_id:03d}"
    real = affix_path(lambda a: cache / f"{sub}_{medium}_medium_stats{a}.json", affix)
    affix = sanitize_affix(affix)
    if not real.is_file():
        matches = sorted(cache.glob(f"{sub}_*_medium_stats*.json"), key=lambda p: p.stat().st_mtime)
        real = matches[-1] if matches else None
    pred = cache / f"{sub}_{medium}_pct_prediction{affix}.json"
    return real, pred if pred.is_file() else None


def stats_title(stats: dict, path: Path, own: bool = True) -> str:
    """Title line: whether the histograms are the default run's mapped skull or a prediction;
    own=False: the file is another run's (no stats file for this config's affix)."""
    n = stats.get("n_voxels", "?")
    if isinstance(n, float) and n.is_integer():
        n = int(n)                                     # JSON from MATLAB: 1606731.0
    if stats.get("predicted"):
        source = stats.get("source", "")
        what = f"predicted from the CT ({SOURCE_TEXT.get(source, f'{source} voxels')})"
    elif stats.get("skull_model") == "layered":
        what = "layered skull of the default run"
    else:
        what = "mapped skull of the default run"
    if not own and not stats.get("predicted"):
        what += ", from another run (none for this config's output affix)"
    return f"In-skull properties, {what} (n={n} skull voxels)"      # the file: the tooltip


class PctPage(QWidget):
    """T1 + CT overlay, plus the in-skull property histograms of the default run or a
    prediction (if any) against the bound variants' homogeneous skull values."""

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self._loaded_for: tuple | None = None

        outer = QVBoxLayout(self)                    # the sub-tab label names the page
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
        if getattr(getattr(self.viewer, "canvas", None), "export_button", None) is not None:
            self.viewer.canvas.export_button.view = "pct_slices"
            self.viewer.canvas.export_button.filename_hint = subject_stem(lambda: self.model, "pct_slices")
        v.addWidget(self.viewer, 2)

        # The histograms appear only with a cached distribution; else one line says how
        # to get one (no four empty axes).
        head = QHBoxLayout()
        self.stats_title = QLabel("In-skull property distributions")
        self.stats_title.setObjectName("SubGroupTitle")
        self.stats_title.setWordWrap(True)
        head.addWidget(self.stats_title, 1)
        self.predict_btn = QPushButton("Predict skull properties")
        self.predict_btn.setToolTip("Apply the pct.mapping_* settings to the CT's skull voxels, "
                                    "without simulating (pct_predict_skull_properties in MATLAB)")
        self.predict_btn.clicked.connect(self.predict)
        head.addWidget(self.predict_btn)
        # MATLAB starts only when a prediction is asked for: showing this page reads files only.
        bridge = window.bridge()
        self.matlab_status = MatlabStatus(bridge)
        head.addWidget(self.matlab_status)
        bridge.state_changed.connect(self._idle_status)     # after MatlabStatus._update
        self._idle_status(getattr(bridge, "state", ""))
        self._predicting = False
        self._has_pct = False
        self.predict_btn.setEnabled(False)
        self.stats_file: Path | None = None
        self.stats_kind = ""
        v.addLayout(head)
        self.stats_note = QLabel()
        self.stats_note.setObjectName("CardSub")
        self.stats_note.setWordWrap(True)
        v.addWidget(self.stats_note)
        if HAVE_MPL:
            # constrained: the legend placed 'outside' above the panels gets its own room
            self.fig = Figure(figsize=(9, 3.2), layout="constrained")
            self.canvas = FigureCanvasQTAgg(self.fig, export_name="pct_properties",
                                            filename_hint=subject_stem(lambda: self.model, "pct_properties"))
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
        key = (sid, anat, pattern, pct, seg, m.get("path.ct_pattern"), m.get("path.sim"), m.get("simulation.medium"),
               m.get("io.output_affix"))
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
        measured = not is_empty(self.model.get("path.ct_pattern"))
        if measured:
            pct_file, legacy = measured_ct_file(self.model), False
        else:
            pct_file, legacy = pseudoct_file(self.model)
        sid = self.model.get("subject_id")
        sub = f"sub-{int(sid):03d}" if isinstance(sid, (int, float)) and not isinstance(sid, bool) else "this subject"
        # One panel for what is missing, most basic first; the page shows content otherwise.
        if not enabled and pct_file is None:
            self.stack.show_empty("CT-based skull mapping is off for this study",
                                  "Skull properties come from the tissue segmentation. Turn on "
                                  "pct.enabled to map them from a measured CT, or from a pseudo-CT "
                                  "generated from a UTE/PETRA image.",
                                  "Turn on CT mapping…", lambda: self.window.goto("pct.enabled"))
        elif t1 is None:
            anat = self.model.get("path.anat")
            self.stack.show_empty("No T1 image yet" if is_empty(anat) else f"No T1 image found for {sub}",
                                  "The CT is shown over the subject's T1 image. Set the "
                                  "anatomy folder (path.anat) and path.t1_pattern on Setup ▸ Paths.",
                                  "Set T1 image…", lambda: self.window.goto("path.anat"))
        elif pct_file is None and measured:
            self.stack.show_empty(f"No measured CT found for {sub}",
                                  "path.ct_pattern matches no file in the anatomy folder (path.anat).",
                                  "Set CT file…", lambda: self.window.goto("path.ct_pattern"))
        elif pct_file is None:
            self.stack.show_empty(f"No pseudo-CT for {sub} yet",
                                  "Running the pipeline with pct.enabled on writes pseudoCT.nii.gz"
                                  + (f" to {pct_dir}" if pct_dir is not None else "") + ".",
                                  "Go to Run", self._goto_run)
        else:
            self.stack.show_content()
        self.viewer.set_anatomy(t1)
        if not enabled:
            msgs.append("pct.enabled is off: this study does not use the CT shown here.")
        if pct_file is not None:
            self.viewer.set_overlay(pct_file)
            msgs.append(f"{'measured CT (shown at its own resolution; the run resamples it to the segmentation grid)' if measured else 'pseudo-CT'}: {pct_file}")
            if legacy:
                msgs.append(f"<span style='color:{style.active()['red']}'>This pseudo-CT is in the SimNIBS folder, "
                            "where PRESTUS no longer writes. Move it (with tissues_mask.nii.gz) to "
                            f"{pct_dir}.</span>")
        else:
            self.viewer.set_overlay(None)
        self.banner.setText("<br>".join(msgs))
        self._has_pct = pct_file is not None
        self._update_predict_btn()

    # ── prediction (MATLAB) ─────────────────────────────────────────────
    def _idle_status(self, state: str) -> None:
        """While MATLAB is off, say which action starts it (nothing starts it on show)."""
        if state == "off":
            c = style.active()
            self.matlab_status.label.setText(f"<span style='color:{c['ink3']}'>●</span> MATLAB starts "
                                             "when you predict")

    def _update_predict_btn(self) -> None:
        has_matlab = self.window.matlab_path() is not None
        self.predict_btn.setEnabled(self._has_pct and has_matlab and not self._predicting)
        if not self._has_pct:
            self.predict_btn.setToolTip("Needs a CT for this subject: a measured CT (path.ct_pattern) "
                                        "or a pseudoCT.nii.gz.")
        elif not has_matlab:
            self.predict_btn.setToolTip("MATLAB not found; set its path on the Run tab.")
        else:
            self.predict_btn.setToolTip("Apply the pct.mapping_* settings to the CT's skull voxels, "
                                        "without simulating (pct_predict_skull_properties in MATLAB)")

    def _saved_config(self) -> str | None:
        """The prediction reads the saved study config, as a run would."""
        if (self.window.dirty or self.model.path is None) and not self.window.save():
            return None
        return str(self.model.path)

    def predict(self) -> None:
        cfg = self._saved_config()
        if cfg is None:
            return
        self._predicting = True
        self._update_predict_btn()
        self._note("⟳ Predicting skull properties from the CT…", "accent")
        self.window.bridge().call("nexus_pct", {"action": "predict", "config": cfg},
                                  on_done=self._predicted, on_error=self._predict_failed)

    def _predicted(self, out: dict) -> None:
        self._predicting = False
        self._update_predict_btn()
        self._load_stats()
        if self.stats_file is not None and self.stats_kind == "real":
            self._note(f"Prediction written to {Path(str(out.get('file', ''))).name}; the default run's "
                       f"own stats are shown ({self.stats_file.name}).", "ink2")

    def _predict_failed(self, err: BridgeError) -> None:
        self._predicting = False
        self._update_predict_btn()
        self._note(f"✗ Prediction failed: {err.message}", "red")

    def _note(self, text: str, colour: str = "ink2") -> None:
        self.stats_note.setText(f"<span style='color:{style.active()[colour]}'>{text}</span>")
        self.stats_note.show()

    def _stats_files(self) -> tuple[Path | None, Path | None]:
        """(default run's medium_stats, prediction) for this subject, see stats_candidates."""
        base = dir_output(self.model)
        sid = self.model.get("subject_id")
        if base is None or not isinstance(sid, (int, float)) or isinstance(sid, bool):
            return None, None
        medium = str(self.model.get("simulation.medium") or "")
        return stats_candidates(base / "cache", int(sid), medium, self.model.get("io.output_affix") or "")

    def _bounds(self) -> dict[str, dict]:
        """{variant: medium_properties of its bounds, edits included} for each readable variant."""
        return dict(property_sets(self.window.model, self.window.prestus_root).listed)

    def _load_stats(self) -> None:
        if self.canvas is None:
            return
        for ax in self.axes:
            ax.clear()
        for legend in list(self.fig.legends):
            legend.remove()
        real, pred = self._stats_files()
        stats_file = real or pred
        self.stats_file, self.stats_kind = stats_file, "real" if real is not None else "predicted"
        self.stats_title.setText("In-skull property distributions")
        if stats_file is None:
            self.canvas.hide()
            self.stats_note.setText("None for this subject yet: predict them from the CT, or "
                                    "run with pct.enabled on (the default run writes "
                                    "sub-XXX_&lt;medium&gt;_medium_stats.json to the output's cache folder).")
            self.stats_note.show()
            return
        try:
            stats = json.loads(stats_file.read_text())
        except Exception as exc:
            self.canvas.hide()
            self.stats_note.setText(f"Could not read {stats_file}: {exc}")
            self.stats_note.show()
            return
        if real is not None and pred is not None:
            self._note(f"Showing the default run's mapped skull; a prediction also exists ({pred.name}).")
        else:
            self.stats_note.hide()
        self.canvas.show()
        bounds = self._bounds()
        reference = stats.get("reference", {}).get("skull", {})
        for ax, (key, label, bound_key) in zip(self.axes, PROPERTIES):
            prop = stats.get("properties", {}).get(key)
            ax.set_title(label.replace(" [", "\n[", 1), fontsize=style.PLOT_PT["title"])   # unit below: narrow panels
            bars = histogram_bars(prop) if isinstance(prop, dict) else None
            if bars is None:
                ax.set_xticks([]); ax.set_yticks([])
                ax.text(0.5, 0.5, "n/a", ha="center", va="center", transform=ax.transAxes,
                        fontsize=style.PLOT_PT["normal"], color=style.active()["ink3"])
                continue
            centers, counts, widths = bars
            band = skull_band(bounds, bound_key) if bound_key is not None else None
            if len(centers) == 1 and widths[0] == 0:
                # a constant mapping: one narrow bar, wide enough to see next to the band
                value = centers[0]
                span = max((band[1] - band[0]) if band else 0.0, abs(value) * 0.05, 1e-6)
                widths = [span * 0.03]
                lo_x = min(value, band[0]) if band else value
                hi_x = max(value, band[1]) if band else value
                ax.set_xlim(lo_x - span * 0.25, hi_x + span * 0.25)
                ax.text(0.02, 0.84, f"constant {value:g}", ha="left", va="top",   # below 'outside' 
                        transform=ax.transAxes, fontsize=style.PLOT_PT["normal"], color=style.active()["ink2"])
            ax.bar(centers, counts, width=widths, color="#8d9db6", edgecolor="none")
            if bound_key is not None:
                # each ITRUSST property has a lower and an upper bound, so per-variant lines
                # would sit on the band edges: the band alone shows them
                if band is not None:
                    lo, hi = band
                    label_, colour, alpha = BAND_STYLE
                    ax.axvspan(lo, hi, color=colour, alpha=alpha, linewidth=0, label=label_, zorder=0)
                    frac = outside_share(prop, lo, hi)
                    if frac is not None:
                        ax.text(0.98, 0.96, f"{frac:.0%} outside", ha="right", va="top",
                                transform=ax.transAxes, fontsize=style.PLOT_PT["normal"], color=style.active()["ink2"])
                run = reference.get(bound_key)
                if run is not None:
                    label_, role, ls = RUN_STYLE
                    ax.axvline(run, color=style.active()[role], linestyle=ls, linewidth=1.2, label=label_)
            ax.tick_params(labelsize=style.PLOT_PT["small"])
            # four panels share a narrow pane: the shape counts, not the voxel counts, and at
            # most three values per axis (no run-together tick labels or 1e6 offsets)
            from matplotlib.ticker import MaxNLocator
            ax.set_yticks([])
            ax.xaxis.set_major_locator(MaxNLocator(3))
        handles, labels = self.axes[1].get_legend_handles_labels()
        if handles:
            self.fig.legend(handles, labels, loc="outside upper center", ncol=len(handles),
                            fontsize=style.PLOT_PT["small"])
        self.canvas.draw_idle()
        raw = self.model.get("io.output_affix") or ""
        own_names = {f"_medium_stats{a}.json" for a in (sanitize_affix(raw), raw) if isinstance(a, str)}
        own = stats_file is not real or any(stats_file.name.endswith(n) for n in own_names)
        self.stats_title.setText(stats_title(stats, stats_file, own))
        self.stats_title.setToolTip(str(stats_file))
