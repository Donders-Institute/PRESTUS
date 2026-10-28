"""Setup tab side panel: the input files PRESTUS will look for and the output tree it will
write, resolved from path.*, subject_id and simulation.* exactly as get_output_dir.m /
resolve_io_dirs.m / predict_outputs.m do, so it stays an example rather than a promise.
Refreshed on every edit (see MainWindow._changed)."""

from __future__ import annotations

import os

from PyQt6.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout, QWidget

from ..model import ConfigModel, is_empty
from ..paths import m2m_dir, matlab_sprintf
from ..predict_outputs import predict
from . import style

_COMMON_NOTES = {
    "nii/": "property/acoustic/thermal maps (T1w + MNI)",
    "img/": "intensity/temperature overlay PNGs",
    "log/": "pipeline diary log",
    "cache/": "regenerable intermediates (.mat)",
    "debug/": "preproc/ medium/ source/",
}


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


class IOPreviewPanel(QWidget):
    """A live example of the files this study will read and write."""

    def __init__(self, model: ConfigModel, parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model

        v = QVBoxLayout(self)
        title = QLabel("I/O setup")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        sub = QLabel("Example paths for the current subject ID, filled in as you edit the "
                     "Setup, Transducer and Simulation tabs. Detects async multi-transducer, "
                     "multi-ISPPA and uncertainty mode the same way prestus_pipeline_start.m does.")
        sub.setObjectName("CardSub")
        sub.setWordWrap(True)
        v.addWidget(sub)
        self.text = QPlainTextEdit(readOnly=True)
        self.text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        mono = style.mono_font()
        self.text.setFont(mono)
        v.addWidget(self.text, 1)

    def refresh(self) -> None:
        self.text.setPlainText("\n".join(self._lines()))

    def _lines(self) -> list[str]:
        g = self.model.get
        sid_raw = g("subject_id")
        sid = int(sid_raw) if _is_num(sid_raw) and sid_raw == int(sid_raw) and 0 <= sid_raw <= 999 else None
        tag = f"sub-{sid:03d}" if sid is not None else "sub-XXX"

        medium = g("simulation.medium") or "layered"
        layered = medium == "layered"
        needs_seg = layered or g("placement.mode") in ("mni", "heuristic") or \
            g("placement.manual.coordinate_space") == "mni"

        anat = g("path.anat")
        seg = g("path.seg")
        seg_root = anat if is_empty(seg) else seg
        pct = g("path.pct")
        sim = g("path.sim")

        lines = [f"INPUTS  ({tag})"]
        if layered:
            if not is_empty(anat) and sid is not None:
                t1 = os.path.join(anat, matlab_sprintf(g("path.t1_pattern") or "", sid))
                lines.append(f"  T1w           {t1}")
                t2_pattern = g("path.t2_pattern")
                if not is_empty(t2_pattern):
                    t2 = os.path.join(anat, matlab_sprintf(t2_pattern, sid))
                    lines.append(f"  T2w           {t2}")
            else:
                lines.append("  T1w           (set path.anat)")
        else:
            lines.append(f"  ({medium} medium: no anatomical scan needed)")
        if needs_seg:
            if not is_empty(seg_root) and sid is not None:
                lines.append(f"  segmentation  {m2m_dir(seg_root, sid)}/")
            else:
                lines.append("  segmentation  (set path.anat or path.seg)")
        if bool(g("pct.enabled")):
            if is_empty(pct):
                lines.append("  pseudo-CT     written into the segmentation m2m folder (path.pct unset)")
            elif sid is not None:
                lines.append(f"  pseudo-CT     {os.path.join(pct, tag)}/")

        plan = predict(self.model)
        if plan.mode == "single":
            mode_note = "single-run mode"
        elif plan.ambiguous:
            mode_note = f"{plan.mode.replace('_', '-')} mode assumed ({plan.platform_used} dispatch — see warning)"
        else:
            mode_note = f"{plan.mode.replace('_', '-')} mode"
        lines += ["", f"OUTPUTS  ({mode_note})"]
        for warning in plan.warnings:
            lines.append(f"  ! {warning}")
        if is_empty(sim):
            lines.append("  (set path.sim)")
        else:
            out = os.path.join(sim, tag)
            lines.append(f"  {out}/")
            entries = [(a.name, a.note) for a in plan.artifacts]
            entries += [(d, _COMMON_NOTES.get(d, "")) for d in plan.common_dirs]
            width = max(len(name) for name, _ in entries) + 2
            lines += [f"      {name:<{width}}{desc}" for name, desc in entries]
            if not os.path.isabs(sim):
                lines.append("")
                lines.append("  path.sim is relative: PRESTUS resolves it against MATLAB's working folder.")
        return lines
