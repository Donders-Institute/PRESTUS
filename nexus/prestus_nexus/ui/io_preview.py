"""Setup tab side panel: the input files PRESTUS will look for and the output tree it will
write, resolved from path.*, subject_id and simulation.* exactly as get_output_dir.m /
resolve_io_dirs.m / predict_outputs.m do, so it stays an example rather than a promise.
Refreshed on every edit (see MainWindow._changed)."""

from __future__ import annotations

import html
import os

from PyQt6.QtWidgets import QLabel, QTextBrowser, QVBoxLayout, QWidget

from ..model import ConfigModel, is_empty
from ..paths import m2m_dir, matlab_sprintf
from ..predict_outputs import predict
from . import icons, style

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
        v.setContentsMargins(0, 8, 0, 0)
        title = QLabel("I/O preview")
        title.setObjectName("SubGroupTitle")
        v.addWidget(title)
        sub = QLabel("Example paths for the current subject ID, filled in as you edit the "
                     "settings. Detects async multi-transducer, multi-ISPPA and uncertainty "
                     "mode the same way prestus_pipeline_start.m does.")
        sub.setObjectName("CardSub")
        sub.setWordWrap(True)
        v.addWidget(sub)
        self.text = QTextBrowser()
        self.text.setOpenLinks(False)
        self.text.setLineWrapMode(QTextBrowser.LineWrapMode.NoWrap)
        self.text.setMinimumHeight(300)
        v.addWidget(self.text, 1)
        self._shown = ""

    def refresh(self) -> None:
        text = self._html()
        if text != self._shown:
            self._shown = text
            self.text.setHtml(text)

    def plain_text(self) -> str:
        """The preview without markup (tests, copy)."""
        return "\n".join(f"{label}  {path}  {note}".rstrip() for _k, label, path, note in self._rows()
                         if label or path)

    def _html(self) -> str:
        c = style.active()
        doc = self.text.document()
        mono = style.mono_font().family()
        out = ["<table cellspacing='0' cellpadding='2'>"]
        for kind, label, path, note in self._rows():
            if kind == "head":
                out.append(f"<tr><td colspan='3' style='padding-top:6px'><b>{html.escape(label)}</b></td></tr>")
                continue
            ic = icons.html_img(doc, kind, c["ink2"]) if kind in ("file", "folder") else \
                icons.html_img(doc, "warning", c["amber"]) if kind == "warning" else ""
            indent = "padding-left:18px;" if kind in ("file", "folder") and not label else ""
            out.append(f"<tr><td style='{indent}'>{ic}</td>"
                       f"<td style='color:{c['ink3']}'>{html.escape(label)}</td>"
                       f"<td><span style='font-family:\"{mono}\"'>{html.escape(path)}</span>"
                       + (f" <span style='color:{c['ink3']}'>{html.escape(note)}</span>" if note else "")
                       + "</td></tr>")
        out.append("</table>")
        return "".join(out)

    def _rows(self) -> list[tuple[str, str, str, str]]:
        """(kind, label, path, note); kind: head | file | folder | warning | text."""
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

        rows = [("head", f"Inputs ({tag})", "", "")]
        if layered:
            if not is_empty(anat) and sid is not None:
                rows.append(("file", "T1w", os.path.join(anat, matlab_sprintf(g("path.t1_pattern") or "", sid)), ""))
                t2_pattern = g("path.t2_pattern")
                if not is_empty(t2_pattern):
                    rows.append(("file", "T2w / UTE", os.path.join(anat, matlab_sprintf(t2_pattern, sid)), ""))
            else:
                rows.append(("text", "T1w", "(set path.anat)", ""))
        else:
            rows.append(("text", "", f"({medium} medium: no anatomical scan needed)", ""))
        if needs_seg:
            if not is_empty(seg_root) and sid is not None:
                rows.append(("folder", "segmentation", f"{m2m_dir(seg_root, sid)}/", ""))
            else:
                rows.append(("text", "segmentation", "(set path.anat or path.seg)", ""))
        if layered and bool(g("pct.enabled")):
            ct = g("path.ct_pattern")
            if not is_empty(ct):
                where = os.path.join(anat, matlab_sprintf(ct, sid)) if not is_empty(anat) and sid is not None \
                    else "(set path.anat)"
                rows.append(("file", "measured CT", where, "co-registered to the T1w"))
            elif sid is not None and not (is_empty(pct) and is_empty(seg_root)):
                base = pct if not is_empty(pct) else os.path.join(os.path.dirname(str(seg_root).rstrip("/\\")), "pct")
                rows.append(("folder", "pseudo-CT", f"{os.path.join(base, tag)}/", "generated from the UTE if missing"))

        plan = predict(self.model)
        if plan.mode == "single":
            mode_note = "single-run mode"
        else:
            mode_note = f"{plan.mode.replace('_', '-')} mode"
            if plan.nested:
                mode_note += f", {plan.nested.replace('_', '-')} per variant"
        rows.append(("head", f"Outputs ({mode_note})", "", ""))
        for warning in plan.warnings:
            rows.append(("warning", "", warning, ""))
        if is_empty(sim):
            rows.append(("text", "", "(set path.sim)", ""))
        else:
            rows.append(("folder", "output", f"{os.path.join(sim, tag)}/", ""))
            for a in plan.artifacts:
                rows.append(("file", "", a.name, a.note))
            for d in plan.common_dirs:
                rows.append(("folder", "", d, _COMMON_NOTES.get(d, "")))
            if not os.path.isabs(sim):
                rows.append(("warning", "", "path.sim is relative: PRESTUS resolves it against MATLAB's "
                             "working folder.", ""))
        return rows
