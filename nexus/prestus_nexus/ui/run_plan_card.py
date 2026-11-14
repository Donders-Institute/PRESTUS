"""Run tab 'Plan' sub-tab: what a run of the current config will do — the pipeline mode, each
stage (runs / reuses / skipped / stops), the jobs, the outputs it always writes, the optional
outputs it saves or leaves out (and why), and the resources earlier runs used. Built from
run_plan.plan, predict_outputs.predict / save_items and estimate.estimate; cheap enough to
rebuild on every config edit."""

from __future__ import annotations

import html

from PyQt6.QtGui import QTextDocument
from PyQt6.QtWidgets import QFrame, QLabel, QTextBrowser, QVBoxLayout

from .. import run_plan as rp
from ..estimate import estimate
from ..model import ConfigModel, is_empty
from ..outputs import existing_outputs
from ..predict_outputs import predict, save_items
from . import icons, style

MODE_TEXT = {
    "single": "One pipeline run, start to finish.",
    "async": "Each transducer is simulated on its own; their intensities are combined for one "
             "joint thermal run.",
    "multi_isppa": "One acoustic simulation, then one thermal run per target ISPPA.",
    "uncertainty": "The default run plus four variants at the ITRUSST property bounds "
                   "(Temp-H/L for heating, Trans-H/L for transmission).",
}
NESTED_TEXT = {"async": " Each variant runs the async multi-transducer pipeline.",
               "multi_isppa": " Each variant runs the multi-ISPPA sweep."}
# stage state -> (icon, palette colour, word)
STATE = {rp.RUN: ("ok", "green", "runs"), rp.REUSE: ("reset", "info", "reuses existing"),
         rp.SKIP: ("minus", "ink4", "skipped"), rp.STOP: ("warning", "red", "stops the run")}
COMMON_DIRS = {"nii/": "NIfTI maps (T1w, and MNI when saved)", "img/": "figures",
               "log/": "pipeline log", "cache/": "intermediate .mat files", "debug/": "diagnostics"}


def _kind(name: str) -> str:
    return "folder" if name.endswith("/") else "file"


def plan_html(model: ConfigModel, doc: QTextDocument) -> str:
    """HTML of the plan; icons are registered as image resources of `doc`."""
    c = style.active()
    esc = html.escape
    img = lambda name, colour=None: icons.html_img(doc, name, colour)   # noqa: E731
    out = predict(model)
    stages = rp.plan(model)
    sim = model.get("path.sim")
    root = f"{sim.rstrip('/')}/{out.tag}/" if isinstance(sim, str) and not is_empty(sim) else f"<path.sim>/{out.tag}/"

    mode = out.mode.replace("_", "-")
    parts = [f"<h3>Mode: {esc(mode)}</h3>",
             f"<p>{esc(MODE_TEXT.get(out.mode, ''))}{esc(NESTED_TEXT.get(out.nested, ''))}</p>"]

    parts.append("<h3>Stages</h3><table cellspacing='0' cellpadding='3'>")
    for s in stages.stages:
        glyph, colour, word = STATE[s.state]
        ink = c["ink3"] if s.state == rp.SKIP else c["ink"]
        what = img(s.icon, c["ink3"] if s.state == rp.SKIP else c["ink2"]) if s.icon in icons.names() else ""
        parts.append(f"<tr><td>{img(glyph, c[colour])}</td><td>{what}</td>"
                     f"<td style='color:{ink}'><b>{esc(s.name)}</b></td>"
                     f"<td style='color:{c['ink3']}'>{esc(s.detail)}</td></tr>")
    parts.append("</table><p style='color:" + c["ink3"] + "'>"
                 + " &nbsp; ".join(f"{img(g, c[col])} {w}" for g, col, w in STATE.values()) + "</p>")

    parts.append("<h3>Jobs</h3><ul>" + "".join(f"<li>{esc(j)}</li>" for j in out.jobs) + "</ul>")
    for w in out.warnings:
        parts.append(f"<p style='color:{c['amber']}'>{img('warning', c['amber'])} {esc(w)}</p>")

    parts.append(f"<h3>Always written</h3><p>{img('folder')} <tt>{esc(root)}</tt></p>"
                 "<table cellspacing='0' cellpadding='2' style='margin-left:12px'>")
    entries = [(a.name, a.note) for a in out.artifacts] + [(d, COMMON_DIRS.get(d, "")) for d in out.common_dirs]
    for name, note in entries:
        parts.append(f"<tr><td>{img(_kind(name))}</td><td><tt>{esc(name)}</tt></td>"
                     f"<td style='color:{c['ink3']}'>{esc(note)}</td></tr>")
    parts.append("</table>")

    items = save_items(model)
    saved = [i for i in items if i.applies and i.saved]
    not_saved = [i for i in items if not (i.applies and i.saved)]
    sm = model.get("io.save_matrices")
    parts.append("<h3>Optional outputs</h3>"
                 f"<p style='color:{c['ink3']}'>Each is decided by its own <tt>io.save_*</tt> flag. "
                 "A matrix flag left empty follows <tt>io.save_matrices</tt> "
                 f"(now: {'empty' if is_empty(sm) else ('1, save' if bool(sm) else '0, do not save')}), "
                 "then the built-in default.</p>")

    def rows(title: str, group: list, colour: str) -> None:
        if not group:
            return
        parts.append(f"<p><b style='color:{c[colour]}'>{title}</b></p>"
                     "<table cellspacing='0' cellpadding='2' style='margin-left:12px'>")
        for i in group:
            ink = c["ink"] if colour == "green" else c["ink3"]
            what = f"{esc(i.label)}" + (f" — {esc(i.note)}" if i.note else "")
            parts.append(f"<tr><td valign='top'>{img(i.kind)}</td>"
                         f"<td style='color:{ink}'>{what}<br><tt>{esc(i.where)}</tt>"
                         f"<br><span style='color:{c['ink3']}'>{esc(i.reason)}</span></td></tr>")
        parts.append("</table>")

    on_disk = existing_outputs(model)
    if on_disk is not None:
        kinds = {i.label: i.kind for i in items}
        chip = lambda label: f"{img(kinds[label])}&nbsp;{esc(label)}"   # noqa: E731
        have = [k for k, v in on_disk.items() if v]
        missing = [k for k, v in on_disk.items() if not v]
        parts.append(f"<p><b>On disk now</b> <span style='color:{c['ink3']}'>(this subject, medium and "
                     "output affix, from earlier runs)</span></p>"
                     "<table cellspacing='0' cellpadding='2' style='margin-left:12px'>"
                     f"<tr><td valign='top' style='color:{c['green']}'>Exist</td><td>"
                     + (" &nbsp; ".join(chip(k) for k in have) or f"<span style='color:{c['ink3']}'>none</span>")
                     + f"</td></tr><tr><td valign='top' style='color:{c['ink3']}'>Not yet</td>"
                     f"<td style='color:{c['ink3']}'>"
                     + (" &nbsp; ".join(chip(k) for k in missing) or "none") + "</td></tr></table>")

    rows("Saved", saved, "green")
    rows("Not saved", not_saved, "ink3")

    res = [s.strip() for s in estimate(model).lines() if s.strip()]
    if res:
        parts.append("<h3>Resources (earlier runs of this project)</h3><ul>"
                     + "".join(f"<li>{esc(s)}</li>" for s in res) + "</ul>")
    return "".join(parts)


class RunPlanCard(QFrame):
    """'What this run does': the Run tab's Plan sub-tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        v = QVBoxLayout(self)
        title = QLabel("What this run does")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        sub = QLabel("From the current config, including unsaved edits. Runs use the saved file.")
        sub.setObjectName("CardSub")
        sub.setWordWrap(True)
        v.addWidget(sub)
        self.view = QTextBrowser()
        self.view.setOpenLinks(False)
        self.view.setFrameShape(QFrame.Shape.NoFrame)
        v.addWidget(self.view, 1)
        self._shown = ""

    def update_plan(self, model: ConfigModel) -> None:
        try:
            text = plan_html(model, self.view.document())
        except Exception as exc:          # a half-edited config must not break the Run tab
            text = f"<p>Could not summarise this config: {html.escape(str(exc))}</p>"
        if text != self._shown:
            self._shown = text
            bar = self.view.verticalScrollBar()
            pos = bar.value()
            self.view.setHtml(text)
            bar.setValue(pos)
