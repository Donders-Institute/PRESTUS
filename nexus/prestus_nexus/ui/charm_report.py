"""Medium ▸ Head model: SimNIBS charm's report (charm_report.html in the m2m folder) as a card
of its own, below the SimNIBS segmentation card, shown while a layered medium's segmentation
has one.

The report is charm's quality check: the final segmentation over the T1, the initial atlas
registration and the T1–T2 registration, as slice images in a small viewer. The card's fold
starts closed; the web view is built and the page loaded only when it is opened. The page is
loaded from its file: the report is ~10 MB+, above QWebEngineView.setHtml's 2 MB limit.
Without PyQt6-WebEngine the card says so; 'Open charm report' on the segmentation card opens
it in the browser.

A stylesheet injected when the document is ready (report_css) gives the report the Nexus's
look: system UI font at the card text size, headings as card titles, the active palette for
text, links and controls. The slice viewer itself (brainsprite canvas) keeps its layout.
"""

from __future__ import annotations

import json
from pathlib import Path

from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import QFrame, QHBoxLayout, QLabel, QToolButton, QVBoxLayout, QWidget

from .. import segmentation as sg
from . import icons, style

VIEW_HEIGHT = 720                       # px: the report's viewer and slice panels without scrolling
STYLE_ID = "prestus-nexus-style"
# the Nexus's sizes (style.py: QWidget 13px, #CardTitle 14px/600, #CardSub ink3)
FONT = 'system-ui, -apple-system, "Segoe UI", "Helvetica Neue", Arial, sans-serif'


def report_css(c: dict) -> str:
    """CSS that restyles charm_report.html in palette `c` (style.active())."""
    return f"""
html, body {{ background: {c['surface']}; color: {c['ink']}; margin: 0; }}
body, div, form, p {{ font-family: {FONT} !important; font-size: 13px; line-height: 1.45; }}
#content {{ padding: 4px 8px 12px 0; }}
h1, h2, h3 {{ font-family: {FONT}; font-size: 14px; font-weight: 600; color: {c['ink']};
             margin: 16px 0 6px; }}
a {{ color: {c['accent']}; }}
#simnibs_logo {{ max-height: 40px; width: auto; }}
.guicontrols {{ color: {c['ink2']}; margin-top: 8px; }}
select, input, button {{ font-family: {FONT}; font-size: 13px; color: {c['ink']};
                         background: {c['surface2']}; border: 1px solid {c['line2']};
                         border-radius: 4px; padding: 2px 4px; }}
input[type=range] {{ accent-color: {c['accent']}; padding: 0; border: none; background: none; }}
"""


def style_js(c: dict) -> str:
    """JavaScript that puts report_css(c) into the page (replacing an earlier one)."""
    return ("(function(){var s=document.getElementById(%s);if(!s){s=document.createElement('style');"
            "s.id=%s;(document.head||document.documentElement).appendChild(s);}s.textContent=%s;})();"
            % (json.dumps(STYLE_ID), json.dumps(STYLE_ID), json.dumps(report_css(c))))


def report_path(model) -> Path | None:
    """charm_report.html of this subject's segmentation (layered medium), if there is one."""
    if str(model.get("simulation.medium") or "layered") != "layered":
        return None
    m2m = sg.m2m_path(model)
    if m2m is None:
        return None
    r = m2m / "charm_report.html"
    return r if r.is_file() else None


class CharmReportCard(QFrame):
    def __init__(self, model, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.report: Path | None = None
        self.web = None
        self._no_web = False
        self._loaded: tuple | None = None       # (path, mtime) in the view
        v = QVBoxLayout(self)
        v.setContentsMargins(16, 12, 16, 14)          # as the section cards (section.py)
        head = QHBoxLayout()
        self.icon = QLabel()
        self.icon.setFixedSize(22, 22)
        self.retint()
        head.addWidget(self.icon)
        title = QLabel("SimNIBS report")
        title.setObjectName("CardTitle")
        head.addWidget(title)
        head.addStretch(1)
        v.addLayout(head)
        self.sub = QLabel()
        self.sub.setObjectName("CardSub")
        self.sub.setWordWrap(True)
        v.addWidget(self.sub)
        self.button = QToolButton(text="Show report")
        self.button.setCheckable(True)
        self.button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.button.setArrowType(Qt.ArrowType.RightArrow)
        self.button.setAutoRaise(True)
        self.button.toggled.connect(self._toggled)
        v.addWidget(self.button, 0, Qt.AlignmentFlag.AlignLeft)
        self.body = QWidget()
        self._body = QVBoxLayout(self.body)
        self._body.setContentsMargins(0, 4, 0, 0)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        self.note.hide()
        self._body.addWidget(self.note)
        self.body.hide()
        v.addWidget(self.body)
        self.refresh()

    def retint(self) -> None:
        """The card icon, and the report's look, in the active palette (as SectionPage.retint)."""
        self.icon.setPixmap(icons.icon("report", style.active()["accent"], 20).pixmap(20, 20))
        if self.web is not None:
            self._restyle()

    def _restyle(self) -> None:
        """Style the page on every load (DocumentReady script) and the one shown now."""
        from PyQt6.QtWebEngineCore import QWebEngineScript
        c = style.active()
        page = self.web.page()
        page.setBackgroundColor(QColor(c["surface"]))   # no white flash before the report draws
        scripts = page.scripts()
        for old in scripts.find(STYLE_ID):
            scripts.remove(old)
        script = QWebEngineScript()
        script.setName(STYLE_ID)
        script.setInjectionPoint(QWebEngineScript.InjectionPoint.DocumentReady)
        script.setWorldId(QWebEngineScript.ScriptWorldId.MainWorld)
        script.setSourceCode(style_js(c))
        scripts.insert(script)
        page.runJavaScript(style_js(c))

    def refresh(self) -> None:
        self.report = report_path(self.model)
        self.setVisible(self.report is not None)
        if self.report is None:
            return
        self.sub.setText(f"charm's quality check of {self.report.parent.name}: the final segmentation, the "
                         "atlas registration, and the T1–T2 registration.")
        if self.button.isChecked():
            self._load()

    def _toggled(self, on: bool) -> None:
        self.button.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        self.button.setText("Hide report" if on else "Show report")
        self.body.setVisible(on)
        if on:
            self._load()

    def _view(self):
        """The web view, built on first use (QtWebEngine is imported before the QApplication
        in app.py). None without PyQt6-WebEngine."""
        if self.web is None and not self._no_web:
            try:
                from PyQt6.QtWebEngineWidgets import QWebEngineView
            except ImportError:
                self._no_web = True
                return None
            self.web = QWebEngineView()
            self.web.setFixedHeight(VIEW_HEIGHT)
            self._body.addWidget(self.web)
            self._restyle()
        return self.web

    def _load(self) -> None:
        if self.report is None:
            return
        web = self._view()
        if web is None:
            self.note.setText("PyQt6-WebEngine is not installed: use 'Open charm report' (SimNIBS "
                              "segmentation card) to see it in the browser.")
            self.note.show()
            return
        try:
            key = (str(self.report), self.report.stat().st_mtime)
        except OSError:
            key = None
        if key != self._loaded:
            self._loaded = key
            web.load(QUrl.fromLocalFile(str(self.report)))
