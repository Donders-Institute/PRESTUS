"""Main window: a header, one tab per area (setup, transducer, placement, simulation,
study, run), the issues panel and the review."""

from __future__ import annotations

import shutil
from pathlib import Path

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QAction, QColor, QIcon, QKeySequence, QPainter, QPixmap
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFileDialog, QFrame, QHBoxLayout,
                             QLabel, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea,
                             QSplitter, QTabWidget, QTreeWidget, QTreeWidgetItem,
                             QVBoxLayout, QWidget)

from .. import checks
from ..model import TX, ConfigModel, is_empty, iter_leaves
from ..paths import DEFAULT_NAME, base_default_for, toolbox_default
from ..project_defaults import carry_pinned_values, replace_project_default
from ..schema import GENERAL, Schema
from ..session import Session
from . import style
from .placement_viewer import QC_PATH
from .section import SectionPage

REVIEW = "__review__"
VIEWER = "__viewer__"
RUN = "__run__"
STUDY = "__study__"
RESULTS = "__results__"
PCT = "__pct__"
CALIBRATION = "__calibration__"
LOCALITE = "__localite__"
START = "__start__"
SPECIAL = (REVIEW, VIEWER, RUN, STUDY, RESULTS, PCT, CALIBRATION, LOCALITE, START)
SUBTAB_PAGES = (VIEWER, PCT, RESULTS, CALIBRATION, LOCALITE)
LEVEL_ICON = {"error": "✗", "warning": "⚠", "info": "ℹ"}
LEVEL_COLOR_KEY = {"error": "red", "warning": "amber", "info": "info"}
_OPEN_WINDOWS: list = []          # keeps windows opened from File > New / Open alive


def _dot_icon(colour: str, size: int = 10) -> QIcon:
    """A filled dot, for status on tabs (the tab style sheet would override a text colour)."""
    pix = QPixmap(size, size)
    pix.fill(Qt.GlobalColor.transparent)
    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setBrush(QColor(colour))
    p.setPen(Qt.PenStyle.NoPen)
    p.drawEllipse(0, 0, size, size)
    p.end()
    return QIcon(pix)


def study_header(path: Path) -> str:
    return (f"PRESTUS study config: {path.stem}\n"
            f"Written by the PRESTUS Nexus; only values that differ from config_default.yaml.\n"
            f"Load with: parameters = load_parameters('{path.name}', '{path.parent}');")


class ReviewPage(QWidget):
    def __init__(self, window: "MainWindow"):
        super().__init__()
        self.window = window
        v = QVBoxLayout(self)
        title = QLabel("Review & save")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.info)
        self.text = QPlainTextEdit(readOnly=True)
        mono = style.mono_font()
        self.text.setFont(mono)
        v.addWidget(self.text, 1)
        row = QHBoxLayout()
        prune = QPushButton("Drop values equal to defaults")
        prune.setToolTip("Remove keys that only repeat config_default.yaml")
        prune.clicked.connect(window.prune_defaults)
        row.addWidget(prune)
        save = QPushButton("Save")
        save.setProperty("primary", True)
        save.clicked.connect(window.save)
        save_as = QPushButton("Save as…")
        save_as.clicked.connect(window.save_as)
        row.addStretch(1)
        row.addWidget(save_as)
        row.addWidget(save)
        v.addLayout(row)

    def refresh(self) -> None:
        m = self.window.model
        path = m.path or Path("config_study.yaml")
        self.text.setPlainText(m.to_yaml(study_header(path)))
        n_err = sum(1 for i in self.window.issues if i.level == "error")
        n_warn = sum(1 for i in self.window.issues if i.level == "warning")
        self.info.setText(
            f"<b>File:</b> {m.path or '(not saved)'}<br>"
            + (f"<b>Merged onto:</b> {m.schema.toolbox_path} + project defaults "
               f"{m.schema.project_path}<br>" if m.schema.project_path is not None
               else f"<b>Merged onto:</b> {m.schema.default_path}<br>") +
            f"<b>{m.modified_count()}</b> value(s) differ from the defaults · "
            f"{n_err} error(s), {n_warn} warning(s) — see Issues in the header.")


class TabPage(QWidget):
    """One tab: intro, 'Show advanced options', section cards (+ extras), optional side panel."""

    def __init__(self, window: "MainWindow", spec: dict):
        super().__init__()
        self.window = window
        self.key = spec["key"]
        self.sections: list[SectionPage] = []
        self.cards: list[QWidget] = []
        self.side: QWidget | None = None
        self.subtabs: QTabWidget | None = None
        self.settings_page: QWidget | None = None     # the 'Settings' sub-tab, once there are sub-tabs

        outer = QVBoxLayout(self)
        outer.setContentsMargins(14, 10, 14, 6)
        head = QHBoxLayout()
        intro = QLabel(spec.get("intro", ""))
        intro.setObjectName("TabIntro")
        intro.setWordWrap(True)
        head.addWidget(intro, 1)
        self.hidden = QLabel()
        self.hidden.setObjectName("CardSub")
        head.addWidget(self.hidden)
        self.advanced = QCheckBox("Show advanced options")
        self.advanced.setToolTip("Show the settings most studies leave at their defaults. "
                                 "Advanced settings that hold a value are always shown.")
        self.advanced.setChecked(bool(window.session.get(f"advanced_{self.key}", False)))
        self.advanced.toggled.connect(self._toggle_advanced)
        head.addWidget(self.advanced)
        self.tune_label = QLabel("Tune for:")
        self.tune = QComboBox()
        self.tune.addItem("—", None)
        schema = window.model.schema
        for goal, spec in schema.tune_goals.items():
            self.tune.addItem(schema.tune_label(goal), goal)
            self.tune.setItemData(self.tune.count() - 1, spec.get("help", ""),
                                  Qt.ItemDataRole.ToolTipRole)
        self.tune.setToolTip("For optimising a run: show the settings that affect speed, accuracy, "
                             "memory or disk space (advanced ones included) and mark them. "
                             "Applies to all tabs.")
        self.tune_label.setToolTip(self.tune.toolTip())
        self.tune.currentIndexChanged.connect(
            lambda i: self.window.set_tune_goal(self.tune.itemData(i)))
        head.addWidget(self.tune_label)
        head.addWidget(self.tune)
        outer.addLayout(head)

        self.split = QSplitter()
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        body = QWidget()
        self.column = QVBoxLayout(body)
        self.column.setContentsMargins(0, 4, 6, 12)
        self.column.setSpacing(12)
        self.scroll.setWidget(body)
        self.split.addWidget(self.scroll)
        outer.addWidget(self.split, 1)

    def add_card(self, widget: QWidget) -> None:
        self.cards.append(widget)
        self.column.addWidget(widget)

    def add_section(self, page: SectionPage) -> None:
        page.set_show_advanced(self.advanced.isChecked())
        page.advanced_toggled.connect(self.refresh)
        self.sections.append(page)
        self.column.addWidget(page)

    def set_side(self, widget: QWidget, stretch: tuple[int, int] = (3, 2)) -> None:
        self.side = widget
        self.scroll.setMinimumWidth(560)
        self.split.addWidget(widget)
        self.split.setStretchFactor(0, stretch[0])
        self.split.setStretchFactor(1, stretch[1])
        self.split.setChildrenCollapsible(False)
        self.split.setSizes([1000 * stretch[0], 1000 * stretch[1]])

    def set_full_tab(self, widget: QWidget, label: str) -> None:
        """Give `widget` its own full-width sub-tab, next to the settings cards as 'Settings'
        (the same nested-tab pattern the Run tab already uses for Run log / Review YAML), rather
        than squeezing it into a splitter half the width."""
        self.side = widget
        outer = self.layout()
        outer.removeWidget(self.split)
        self.subtabs = QTabWidget()
        self.subtabs.setDocumentMode(True)
        self.subtabs.addTab(self.scroll, "Settings")
        self.settings_page = self.scroll
        self.subtabs.addTab(widget, label)
        self.subtabs.currentChanged.connect(self._sub_tab_changed)
        outer.addWidget(self.subtabs, 1)
        self.split.setParent(None)

    def add_subtab(self, widget: QWidget, label: str) -> None:
        """Another full-width sub-tab. A tab without sub-tabs gets them here, with its
        settings (and side panel, if any) as 'Settings'."""
        if self.subtabs is None:
            outer = self.layout()
            outer.removeWidget(self.split)
            self.subtabs = QTabWidget()
            self.subtabs.setDocumentMode(True)
            self.subtabs.addTab(self.split, "Settings")
            self.settings_page = self.split
            self.subtabs.currentChanged.connect(self._sub_tab_changed)
            outer.addWidget(self.subtabs, 1)
        self.subtabs.addTab(widget, label)

    def _sub_tab_changed(self, _index: int) -> None:
        w = self.subtabs.widget(_index)
        if hasattr(w, "refresh"):
            w.refresh()

    def finish(self) -> None:
        self.column.addStretch(1)
        tunable = any(r.node.tune for page in self.sections for r in page.rows)
        self.tune_label.setVisible(tunable)
        self.tune.setVisible(tunable)

    def set_goal(self, goal: str | None) -> None:
        self.tune.blockSignals(True)
        self.tune.setCurrentIndex(max(0, self.tune.findData(goal)))
        self.tune.blockSignals(False)
        for page in self.sections:
            page.set_goal(goal)
        self.refresh()

    def _toggle_advanced(self, on: bool) -> None:
        self.window.session.set(f"advanced_{self.key}", on)
        for page in self.sections:
            page.set_show_advanced(on)
        self.refresh()

    def refresh(self) -> None:
        show = self.advanced.isChecked()
        schema = self.window.model.schema
        n_hidden = 0
        for page in self.sections:
            page.refresh()
            sec = page.section
            visible = sec == GENERAL or self.window.model.is_visible(sec, 0)
            if schema.section_is_advanced(sec) and not show and not page.has_overrides() \
                    and not page.matches_goal():
                n_hidden += 1
                visible = False
            page.setVisible(visible and any(r.block.isVisibleTo(page) for r in page.rows)
                            or (visible and page.section == TX))
            n_hidden += page.hidden_advanced() if visible else 0
        for card in self.cards:
            if hasattr(card, "refresh"):
                card.refresh()
        self.hidden.setText(f"{n_hidden} advanced hidden" if n_hidden and not show else "")


class MainWindow(QMainWindow):
    def __init__(self, model: ConfigModel, session: Session, prestus_root: Path):
        super().__init__()
        self.model = model
        self.session = session
        self.prestus_root = prestus_root
        self.dirty = False
        self.issues: list = []
        self.pages: dict[str, QWidget] = {}
        self.tab_pages: dict[str, TabPage] = {}
        self._section_tab: dict[str, str] = {}
        self.tune_goal: str | None = None    # 'Tune for' goal, shared by all tabs (not remembered)
        self.resize(1320, 900)

        central = QWidget()
        cv = QVBoxLayout(central)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(0)
        cv.addWidget(self._header())
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        cv.addWidget(self.tabs, 1)
        self.setCentralWidget(central)

        self.issue_popup = self._issue_popup()

        self._build_menu()
        self._build_view_menu()
        self._build_tools_menu()
        self._build_tabs()
        self._refresh_cal_cov()
        if hasattr(self, "io_preview"):
            self.io_preview.refresh()
        self.tabs.currentChanged.connect(self._tab_changed)
        self.tabs.setCurrentIndex(0)

        self._issue_timer = QTimer(self, singleShot=True, interval=200)
        self._issue_timer.timeout.connect(self._refresh_issues)
        self._refresh_issues()
        self._update_title()

    # ── construction ────────────────────────────────────────────────────
    def _header(self) -> QWidget:
        bar = QFrame()
        bar.setObjectName("HeaderBar")
        h = QHBoxLayout(bar)
        h.setContentsMargins(16, 8, 16, 8)
        pix = style.logo_pixmap(self.prestus_root, 36)
        if pix is not None:
            logo = QLabel()
            logo.setPixmap(pix)
            logo.setStyleSheet(f"background: {style.C['surface']}; padding: 2px; border-radius: 4px;")
            h.addWidget(logo)
        texts = QVBoxLayout()
        texts.setSpacing(0)
        title = QLabel("PRESTUS Nexus")
        title.setObjectName("HeaderTitle")
        texts.addWidget(title)
        self.header_file = QLabel()
        self.header_file.setObjectName("HeaderSub")
        texts.addWidget(self.header_file)
        h.addLayout(texts)
        h.addStretch(1)
        self.header_issues = QPushButton()
        self.header_issues.setToolTip("Show the issues (Ctrl+Shift+I)")
        self.header_issues.clicked.connect(self.show_issues)
        h.addWidget(self.header_issues)
        open_btn = QPushButton("Open…")
        open_btn.setToolTip("Open a study config (File ▸ Open…)")
        open_btn.clicked.connect(lambda: self.open_study())
        h.addWidget(open_btn)
        save = QPushButton("Save")
        save.setProperty("primary", True)
        save.clicked.connect(self.save)
        h.addWidget(save)
        return bar

    def _issue_popup(self) -> QFrame:
        """The issues list, opened from the header counts. A popup rather than a dock: it
        takes no room until asked for and opens at a size set by the window, not by
        whatever space the tabs leave over. Clicking outside or Esc closes it."""
        popup = QFrame(self, Qt.WindowType.Popup)
        popup.setFrameShape(QFrame.Shape.StyledPanel)
        v = QVBoxLayout(popup)
        v.setContentsMargins(8, 8, 8, 8)
        self.issue_counts_lbl = QLabel()
        v.addWidget(self.issue_counts_lbl)
        self.issue_tree = QTreeWidget()
        self.issue_tree.setHeaderLabels(["", "Setting", "Problem"])
        self.issue_tree.setRootIsDecorated(False)
        self.issue_tree.setWordWrap(True)
        self.issue_tree.setColumnWidth(0, 28)
        self.issue_tree.setColumnWidth(1, 300)
        self.issue_tree.itemActivated.connect(self._goto_issue)
        self.issue_tree.itemClicked.connect(self._goto_issue)
        v.addWidget(self.issue_tree, 1)
        hint = QLabel("Click an issue to open its field.")
        hint.setObjectName("HeaderSub")
        v.addWidget(hint)
        return popup

    def show_issues(self) -> None:
        """Open the issues popup under the header counts, right-aligned with them, at up
        to 70 % of the window width and 60 % of its height (less if the list is short)."""
        self._refresh_issues()
        rows = max(self.issue_tree.topLevelItemCount(), 1)
        row_h = max(self.issue_tree.sizeHintForRow(0), self.fontMetrics().height() + 6)
        want_h = 90 + rows * row_h
        w = max(min(int(self.width() * 0.7), 1100), 560)
        h = max(min(want_h, int(self.height() * 0.6)), 200)
        self.issue_popup.resize(w, h)
        anchor = self.header_issues.mapToGlobal(self.header_issues.rect().bottomRight())
        screen = self.screen().availableGeometry()
        x = max(min(anchor.x() - w, screen.right() - w), screen.left())
        y = min(anchor.y() + 4, screen.bottom() - h)
        self.issue_popup.move(x, y)
        self.issue_popup.show()
        self.issue_tree.setFocus()

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        for text, slot, key in [("New study…", self.new_study, QKeySequence.StandardKey.New),
                                ("Open…", lambda: self.open_study(), QKeySequence.StandardKey.Open)]:
            action = QAction(text, self)
            action.triggered.connect(slot)
            action.setShortcut(key)
            file_menu.addAction(action)
        self.recent_menu = file_menu.addMenu("Open recent")
        self.recent_menu.aboutToShow.connect(self._fill_recent)
        file_menu.addSeparator()
        for text, slot, key in [("Save", self.save, QKeySequence.StandardKey.Save),
                                ("Save as…", self.save_as, QKeySequence.StandardKey.SaveAs),
                                ("Reload from disk", self.reload, None),
                                ("Quit", self.close, QKeySequence.StandardKey.Quit)]:
            action = QAction(text, self)
            action.triggered.connect(slot)
            if key is not None:
                action.setShortcut(key)
            file_menu.addAction(action)

    def _build_tools_menu(self) -> None:
        tools = self.menuBar().addMenu("&Tools")
        for text, slot in [("Localite ingestion", lambda: self.show_special(LOCALITE)),
                           ("Transducer calibration", lambda: self.show_special(CALIBRATION)),
                           (None, None),
                           ("neuronav_gui in MATLAB…", self.open_neuronav_gui),
                           ("calibration_gui in MATLAB…", self.open_calibration_gui),
                           (None, None),
                           ("Update project defaults…", self.update_project_defaults)]:
            if text is None:
                tools.addSeparator()
                continue
            action = QAction(text, self)
            action.triggered.connect(slot)
            tools.addAction(action)

    def _build_view_menu(self) -> None:
        view = self.menuBar().addMenu("&View")
        self.dark_mode_action = QAction("Dark mode", self, checkable=True)
        self.dark_mode_action.setChecked(style.is_dark())
        self.dark_mode_action.toggled.connect(self._set_dark_mode)
        view.addAction(self.dark_mode_action)
        issues = QAction("Issues…", self)
        issues.setShortcut(QKeySequence("Ctrl+Shift+I"))
        issues.triggered.connect(self.show_issues)
        view.addAction(issues)

    def _set_dark_mode(self, dark: bool) -> None:
        style.apply(QApplication.instance(), self.prestus_root, dark=dark)
        self.session.set("dark_mode", dark)
        self.session.save()

    def _tab_specs(self) -> list[dict]:
        schema = self.model.schema
        specs = [dict(t) for t in schema.tabs] or [{"key": "setup", "label": "Setup", "sections": []}]
        listed = {s for t in specs for s in t.get("sections", [])}
        rest = [s for s in schema.section_order() if s not in listed]
        specs[0]["sections"] = list(specs[0].get("sections", [])) + rest
        return specs

    def _build_tabs(self) -> None:
        schema = self.model.schema
        from .start_page import StartPage
        start = StartPage(self)
        start.goto_tab.connect(self._start_goto)
        self.pages[START] = start
        self.tabs.addTab(start, "Start")
        self.tabs.setTabToolTip(0, "A first simulation in four steps; everything else keeps its default")
        for spec in self._tab_specs():
            tab = TabPage(self, spec)
            key = spec["key"]
            if key == "study":
                from .study_page import StudyPage
                study = StudyPage(self)
                study.changed.connect(self._changed)
                self.pages[STUDY] = study
                self._section_tab[STUDY] = key
            for section in spec.get("sections", []):
                if section != GENERAL and schema.get(section) is None:
                    continue
                nodes = schema.display_nodes(section)
                if not nodes:
                    continue
                page = SectionPage(self.model, section, nodes, schema.section_label(section))
                page.changed.connect(self._changed)
                self.pages[section] = page
                self._section_tab[section] = key
                if section == "placement":
                    page.add_extra(LocaliteActions(self))
                if section == "timing":
                    from .timing_view import TimingPanel
                    self.timing_panel = TimingPanel(self.model)
                    self.timing_panel.changed.connect(self._changed)
                    page.add_extra(self.timing_panel)
                if section == "hpc":
                    from .run_profile import RunProfilePanel
                    self.run_profile_panel = RunProfilePanel(self)
                    self.run_profile_panel.changed.connect(self._changed)
                    self.run_profile_panel.changed.connect(page.refresh)
                    page.add_extra(self.run_profile_panel)
                tab.add_section(page)
                if section == "study":
                    tab.add_card(self.pages[STUDY])
            if key == "setup":
                from .io_preview import IOPreviewPanel
                self.io_preview = IOPreviewPanel(self.model)
                tab.set_side(self.io_preview, (3, 2))
            elif key == "transducer":
                from .transducer_view import TransducerSchematic
                self.schematic = TransducerSchematic(self.model)
                card = self._card("Transducer geometry", self.schematic,
                                  "Cross-section along the beam: elements, exit plane, geometric focus, the "
                                  "focal depth setting (✕) and the device's steering range.")
                from .calibration_plots import CalibrationCoverage
                self.cal_cov = CalibrationCoverage(height=220)
                cov_card = self._card("Calibration coverage", self.cal_cov,
                                      "What the library entry covers, from its YAML: fit quality, drive and "
                                      "phase corrections per calibrated depth; this study's depths as lines.")
                # Both cards at their natural height in one scrolling column, so a short
                # window scrolls instead of squeezing the plot under its caption
                column = QWidget()
                cl = QVBoxLayout(column)
                cl.setContentsMargins(0, 0, 4, 0)
                cl.setSpacing(12)
                card.setMinimumHeight(300)
                cl.addWidget(card, 1)
                cl.addWidget(cov_card)
                side = QScrollArea()
                side.setWidgetResizable(True)
                side.setWidget(column)
                tab.set_side(side)
                from .calibration_page import CalibrationPage
                calib = CalibrationPage(self)
                self.pages[CALIBRATION] = calib
                self._section_tab[CALIBRATION] = key
                tab.add_subtab(calib, "Calibration")
            elif key == "placement":
                from .placement_viewer import PlacementPage
                viewer = PlacementPage(self.model)
                viewer.changed.connect(self._changed)
                viewer.navigate.connect(self.goto)          # empty state: 'Set T1 image…' etc.
                self.pages[VIEWER] = viewer
                self._section_tab[VIEWER] = key
                tab.set_full_tab(viewer, "Viewer")
                from .localite_page import LocalitePage
                loc = LocalitePage(self)
                self.pages[LOCALITE] = loc
                self._section_tab[LOCALITE] = key
                tab.add_subtab(loc, "Localite")
            elif key == "simulation":
                from .pct_page import PctPage
                pct_page = PctPage(self)
                self.pages[PCT] = pct_page
                self._section_tab[PCT] = key
                tab.set_full_tab(pct_page, "Pseudo-CT")
                from .uncertainty_view import UncertaintyCard
                self.uncertainty_card = UncertaintyCard(self.model, self.prestus_root)
                tab.add_card(self.uncertainty_card)
            elif key == "run":
                from .run_page import RunPage
                side = QTabWidget()
                run = RunPage(self)
                review = ReviewPage(self)
                self.pages[RUN], self.pages[REVIEW] = run, review
                self._section_tab[RUN] = self._section_tab[REVIEW] = key
                side.addTab(run, "Run log")
                side.addTab(review, "Review YAML")
                side.addTab(run.jobs_view, "Jobs")
                side.currentChanged.connect(lambda _i, s=side: s.currentWidget().refresh())
                self.run_side = side
                tab.set_side(side, (2, 3))
            elif key == "results":
                from .results_page import ResultsPage
                results = ResultsPage(self)
                self.pages[RESULTS] = results
                self._section_tab[RESULTS] = key
                tab.set_full_tab(results, "Results")
                tab.subtabs.blockSignals(True)          # scanned when the tab is first shown
                tab.subtabs.setCurrentWidget(results)   # skip the empty 'Settings' pane
                tab.subtabs.blockSignals(False)
            tab.finish()
            self.tab_pages[key] = tab
            self.tabs.addTab(tab, spec.get("label", key))
        if RUN in self.pages and RESULTS in self.pages:
            self.pages[RUN].finished.connect(
                lambda _task, ok: self.pages[RESULTS].refresh() if ok else None)

    def _refresh_cal_cov(self) -> None:
        """Calibration coverage of the transducer being edited (library YAML only, no MATLAB)."""
        if not hasattr(self, "cal_cov"):
            return
        from .. import calibration as cal
        m = self.model
        tx = self.pages.get(TX)
        ti = tx.ti if isinstance(tx, SectionPage) else 0
        serial, ds, lib = (m.get(f"{TX}.serial", ti), m.get(f"{TX}.combo.ds_serial", ti),
                           m.get(f"{TX}.library_path", ti))
        serial = "" if is_empty(serial) else str(serial)
        ds = "" if is_empty(ds) else str(ds)
        lib = "" if is_empty(lib) else str(lib)
        depths = cal.study_depths(m, serial, ds) if serial else []
        self.cal_cov.set_device(m.library, serial, ds, lib, depths)

    @staticmethod
    def _card(title: str, widget: QWidget, sub: str = "") -> QFrame:
        card = QFrame()
        card.setObjectName("Card")
        v = QVBoxLayout(card)
        v.setContentsMargins(14, 10, 14, 12)
        t = QLabel(title)
        t.setObjectName("CardTitle")
        v.addWidget(t)
        if sub:
            s = QLabel(sub)
            s.setObjectName("CardSub")
            s.setWordWrap(True)
            v.addWidget(s)
        v.addWidget(widget, 1)
        return card

    # ── MATLAB GUIs ─────────────────────────────────────────────────────
    def _matlab(self) -> str | None:
        from .. import launcher
        run = self.pages.get(RUN)
        path = (run.matlab.text() if run is not None else "") or launcher.find_matlab(
            self.session.get("matlab_path"))
        if not path:
            QMessageBox.warning(self, "MATLAB", "MATLAB not found; set its path on the Run tab.")
        return path or None

    def matlab_path(self) -> str | None:
        """MATLAB binary without asking: the Run tab's field, the session, or a standard install."""
        from .. import launcher
        run = self.pages.get(RUN)
        path = (run.matlab.text() if run is not None else "") or launcher.find_matlab(
            self.session.get("matlab_path"))
        return path or None

    def bridge(self):
        """The shared MATLAB session for interactive steps (Localite, Calibration)."""
        from .. import bridge
        return bridge.shared_bridge(self.prestus_root, self.matlab_path)

    def _start_matlab_gui(self, call: str, what: str) -> None:
        from PyQt6.QtCore import QProcess

        from .. import launcher
        matlab = self._matlab()
        if matlab is None:
            return
        argv = launcher.gui_command(matlab, self.prestus_root, call)
        ok, _pid = QProcess.startDetached(argv[0], argv[1:])
        self.statusBar().showMessage(
            f"Starting MATLAB with {what}…" if ok else f"Could not start {matlab}", 8000)

    def open_neuronav_gui(self) -> None:
        """Localite ingestion, pre-filled from the saved study config."""
        from .. import launcher
        if (self.dirty or self.model.path is None) and not self.save():
            return
        self._start_matlab_gui(launcher.neuronav_gui_call(self.model.path), "neuronav_gui")

    def open_calibration_gui(self) -> None:
        """Library fitting; entries go to transducer.library_path when it is set."""
        from .. import launcher
        lib = self.model.get("transducer.library_path")
        if not lib:
            QMessageBox.information(
                self, "Calibration GUI",
                "Tip: set transducer.library_path to a folder of your own; calibration_gui then "
                "saves there and the Nexus's calibration coverage picks the entries up.")
        self._start_matlab_gui(launcher.calibration_gui_call(lib or None), "calibration_gui")

    # ── navigation ──────────────────────────────────────────────────────
    def _page(self, key: str) -> QWidget:
        return self.pages[key]

    def _start_goto(self, key: str) -> None:
        if key == "run":
            self.show_special(RUN)
        else:
            self._select_tab(key)

    def on_model_edited(self) -> None:
        """A page outside the form tabs (Start) changed the model: mark it modified and
        check it now rather than after the debounce, so the page's status follows at once."""
        self._changed()
        self._issue_timer.stop()
        self._refresh_issues()

    def _tab_changed(self, _index: int) -> None:
        if self.tabs.currentWidget() is self.pages.get(START):
            self.pages[START].refresh()
            return
        tab = self.tabs.currentWidget()
        if isinstance(tab, TabPage):
            tab.refresh()
            if tab.key == "placement" and VIEWER in self.pages:
                self.pages[VIEWER].refresh()
                self._issue_timer.start()           # the viewer's checks, now its images are loaded
            if tab.key == "simulation" and PCT in self.pages:
                self.pages[PCT].refresh()
            if tab.key == "run":
                if VIEWER in self.pages and self.model.get("placement.mode") in ("manual", None, ""):
                    QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                    try:                            # check the placement before a run
                        self.pages[VIEWER].ensure_images()
                    finally:
                        QApplication.restoreOverrideCursor()
                    self._refresh_issues()
                self.run_side.currentWidget().refresh()
            if tab.key == "results" and RESULTS in self.pages:
                self.pages[RESULTS].refresh()
            if tab.subtabs is not None and tab.subtabs.currentWidget() in (
                    self.pages.get(CALIBRATION), self.pages.get(LOCALITE)):
                tab.subtabs.currentWidget().refresh()

    def _select_tab(self, key: str) -> TabPage | None:
        tab = self.tab_pages.get(key)
        if tab is not None:
            self.tabs.setCurrentWidget(tab)
        return tab

    def show_special(self, key: str) -> QWidget:
        """Show a special page (Run, Review, Study manifest, viewer) and return it."""
        tab = self._select_tab(self._section_tab.get(key, ""))
        page = self.pages[key]
        if key in (RUN, REVIEW):
            self.run_side.setCurrentWidget(page)
        elif key in SUBTAB_PAGES and tab is not None and tab.subtabs is not None:
            tab.subtabs.setCurrentWidget(page)
        page.refresh()
        return page

    def display_section(self, path: str) -> str:
        schema = self.model.schema
        for node in reversed(schema.ancestors(path)):
            if node.group:
                return node.group
        section = path.split(".", 1)[0] if "." in path else GENERAL
        node = schema.get(section)
        if node is None or not node.is_group:
            section = GENERAL
        return section

    def set_tune_goal(self, goal: str | None) -> None:
        """Show (and mark) the settings tagged with `goal` in overlay.yaml 'tuning', on every tab."""
        self.tune_goal = goal
        for tab in self.tab_pages.values():
            tab.set_goal(goal)

    def goto(self, path: str, ti: int | None = None) -> None:
        section = self.display_section(path)
        page = self.pages.get(section)
        if not isinstance(page, SectionPage):
            return
        tab = self._select_tab(self._section_tab.get(section, ""))
        if tab is not None:
            if tab.subtabs is not None:
                tab.subtabs.setCurrentWidget(tab.settings_page)   # the field lives under 'Settings'
            tab.refresh()
            row = next((r for r in page.rows if r.node.path == path), None)
            if not page.isVisibleTo(tab):
                tab.advanced.setChecked(True)       # the target is in an advanced section
            elif row is not None and not row.block.isVisibleTo(page):
                page.set_show_advanced(True)        # an advanced field: open only its card
                tab.refresh()
        page.focus_key(path, ti)

    def _goto_issue(self, item: QTreeWidgetItem, _col: int = 0) -> None:
        issue = item.data(0, Qt.ItemDataRole.UserRole)
        if issue is not None:
            self.goto_issue(issue)

    def goto_issue(self, issue) -> None:
        """Open the field (or the placement viewer) an issue points to."""
        self.issue_popup.hide()
        if issue.path == QC_PATH:
            self.show_special(VIEWER).show_transducer(issue.transducer or 0)
        else:
            self.goto(issue.path, issue.transducer)

    # ── state ───────────────────────────────────────────────────────────
    def _section_count(self, section: str) -> int:
        m = self.model
        if section in SPECIAL:
            return 0
        if section == GENERAL:
            return sum(1 for k, v in m.overrides.items() if not isinstance(v, dict))
        if section == TX:
            return sum(1 for t in m.tx for _ in iter_leaves(t))
        value = m.overrides.get(section)
        return sum(1 for _ in iter_leaves(value)) if isinstance(value, dict) else int(value is not None)

    def _update_tab_counts(self) -> None:
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            if not isinstance(tab, TabPage):
                continue
            n = sum(self._section_count(p.section) for p in tab.sections)
            label = next((t.get("label", t["key"]) for t in self._tab_specs() if t["key"] == tab.key), tab.key)
            # The tab shows what needs attention (errors, else warnings) as a coloured dot and
            # count; the number of changed values, which needs none, moves to the tooltip.
            mine = [iss for iss in self.issues if self._issue_tab(iss) == tab.key]
            n_err = sum(1 for iss in mine if iss.level == "error")
            n_warn = sum(1 for iss in mine if iss.level == "warning")
            level, count = ("red", n_err) if n_err else ("amber", n_warn) if n_warn else (None, 0)
            self.tabs.setTabText(i, f"{label}  {count}" if count else label)
            self.tabs.setTabIcon(i, _dot_icon(style.active()[level]) if level else QIcon())
            tip = [f"{n_err} error(s), {n_warn} warning(s)" if n_err or n_warn else "No problems found"]
            if n:
                tip.append(f"{n} value(s) changed from the defaults")
            self.tabs.setTabToolTip(i, " · ".join(tip))

    def _issue_tab(self, issue) -> str | None:
        """Key of the tab that holds the field an issue points to."""
        if issue.path == QC_PATH:
            return self._section_tab.get(VIEWER)
        return self._section_tab.get(self.display_section(issue.path))

    def _changed(self) -> None:
        self.dirty = True
        tab = self.tabs.currentWidget()
        if isinstance(tab, TabPage):
            tab.refresh()
        if hasattr(self, "timing_panel"):
            self.timing_panel.refresh()
        if hasattr(self, "schematic"):
            tx = self.pages.get(TX)
            self.schematic.set_transducer(tx.ti if isinstance(tx, SectionPage) else 0)
            self._refresh_cal_cov()
        if hasattr(self, "io_preview"):
            self.io_preview.refresh()
        run = self.pages.get(RUN)
        if run is not None:
            run.update_deployment_summary()
        self._update_title()
        self._update_tab_counts()
        self._issue_timer.start()

    def _refresh_issues(self) -> None:
        self.issues = checks.all_issues(self.model)
        if VIEWER in self.pages:
            self.issues += self.pages[VIEWER].qc_issues()
        self.issue_tree.clear()
        palette = style.active()
        for issue in self.issues:
            item = QTreeWidgetItem([LEVEL_ICON.get(issue.level, ""), issue.where(), issue.message])
            item.setData(0, Qt.ItemDataRole.UserRole, issue)
            item.setToolTip(2, issue.message)
            key = LEVEL_COLOR_KEY.get(issue.level)
            if key:
                item.setForeground(0, QColor(palette[key]))
            self.issue_tree.addTopLevelItem(item)
        n_err = sum(1 for i in self.issues if i.level == "error")
        n_warn = sum(1 for i in self.issues if i.level == "warning")
        counts = f"{n_err} error(s) · {n_warn} warning(s)"
        level = "red" if n_err else "amber" if n_warn else "green"
        self.header_issues.setText(f"●  {counts}  ▾")       # QPushButton text is plain, not rich text
        self.header_issues.setStyleSheet(f"QPushButton {{ color: {style.C[level]}; font-weight: 600; }}")
        self.issue_counts_lbl.setText(
            f"<b>Issues</b>  <span style='color:{palette[level]}'>●</span> {counts}")
        self._update_tab_counts()
        review = self.pages.get(REVIEW)
        if review is not None and review.isVisible():
            review.refresh()
        start = self.pages.get(START)
        if start is not None and self.tabs.currentWidget() is start:
            start.refresh()

    def _update_title(self) -> None:
        name = self.model.path.name if self.model.path else "untitled"
        self.setWindowTitle(f"PRESTUS Nexus — {name}{' *' if self.dirty else ''}")
        self.header_file.setText(f"{self.model.path or 'Unsaved study config'}{'  (modified)' if self.dirty else ''}")

    # ── file actions ────────────────────────────────────────────────────
    def _fill_recent(self) -> None:
        self.recent_menu.clear()
        for path in self.session.recent():
            action = QAction(path, self)
            action.triggered.connect(lambda _=False, p=path: self.open_study(p))
            self.recent_menu.addAction(action)
        if not self.recent_menu.actions():
            self.recent_menu.addAction("(none)").setEnabled(False)

    def _discard_ok(self) -> bool:
        """True when there are no unsaved changes, or the user saved or discarded them."""
        if not self.dirty:
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes", "Save changes to the current study config first?",
            QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel)
        if answer == QMessageBox.StandardButton.Cancel:
            return False
        return answer == QMessageBox.StandardButton.Discard or self.save()

    def _switch_to(self, model: ConfigModel) -> "MainWindow":
        """Replace this window by one for another study config (its defaults may differ)."""
        window = MainWindow(model, self.session, self.prestus_root)
        window.resize(self.size())
        _OPEN_WINDOWS.append(window)
        window.show()
        if model.path is not None:
            self.session.add_recent(model.path)
            self.session.save()
        self.dirty = False
        self.close()
        return window

    def open_study(self, path: str | Path | None = None) -> "MainWindow | None":
        if not self._discard_ok():
            return None
        if path is None:
            start = str(self.model.path.parent if self.model.path else Path.home())
            path, _ = QFileDialog.getOpenFileName(self, "Open study config", start, "YAML (*.yaml *.yml)")
            if not path:
                return None
        path = Path(path).resolve()
        model = ConfigModel.from_file(Schema(base_default_for(path, self.prestus_root), toolbox_path=toolbox_default(self.prestus_root)), path,
                                      self.model.library)
        return self._switch_to(model)

    def new_study(self, target: str | Path | None = None, copy_default: bool | None = None) -> "MainWindow | None":
        """Create a study config; optionally pin config_default.yaml in its folder (prestus_config_init)."""
        if not self._discard_ok():
            return None
        if target is None:
            start = str((self.model.path.parent if self.model.path else Path.home()) / "config_study.yaml")
            target, _ = QFileDialog.getSaveFileName(self, "New study config", start, "YAML (*.yaml)")
            if not target:
                return None
        target = Path(target).resolve()
        local = target.parent / DEFAULT_NAME
        if not local.exists():
            if copy_default is None:
                copy_default = QMessageBox.question(
                    self, "Project defaults",
                    f"Copy {DEFAULT_NAME} into {target.parent}?\n\nThis pins the defaults for this study, "
                    "as prestus_config_init does; load_parameters then merges onto that copy.") \
                    == QMessageBox.StandardButton.Yes
            if copy_default:
                shutil.copy2(toolbox_default(self.prestus_root), local)
        model = ConfigModel(Schema(base_default_for(target, self.prestus_root), toolbox_path=toolbox_default(self.prestus_root)), self.model.library)
        model.save(target, header=study_header(target))
        return self._switch_to(model)

    def update_project_defaults(self) -> "MainWindow | None":
        """Refresh the project's config_default.yaml from the toolbox (prestus_config_upgrade)."""
        schema = self.model.schema
        if schema.project_path is None:
            QMessageBox.information(self, "Project defaults",
                                    "This study is merged onto the toolbox config_default.yaml directly; "
                                    "there is no project copy to update.")
            return None
        if not schema.stale_keys and not schema.pinned:
            QMessageBox.information(self, "Project defaults",
                                    f"{schema.project_path} matches the toolbox default.")
            return None
        if self.model.path is None:
            QMessageBox.information(self, "Project defaults", "Save the study config first.")
            return None
        pinned = [f"{k}: {pv!r} (toolbox: {tv!r})" for k, (pv, tv) in schema.pinned.items()]
        box = QMessageBox(self)
        box.setWindowTitle("Update project defaults")
        box.setText(
            f"<b>{schema.project_path}</b><br><br>"
            f"{len(schema.stale_keys)} key(s) added to PRESTUS are missing from this copy; "
            f"{len(pinned)} value(s) are pinned differently from the toolbox default.<br><br>"
            "Replacing the copy with the toolbox default keeps a dated backup. "
            "<i>Keep pinned values</i> first writes the pinned values into the study config, "
            "so the simulation settings stay the same. The study config is saved.")
        box.setDetailedText("Pinned values:\n" + ("\n".join(pinned) or "(none)") +
                            "\n\nMissing keys:\n" + ("\n".join(schema.stale_keys) or "(none)"))
        keep = box.addButton("Keep pinned values", QMessageBox.ButtonRole.AcceptRole)
        plain = box.addButton("Use toolbox values", QMessageBox.ButtonRole.DestructiveRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        clicked = box.clickedButton()
        if clicked not in (keep, plain):
            return None
        if clicked is keep:
            carry_pinned_values(self.model, schema.pinned)
        self.model.save(header=study_header(self.model.path))
        backup = replace_project_default(schema.project_path, toolbox_default(self.prestus_root))
        self.dirty = False
        model = ConfigModel.from_file(
            Schema(base_default_for(self.model.path, self.prestus_root),
                   toolbox_path=toolbox_default(self.prestus_root)),
            self.model.path, self.model.library)
        window = self._switch_to(model)
        window.statusBar().showMessage(f"Project defaults refreshed; backup: {backup.name}", 10000)
        return window

    def save(self) -> bool:
        if self.model.path is None:
            return self.save_as()
        self.model.save(header=study_header(self.model.path))
        self.dirty = False
        self.session.add_recent(self.model.path)
        self.session.save()
        self._update_title()
        self.statusBar().showMessage(f"Saved {self.model.path}", 5000)
        return True

    def save_as(self) -> bool:
        start = str(self.model.path or Path.home() / "config_study.yaml")
        chosen, _ = QFileDialog.getSaveFileName(self, "Save study config", start, "YAML (*.yaml)")
        if not chosen:
            return False
        target = Path(chosen)
        new_base = base_default_for(target, self.prestus_root)
        if new_base.resolve() != self.model.schema.default_path.resolve():
            QMessageBox.information(
                self, "Different defaults",
                f"In this folder, load_parameters will merge onto\n{new_base}\n"
                f"instead of\n{self.model.schema.default_path}.\n\n"
                "Reopen the file to edit it against those defaults.")
        self.model.path = target
        return self.save()

    def prune_defaults(self) -> None:
        n = self.model.prune_defaults()
        for page in self.pages.values():
            page.refresh()
        if n:
            self._changed()
        self.statusBar().showMessage(f"Dropped {n} value(s) equal to the defaults", 5000)

    def reload(self) -> None:
        if self.model.path is None or not self.model.path.is_file():
            return
        if self.dirty and QMessageBox.question(self, "Reload", "Discard unsaved changes?") \
                != QMessageBox.StandardButton.Yes:
            return
        self.model.load_dict(ConfigModel.from_file(self.model.schema, self.model.path, self.model.library).override_dict())
        self.dirty = False
        for page in self.pages.values():
            page.refresh()
        self._changed()
        self.dirty = False
        self._update_title()

    def closeEvent(self, event) -> None:
        run = self.pages.get(RUN)
        if run is not None and run.proc is not None:
            if QMessageBox.question(self, "MATLAB running", "Stop the running MATLAB job and quit?") \
                    != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            run.stop()
        if self.dirty:
            answer = QMessageBox.question(
                self, "Unsaved changes", "Save changes before closing?",
                QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard
                | QMessageBox.StandardButton.Cancel)
            if answer == QMessageBox.StandardButton.Cancel or \
                    (answer == QMessageBox.StandardButton.Save and not self.save()):
                event.ignore()
                return
        self.session.save()
        event.accept()


class LocaliteActions(QWidget):
    """Shortcuts shown with placement.mode = localite: ingest in MATLAB, build study runs."""

    def __init__(self, window: MainWindow):
        super().__init__()
        self.window = window
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 6)
        note = QLabel("Localite: ingest recordings on the Localite page, then use them for single runs "
                      "(session / target below) or turn them into study runs.")
        note.setObjectName("CardSub")
        note.setWordWrap(True)
        h.addWidget(note, 1)
        page = QPushButton("Localite page")
        page.setToolTip("Detect, name and write the recorded series of a session")
        page.clicked.connect(lambda: window.show_special(LOCALITE))
        h.addWidget(page)
        imp = QPushButton("Import into study…")
        imp.setToolTip("Build study manifest rows from the ingested sessions and targets")
        imp.clicked.connect(self._import)
        h.addWidget(imp)
        self.refresh()

    def _import(self) -> None:
        study = self.window.show_special(STUDY)
        study.import_neuronav()

    def refresh(self) -> None:
        self.setVisible(self.window.model.get("placement.mode") == "localite")
