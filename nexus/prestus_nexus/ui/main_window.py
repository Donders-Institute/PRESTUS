"""Main window: a header, one tab per area (setup, transducer, placement, simulation,
study, run), the issues panel and the review."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QByteArray, Qt, QTimer
from PyQt6.QtGui import QAction, QColor, QIcon, QKeySequence, QPainter, QPixmap
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFileDialog, QFrame, QHBoxLayout,
                             QLabel, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea,
                             QSplitter, QTabWidget, QTreeWidget, QTreeWidgetItem,
                             QVBoxLayout, QWidget)

from .. import checks
from ..fsio import copy_atomic
from ..telemetry import show_consent_dialog
from ..undo import History
from ..model import TX, ConfigModel, is_empty, iter_leaves
from ..paths import DEFAULT_NAME, base_default_for, toolbox_default
from ..project_defaults import carry_pinned_values, replace_project_default
from ..schema import GENERAL, Schema
from ..session import Session
from . import icons, style
from .errors import guarded
from .placement_viewer import QC_PATH
from .section import SectionPage
from .tabbar import tab_widget

REVIEW = "__review__"
VIEWER = "__viewer__"
RUN = "__run__"
STUDY = "__study__"
RESULTS = "__results__"
PCT = "__pct__"
CALIBRATION = "__calibration__"
LOCALITE = "__localite__"
PLANTUS = "__plantus__"
START = "__start__"
TUNE_FOR = False                     # 'Tune for' control: off until it suggests values (overlay 'tuning' kept)
MEDIUM_SELECT = "__medium__"         # overlay group of simulation.medium: the Medium tab header row
PLACEMENT_SELECT = "__placement_method__"   # overlay group of placement.mode: the Placement tab header row
TUS_ENTRY = "__tus_entry__"
HEURISTIC = "__heuristic__"
SPECIAL = (REVIEW, VIEWER, RUN, STUDY, RESULTS, PCT, CALIBRATION, LOCALITE, PLANTUS, TUS_ENTRY, HEURISTIC, START,
           MEDIUM_SELECT, PLACEMENT_SELECT)
SELECTORS = (MEDIUM_SELECT, PLACEMENT_SELECT)   # header rows that own their key (goto)
# Sub-tab label -> icon (ui/icons.py)
SUBTAB_ICONS = {
    "Settings": "advanced", "Calibration": "chart", "Viewer": "mri", "Beam axis": "wave", "3D": "cube",
    "Transducers": "tx_annular", "Localite": "place_localite", "Heuristic": "place_heuristic",
    "PlanTUS": "place_plantus", "TUS entry": "place_tus_entry", "Head model": "segmentation",
    "Acoustic properties": "table", "Uncertainty": "range", "Skull CT": "ct", "Stages": "stages", "Where to run": "server",
    "Plan": "report", "Review YAML": "file", "Launch": "play", "Jobs": "jobs", "Results": "chart",
}
SUBTAB_PAGES = (VIEWER, PCT, RESULTS, CALIBRATION, LOCALITE, PLANTUS, TUS_ENTRY, HEURISTIC)
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


MAX_FORM_WIDTH = 1240     # cards of a tab without a side panel: wider lines are hard to read


def _centred(body: QWidget) -> QWidget:
    """`body`, at most MAX_FORM_WIDTH wide and centred in the space it is given."""
    body.setMaximumWidth(MAX_FORM_WIDTH)
    holder = QWidget()
    h = QHBoxLayout(holder)
    h.setContentsMargins(0, 0, 0, 0)
    h.addStretch(1)
    h.addWidget(body, 100)
    h.addStretch(1)
    return holder


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
        self.headers: list[QWidget] = []               # rows between the intro and the cards
        self.columns: list[QVBoxLayout] = []           # card columns of further sub-tabs
        self.section_columns: dict[str, QVBoxLayout] = {}   # section -> column other than 'Settings'
        self.on_refresh = None                          # tab-specific refresh (e.g. sub-tab visibility)

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
        self.advanced = QPushButton(icons.icon("advanced"), "Advanced")
        self.advanced.setObjectName("AdvancedToggle")
        self.advanced.setCheckable(True)
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
        self._body = body
        self.scroll.setWidget(_centred(body))
        self.split.addWidget(self.scroll)
        outer.addWidget(self.split, 1)

    def add_card(self, widget: QWidget) -> None:
        self.cards.append(widget)
        self.column.addWidget(widget)

    def add_section(self, page: SectionPage, column: QVBoxLayout | None = None) -> None:
        page.set_show_advanced(self.advanced.isChecked())
        page.advanced_toggled.connect(self.refresh)
        self.sections.append(page)
        if column is None:                  # (an empty layout is falsy: no 'or' chain)
            column = self.section_columns.get(page.section, self.column)
        column.addWidget(page)

    def add_header(self, widget: QWidget) -> None:
        """A row under the intro, above the cards and sub-tabs (refreshed with the tab)."""
        self.headers.append(widget)
        self.layout().insertWidget(len(self.headers), widget)

    @staticmethod
    def new_column() -> tuple[QScrollArea, QVBoxLayout]:
        """A scrolling column of cards, as the 'Settings' one, for another sub-tab."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        column = QVBoxLayout(body)
        column.setContentsMargins(0, 4, 6, 12)
        column.setSpacing(12)
        scroll.setWidget(_centred(body))
        return scroll, column

    def subtab_of(self, widget: QWidget) -> QWidget | None:
        """The sub-tab page that holds `widget`."""
        if self.subtabs is None:
            return None
        for i in range(self.subtabs.count()):
            page = self.subtabs.widget(i)
            if page is widget or page.isAncestorOf(widget):
                return page
        return None

    def set_side(self, widget: QWidget, stretch: tuple[int, int] = (3, 2)) -> None:
        self.side = widget
        self.scroll.setMinimumWidth(560)
        self._body.setMaximumWidth(16_777_215)      # beside a panel the cards take their share
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
        self.subtabs = tab_widget()
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
            self.subtabs = tab_widget()
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
        for column in self.columns:
            column.addStretch(1)
        tunable = TUNE_FOR and any(r.node.tune for page in self.sections for r in page.rows)
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
            page.setVisible(visible and any(page._row_shown(r) for r in page.rows)   # also while collapsed
                            or (visible and (page.section == TX or page.table is not None)))
            n_hidden += page.hidden_advanced() if visible else 0
        for card in [*self.cards, *self.headers]:
            if hasattr(card, "refresh"):
                card.refresh()
        if self.on_refresh is not None:
            self.on_refresh(self)
        # the count of hidden settings sits in the toggle itself
        self.hidden.setText("")
        self.advanced.setText(f"Advanced ({n_hidden})" if n_hidden and not show else "Advanced")
        self.advanced.setToolTip(
            (f"{n_hidden} advanced setting{'s' if n_hidden != 1 else ''} hidden on this tab. " if n_hidden and not show
             else "") + "Show the settings most studies leave at their defaults. "
            "Advanced settings that hold a value are always shown.")


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
        self.history = History(model, parent=self)      # Edit ▸ Undo / Redo
        self.history.restored.connect(self._history_restored)
        self._restore_geometry()

        central = QWidget()
        cv = QVBoxLayout(central)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(0)
        cv.addWidget(self._header())
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.tabBar().setExpanding(True)       # the workflow steps span the window
        self.tabs.tabBar().setUsesScrollButtons(False)
        cv.addWidget(self.tabs, 1)
        self.setCentralWidget(central)

        self.issue_popup = self._issue_popup()

        self._build_menu()
        self._build_edit_menu()
        self._build_view_menu()
        self._build_tools_menu()
        self.update_context()
        self._build_tabs()
        self._refresh_cal_cov()
        self._apply_subtab_icons()
        self._build_help_menu()
        self._build_shortcuts()
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
        open_btn = QPushButton(icons.icon("folder", style.C["ink2"]), "Open…")
        open_btn.setToolTip("Open a study config (File ▸ Open…)")
        open_btn.clicked.connect(lambda: self.open_study())
        h.addWidget(open_btn)
        save = QPushButton(icons.icon("save", "#ffffff"), "Save")
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
        for text, slot in [("Localite", lambda: self.show_special(LOCALITE)),
                           ("Calibration", lambda: self.show_special(CALIBRATION)),
                           (None, None),
                           ("neuronav_gui in MATLAB…", self.open_neuronav_gui),
                           ("calibration_gui in MATLAB…", self.open_calibration_gui),
                           (None, None),
                           ("Update project defaults…", self.update_project_defaults),
                           (None, None),
                           ("Usage statistics…", lambda: show_consent_dialog(self, self.prestus_root))]:
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

    def _build_edit_menu(self) -> None:
        edit = self.menuBar().addMenu("&Edit")
        self.undo_action = QAction("Undo", self)
        self.undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        self.undo_action.triggered.connect(self.history.undo)
        self.redo_action = QAction("Redo", self)
        self.redo_action.setShortcuts([QKeySequence.StandardKey.Redo, QKeySequence("Ctrl+Shift+Z")])
        self.redo_action.triggered.connect(self.history.redo)
        edit.addAction(self.undo_action)
        edit.addAction(self.redo_action)
        self.history.changed.connect(self._update_undo_actions)
        self._update_undo_actions()

    def _update_undo_actions(self) -> None:
        self.undo_action.setEnabled(self.history.can_undo())
        self.undo_action.setText(self.history.undo_text())
        self.redo_action.setEnabled(self.history.can_redo())
        self.redo_action.setText(self.history.redo_text())

    def _history_restored(self, _label: str) -> None:
        for page in self.pages.values():
            page.refresh()
        self._changed()

    def _current_tab(self):
        """(tab key, sub-tab label) shown now, for Help ▸ Docs for this tab."""
        w = self.tabs.currentWidget()
        for key, tab in self.tab_pages.items():
            if tab is w:
                subs = getattr(tab, "subtabs", None)
                return key, (subs.tabText(subs.currentIndex()) if subs is not None else None)
        return "quickstart"

    def _build_help_menu(self) -> None:
        from .. import docs
        help_menu = self.menuBar().addMenu("&Help")
        tab_docs = docs.build_help_menu(help_menu, self._current_tab, parent=self, session=self.session)
        tab_docs.setShortcut(QKeySequence("F1"))

    def _build_shortcuts(self) -> None:
        """Ctrl+1…9 jump to the workflow tabs; Ctrl+Shift+V validates, Ctrl+R runs locally."""
        for i in range(min(9, self.tabs.count())):
            action = QAction(self)
            action.setShortcut(QKeySequence(f"Ctrl+{i + 1}"))
            action.triggered.connect(lambda _c=False, i=i: self.tabs.setCurrentIndex(i))
            self.addAction(action)
        run = self.pages.get(RUN)
        if run is not None:
            for keys, slot in (("Ctrl+Shift+V", run.validate), ("Ctrl+R", run.run)):
                action = QAction(self)
                action.setShortcut(QKeySequence(keys))
                action.triggered.connect(lambda _c=False, slot=slot: (self.show_special(RUN), slot()))
                self.addAction(action)

    def _set_dark_mode(self, dark: bool) -> None:
        style.apply(QApplication.instance(), self.prestus_root, dark=dark)
        self.session.set("dark_mode", dark)
        self.session.save()
        self._apply_subtab_icons()
        for page in self.pages.values():
            if isinstance(page, SectionPage):
                page.retint()

    def update_context(self) -> None:
        """Facts about this computer that visibility rules read (ConfigModel.context):
        '__scheduler' when sbatch/qsub is here or runs are submitted over SSH."""
        from .. import jobs as hpc_jobs
        self.model.context["__scheduler"] = bool(hpc_jobs.scheduler()) or \
            bool(self.session.get("remote_enabled", False))

    def _apply_subtab_icons(self) -> None:
        """An icon per sub-tab (SUBTAB_ICONS), in the active palette."""
        for i in range(self.tabs.count()):
            tab = self.tabs.widget(i)
            subtabs = getattr(tab, "subtabs", None)
            if subtabs is None:
                continue
            for j in range(subtabs.count()):
                name = SUBTAB_ICONS.get(subtabs.tabText(j))
                if name:
                    subtabs.setTabIcon(j, icons.icon(name))

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
        self.tabs.addTab(start, "QuickStart")
        self.tabs.setTabToolTip(0, "A first simulation in four steps; everything else keeps its default")
        for spec in self._tab_specs():
            tab = TabPage(self, spec)
            key = spec["key"]
            if key == "medium":
                self._build_medium_tab(tab)
            if key == "run":
                # sub-tabs in workflow order: Stages, Where to run, Plan, Review YAML, Launch, Jobs
                where, where_col = tab.new_column()
                tab.add_subtab(where, "Where to run")
                tab.subtabs.setTabText(0, "Stages")
                tab.columns.append(where_col)
                for sec in ("hpc", "snellius"):
                    tab.section_columns[sec] = where_col
            if key == "study":
                from .study_page import StudyPage
                study = StudyPage(self)
                study.changed.connect(self._changed)
                self.pages[STUDY] = study
                self._section_tab[STUDY] = key
            for section in spec.get("sections", []):
                if section != GENERAL and schema.get(section) is None \
                        and not schema.section_is_virtual(section):
                    continue
                nodes = schema.display_nodes(section)
                if not nodes:
                    continue
                page = SectionPage(self.model, section, nodes, schema.section_label(section),
                                   subtitle=schema.section_subtitle(section),
                                   locked=schema.section_is_locked(section),
                                   icon=schema.sections.get(section, {}).get("icon"))
                page.lock_note = str(schema.sections.get(section, {}).get("lock_note") or "")
                page.changed.connect(self._changed)
                self.pages[section] = page
                self._section_tab[section] = key
                if section == "babelbrain":
                    # BabelBrain's HU mapping and internal air sit with the other CT settings
                    bb_ct_link = QLabel(
                        "CT-based skull properties (BabelBrain HU mapping, internal air) are set under "
                        "<a href='#pct'>Medium › CT-based skull mapping</a>, once a CT or pseudo-CT is used.")
                    bb_ct_link.setObjectName("CardSub")
                    bb_ct_link.setWordWrap(True)
                    bb_ct_link.linkActivated.connect(lambda _l: self.goto("pct.enabled"))
                    page.add_extra(bb_ct_link)
                if section == "segmentation":
                    from .segmentation_panel import SegmentationPanel
                    self.segmentation_panel = SegmentationPanel(self)
                    page.add_extra(self.segmentation_panel)
                if section == "placement":
                    page.add_extra(LocaliteActions(self))
                    page.add_extra(PlantusActions(self))
                if section == "segmentation":
                    from .section import LABEL_WIDTH
                    from .simnibs_panel import SimnibsPanel
                    from .. import remote as rmt
                    from .. import simnibs as sn

                    def _ssh_host(session=self.session):
                        h = rmt.host_from_session(session)
                        return h.host if h else None
                    sn.remote_host = _ssh_host          # checks: an install may be on the cluster
                    self.simnibs_panel = SimnibsPanel(self.model, LABEL_WIDTH, session=self.session)
                    self.simnibs_panel.page = page
                    self.simnibs_panel.changed.connect(self._changed)
                    page.add_extra(self.simnibs_panel)
                if section == "medium_properties":
                    from .medium_table import MediumPanel
                    self.medium_panel = MediumPanel(self.model, self.prestus_root)
                    self.medium_panel.goto.connect(self.goto)
                    page.set_table(self.medium_panel)
                if section == "coupling":
                    from .coupling_view import CouplingPrediction
                    self.coupling_view = CouplingPrediction(self.model, lambda: self.pages.get(VIEWER))
                    self.coupling_view.open_placement_requested.connect(lambda: self.show_special(VIEWER))
                    page.add_extra(self.coupling_view, below=True)
                if section == "timing":
                    from .timing_view import TimingPanel
                    self.timing_panel = TimingPanel(self.model)
                    self.timing_panel.changed.connect(self._changed)
                    page.add_extra(self.timing_panel)
                    page.set_fold(TimingPanel.FOLDED, open_=False)      # edited in the diagram
                    self.timing_panel.fold_toggled.connect(
                        lambda shown, pg=page: (setattr(pg, "fold_open", shown), pg.refresh()))
                if section == "hpc":
                    from .run_profile import RunProfilePanel
                    self.run_profile_panel = RunProfilePanel(self)
                    self.run_profile_panel.changed.connect(self._changed)
                    self.run_profile_panel.changed.connect(page.refresh)
                    page.add_extra(self.run_profile_panel)
                    from .pipeline_options import PipelineOptionsPanel
                    self.pipeline_options_panel = PipelineOptionsPanel(self)
                    page.add_extra(self.pipeline_options_panel)
                tab.add_section(page)
                if section == "study":
                    tab.add_card(self.pages[STUDY])
            if key == "setup" and "io" in self.pages:
                from .io_preview import IOPreviewPanel
                self.io_preview = IOPreviewPanel(self.model)
                self.pages["io"].extras_below.addWidget(self.io_preview)   # under the I/O fields
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
            elif key == "medium":
                # cards that follow the sections of their sub-tab
                tab.section_columns["uncertainty"].addWidget(self.uncertainty_card)
                if "uncertainty" in self.pages:          # variant files: in the card's advanced fold
                    self.uncertainty_card.attach_files(self.pages["uncertainty"])
                tab.section_columns["pct"].addWidget(self.pct_flow)
                tab.section_columns["pct"].addWidget(self.pct_formulas)
            elif key == "run":
                from .run_page import RunPage
                run = RunPage(self)
                review = ReviewPage(self)
                self.pages[RUN], self.pages[REVIEW] = run, review
                self._section_tab[RUN] = self._section_tab[REVIEW] = key
                tab.add_subtab(run.plan_card, "Plan")
                tab.add_subtab(review, "Review YAML")
                tab.add_subtab(run, "Launch")
                tab.add_subtab(run.jobs_view, "Jobs")
                self.run_side = tab.subtabs         # RunPage.show_jobs_tab and show_special select in it
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
            if tab.key == "medium" and PCT in self.pages:
                self.pages[PCT].refresh()
            if tab.key == "run":
                if VIEWER in self.pages and self.model.get("placement.mode") in ("manual", None, ""):
                    QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                    try:                            # check the placement before a run
                        self.pages[VIEWER].ensure_images()
                    finally:
                        QApplication.restoreOverrideCursor()
                    self._refresh_issues()
                current = self.run_side.currentWidget()
                if hasattr(current, "refresh"):
                    current.refresh()
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

    def _build_medium_tab(self, tab: TabPage) -> None:
        """Medium tab: the medium selector as a header row, then sub-tabs 'Head model'
        (segmentation, head model, tissue layers), 'Acoustic properties' (the property table),
        'Uncertainty' (the analysis switch and ITRUSST bounds plot; variant table and files in an advanced fold) and
        'Skull CT' (CT-based skull mapping settings, equations and the CT page). Sub-tabs the medium does not use are hidden (_refresh_medium_tab)."""
        from .medium_selector import MediumSelector
        from .pct_flow import PctFlow
        from .pct_formulas import PctFormulas
        from .pct_page import PctPage
        from .uncertainty_view import UncertaintyCard
        self.medium_selector = MediumSelector(self.model)
        self.medium_selector.changed.connect(self._changed)
        self.pages[MEDIUM_SELECT] = self.medium_selector
        self._section_tab[MEDIUM_SELECT] = tab.key
        tab.add_header(self.medium_selector)

        scroll, column = tab.new_column()
        tab.add_subtab(scroll, "Acoustic properties")
        tab.subtabs.setTabText(0, "Head model")
        tab.columns.append(column)
        tab.section_columns["medium_properties"] = column

        # its own sub-tab: switching the analysis on does not reflow the property table
        unc_scroll, unc_column = tab.new_column()
        tab.add_subtab(unc_scroll, "Uncertainty")
        tab.columns.append(unc_column)
        tab.section_columns["uncertainty"] = unc_column
        self.uncertainty_card = UncertaintyCard(self.model, self.prestus_root)
        tab.cards.append(self.uncertainty_card)

        pct_scroll, pct_column = tab.new_column()
        tab.columns.append(pct_column)
        tab.section_columns["pct"] = pct_column
        self.pct_flow = PctFlow(self)
        tab.cards.append(self.pct_flow)
        self.pct_formulas = PctFormulas(self.model, self.prestus_root)
        tab.cards.append(self.pct_formulas)
        pct_page = PctPage(self)
        self.pages[PCT] = pct_page
        self._section_tab[PCT] = tab.key
        split = QSplitter()
        split.addWidget(pct_scroll)
        split.addWidget(pct_page)
        pct_scroll.setMinimumWidth(520)           # the settings and equations keep their width
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setChildrenCollapsible(False)
        split.setSizes([1000, 1000])
        split.refresh = pct_page.refresh          # TabPage._sub_tab_changed
        tab.add_subtab(split, "Skull CT")
        self.medium_subtabs = {"head": tab.settings_page, "acoustic": scroll, "uncertainty": unc_scroll,
                               "pct": split}
        tab.on_refresh = self._refresh_medium_tab

    def _refresh_medium_tab(self, tab: TabPage) -> None:
        """Water uses no head model and no skull CT; a phantom uses tissue layers but no skull CT."""
        medium = str(self.model.get("simulation.medium") or "layered")
        # the variant files bound brain and skull: nothing to vary in water
        shown = {"head": medium != "water", "acoustic": True, "uncertainty": medium != "water" or self.model.get("simulation.uncertainty") is True,
                 "pct": medium == "layered"}
        for name, widget in self.medium_subtabs.items():
            tab.subtabs.setTabVisible(tab.subtabs.indexOf(widget), shown[name])
        layers = self.pages.get("layers")
        if layers is not None and medium == "phantom":
            # a phantom is defined by its layers: shown (still locked) though the card is advanced
            layers.setVisible(True)
        current = next((n for n, w in self.medium_subtabs.items() if w is tab.subtabs.currentWidget()), None)
        if current is not None and not shown[current]:
            tab.subtabs.setCurrentWidget(self.medium_subtabs["acoustic"])

    def goto(self, path: str, ti: int | None = None) -> None:
        section = self.display_section(path)
        page = self.pages.get(section)
        if section in SELECTORS and page is not None:
            self._select_tab(self._section_tab.get(section, ""))
            page.focus_key(path, ti)
            return
        if not isinstance(page, SectionPage):
            return
        tab = self._select_tab(self._section_tab.get(section, ""))
        if tab is not None:
            if tab.subtabs is not None:
                holder = tab.subtab_of(page) or tab.settings_page
                tab.subtabs.setCurrentWidget(holder)      # the sub-tab the field lives on
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
        if run is not None and run.plan_card.isVisible():   # the Plan sub-tab refreshes it when shown
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
        with guarded(self, f"Opening {path.name}") as g:
            model = ConfigModel.from_file(Schema(base_default_for(path, self.prestus_root),
                                                 toolbox_path=toolbox_default(self.prestus_root)),
                                          path, self.model.library)
        if g.failed:
            return None
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
        with guarded(self, "Creating the study config") as g:
            if copy_default:
                copy_atomic(toolbox_default(self.prestus_root), local)
            model = ConfigModel(Schema(base_default_for(target, self.prestus_root),
                                       toolbox_path=toolbox_default(self.prestus_root)), self.model.library)
            model.save(target, header=study_header(target))
        if g.failed:
            return None
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
        with guarded(self, "Updating the project defaults") as g:
            if clicked is keep:
                carry_pinned_values(self.model, schema.pinned)
            self.model.save(header=study_header(self.model.path))
            self.dirty = False
            backup = replace_project_default(schema.project_path, toolbox_default(self.prestus_root))
            model = ConfigModel.from_file(
                Schema(base_default_for(self.model.path, self.prestus_root),
                       toolbox_path=toolbox_default(self.prestus_root)),
                self.model.path, self.model.library)
        if g.failed:
            return None
        window = self._switch_to(model)
        window.statusBar().showMessage(f"Project defaults refreshed; backup: {backup.name}", 10000)
        return window

    def save(self) -> bool:
        if self.model.path is None:
            return self.save_as()
        with guarded(self, "Saving the study config") as g:
            self.model.save(header=study_header(self.model.path))
        if g.failed:
            return False
        self.history.seal()
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
        previous, self.model.path = self.model.path, target
        if not self.save():
            self.model.path = previous
            return False
        return True

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
        with guarded(self, f"Reloading {self.model.path.name}") as g:
            fresh = ConfigModel.from_file(self.model.schema, self.model.path, self.model.library)
        if g.failed:
            return
        self.model.load_dict(fresh.override_dict())
        self.history.clear()
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
        self.session.set("window_geometry", bytes(self.saveGeometry().toBase64()).decode())
        self.session.save()
        event.accept()

    def _restore_geometry(self) -> None:
        """Where the window was last closed; the first time, most of the screen (capped, so a
        large monitor does not get a window of mostly empty margins), centred."""
        self.setMinimumSize(900, 560)
        saved = self.session.get("window_geometry")
        if isinstance(saved, str) and saved and \
                self.restoreGeometry(QByteArray.fromBase64(saved.encode())):
            return
        screen = QApplication.primaryScreen()
        if screen is None:
            self.resize(1320, 900)
            return
        area = screen.availableGeometry()
        w = min(1500, int(area.width() * 0.9))
        h = min(980, int(area.height() * 0.9))
        self.resize(max(w, 900), max(h, 560))
        self.move(area.x() + (area.width() - self.width()) // 2,
                  area.y() + (area.height() - self.height()) // 2)


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


class PlantusActions(QWidget):
    """Shown with placement.mode = plantus: run only the placement, then show it in the viewer."""

    TASK = "PlanTUS placement"

    def __init__(self, window: MainWindow):
        super().__init__()
        self.window = window
        self._waiting = False
        h = QHBoxLayout(self)
        h.setContentsMargins(0, 0, 0, 6)
        note = QLabel("PlanTUS: run only the placement, without simulating. <b>Choose in Workbench</b> "
                      "opens PlanTUS's maps (needs Connectome Workbench and a display); click a skin "
                      "position, close Workbench, then use it or choose again in the review dialog. "
                      "<b>Automatic</b> takes the best composite score. The result opens in the Viewer.")
        note.setObjectName("CardSub")
        note.setWordWrap(True)
        h.addWidget(note, 1)
        pick = QPushButton("Choose in Workbench…")
        pick.setToolTip("Run PlanTUS and pick the transducer position on its criterion maps")
        pick.clicked.connect(lambda: self._run(True))
        h.addWidget(pick)
        auto = QPushButton("Automatic")
        auto.setToolTip("Run PlanTUS and place the transducer at the composite-score optimum")
        auto.clicked.connect(lambda: self._run(False))
        h.addWidget(auto)
        self.refresh()

    def _run(self, interactive: bool) -> None:
        from .. import launcher
        if not self.window.save():
            return
        script = launcher.write_plantus_script(self.window.model.path, self.window.prestus_root,
                                               interactive)
        page = self.window.show_special(RUN)
        if page.run_script(self.TASK, script) and not self._waiting:
            self._waiting = True
            page.finished.connect(self._done)

    def _done(self, task: str, ok: bool) -> None:
        if task != self.TASK:
            return
        self.window.pages[RUN].finished.disconnect(self._done)
        self._waiting = False
        if ok and VIEWER in self.window.pages:
            self.window.show_special(VIEWER)

    def refresh(self) -> None:
        self.setVisible(self.window.model.get("placement.mode") == "plantus")
