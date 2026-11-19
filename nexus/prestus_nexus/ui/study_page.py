"""Study page: edit the study manifest (one row per run) and start study runs / the group report.

The UI calls the manifest the 'table' and its rows 'runs'; config keys and MATLAB names
stay in tooltips."""

from __future__ import annotations

import contextlib
import copy
import re
from pathlib import Path

from PyQt6.QtCore import Qt, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (QAbstractItemView, QFileDialog, QFrame, QHBoxLayout, QInputDialog, QLabel,
                             QLineEdit, QMenu, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem,
                             QToolButton, QVBoxLayout, QWidget)

from ..labels import plural
from .. import launcher
from ..study import (Manifest, from_neuronav_summary, parse_subject_ids, resolve_manifest, study_report_path,
                     subject_list_report_path)
from . import icons, style
from .empty_state import EmptyStack

LEAD = ("Simulate multiple subjects, or several targets per subject, in one go. Each row of the table "
        "is one simulation run; everything not in the table comes from this study config. "
        "For a single subject you do not need this: use the <b>QuickStart</b> tab.")

HOW_TO = (
    "<ol style='margin-left:-20px'>"
    "<li><b>Add the runs.</b> <i>Add subjects…</i> takes IDs such as <tt>1-12, 15</tt>. "
    "<i>Import from Localite…</i> makes one run per recorded target.</li>"
    "<li><b>Name each condition</b> if a subject has more than one run (e.g. the target "
    "<tt>L_M1</tt> and <tt>R_M1</tt>). The condition names the outputs and groups the group report. "
    "Two included rows with the same subject and condition would overwrite each other.</li>"
    "<li><b>Optional: vary a setting per run.</b> <i>Add setting column…</i> adds any config key, "
    "e.g. <tt>transducer.target_isppa_wcm2</tt>. An empty cell keeps the study value.</li>"
    "<li><b>Check, then run.</b> The table is saved next to the study config when you check or run.</li>"
    "</ol>"
    "Example:"
    "<table cellspacing='0' cellpadding='3' border='1' style='margin-top:4px'>"
    "<tr><th>subject_id</th><th>condition</th><th>transducer.target_isppa_wcm2</th></tr>"
    "<tr><td>1</td><td>L_M1</td><td></td></tr>"
    "<tr><td>1</td><td>R_M1</td><td></td></tr>"
    "<tr><td>2</td><td>L_M1</td><td>20</td></tr></table>")

# Hover help for the reserved columns; config-key columns show their config_default.yaml help
RESERVED_HELP = {
    "subject_id": "Required. The subject: 3, 003, or sub-003.",
    "condition": "Optional run label, e.g. the target or protocol. It is added to the output names "
                 "(_desc-<condition>) and groups the group report. Leave empty with one run per subject.",
    "include": "Optional. 1/0, yes/no, or true/false; empty = include. Set 0 to skip a run without deleting it.",
    "config": "Optional YAML file merged over the study config for this run only (path relative to "
              "the manifest's folder). Use it when a run differs in many settings.",
    "notes": "Free text; ignored by PRESTUS.",
}

FIRST_ADD = "Each run takes the study config's settings; add a setting column to vary one per run."


class StudyPage(QFrame):
    changed = pyqtSignal()

    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.manifest = Manifest()
        self._loading = False

        self.setObjectName("Card")
        v = QVBoxLayout(self)
        v.setContentsMargins(16, 12, 16, 14)
        title = QLabel("Study runs")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        lead = QLabel(LEAD)
        lead.setObjectName("CardSub")
        lead.setWordWrap(True)
        v.addWidget(lead)
        self.b_how = QToolButton(text="How to fill the table")
        self.b_how.setCheckable(True)
        self.b_how.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.b_how.setArrowType(Qt.ArrowType.RightArrow)
        self.b_how.setAutoRaise(True)
        self.how = QLabel(HOW_TO)
        self.how.setObjectName("CardSub")
        self.how.setWordWrap(True)
        self.how.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.how.hide()
        self.b_how.toggled.connect(self._toggle_how)
        v.addWidget(self.b_how)
        v.addWidget(self.how)

        # Adding runs first; editing the table's shape second
        edit = QHBoxLayout()
        self.b_add = QPushButton("Add subjects…")
        self.b_add.setProperty("primary", True)
        self.b_add.setToolTip("One run per subject ID, e.g. 1-12, 15")
        self.b_add.clicked.connect(lambda: self.add_subjects())
        self.b_import = QPushButton("Import from Localite…")
        self.b_import.setToolTip("One run per target recorded in Localite (after ingestion on the Localite page)")
        self.b_import.clicked.connect(lambda: self.import_neuronav())
        b_col = QPushButton("Add setting column…")
        b_col.setToolTip("Vary a config setting per run; empty cells keep the study value")
        b_col.clicked.connect(lambda: self.add_column())
        self.b_more = QToolButton(text="More ▾")
        self.b_more.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.b_more.setStyleSheet("QToolButton::menu-indicator { image: none; width: 0; }")
        more = QMenu(self.b_more)
        more.addAction("Add empty row", lambda: self.add_row())
        more.addAction("Remove selected rows", self.remove_rows)
        more.addAction("Remove selected column", self.remove_column)
        more.addSeparator()
        more.addAction("Save table", lambda: self.save_manifest())
        more.addAction("Save table as…", self.save_as)
        self.b_more.setMenu(more)
        for b in (self.b_add, self.b_import, b_col, self.b_more):
            edit.addWidget(b)
        edit.addStretch(1)
        v.addLayout(edit)

        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.itemChanged.connect(self._cell_edited)
        self.table.setMinimumHeight(220)
        self.stack = EmptyStack(self.table)
        v.addWidget(self.stack, 1)

        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.summary)

        # Where the table lives: a detail, so after the table and pre-filled on first save
        row = QHBoxLayout()
        lbl = QLabel("Table file")
        lbl.setToolTip("Stored as study.manifest in the study config; relative to the config's folder")
        row.addWidget(lbl)
        self.path = QLineEdit()
        self.path.setPlaceholderText("set on first save (default: <config name>_runs.tsv next to the config)")
        self.path.editingFinished.connect(self._path_edited)
        browse = icons.PathButton("file", "Choose or create a .tsv/.csv table")
        browse.clicked.connect(self._browse)
        row.addWidget(self.path, 1)
        row.addWidget(browse)
        v.addLayout(row)

        run = QHBoxLayout()
        self.b_dry = QPushButton("Check all runs")
        self.b_dry.setToolTip("Saves the table and loads every run's config in MATLAB without simulating "
                              "(prestus_study_start(config, 'DryRun', true))")
        self.b_dry.clicked.connect(self.dry_run)
        self.b_run = QPushButton("Run all runs locally")
        self.b_run.setToolTip("Writes run_study_<study>.m and runs it here; output streams to the Run tab")
        self.b_run.clicked.connect(self.run_study)
        self.b_launcher = QPushButton("Write launcher script")
        self.b_launcher.setToolTip("Writes run_study_<study>.m to run yourself, e.g. on the cluster")
        self.b_launcher.clicked.connect(self.write_launcher)
        self.b_report = QPushButton("Build group report")
        self.b_report.setToolTip("After the runs finished: group MNI maps + group report from their outputs")
        self.b_report.clicked.connect(self.group_report)
        self.b_open_report = QPushButton("Open group report")
        self.b_open_report.setToolTip("Open the study's group report (path.sim/group/<study>/"
                                      "group_<study>_report.html) in the web browser")
        self.b_open_report.clicked.connect(self.open_group_report)
        for b in (self.b_dry, self.b_run, self.b_launcher, self.b_report, self.b_open_report):
            run.addWidget(b)
        run.addStretch(1)
        v.addLayout(run)

    def _toggle_how(self, on: bool) -> None:
        self.how.setVisible(on)
        self.b_how.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)

    # ── state ───────────────────────────────────────────────────────────
    @property
    def model(self):
        return self.window.model

    # edit history (MainWindow.history tracks this page: the table undoes with the config)
    def snapshot(self) -> tuple:
        m = self.manifest
        return (m.path, list(m.columns), copy.deepcopy(m.rows))

    def restore(self, snap: tuple) -> None:
        path, columns, rows = snap
        self.manifest = Manifest(columns=list(columns), rows=copy.deepcopy(rows), path=path)

    def _edit(self, label: str):
        """One undo step for the table changes made inside."""
        history = getattr(self.window, "history", None)
        return history.group(label) if history is not None else contextlib.nullcontext()

    def manifest_path(self) -> Path | None:
        return resolve_manifest(self.model.path, self.model.get("study.manifest") or None)

    def refresh(self) -> None:
        value = self.model.get("study.manifest") or ""
        if self.path.text() != value:
            self.path.setText(value)
        mp = self.manifest_path()
        if mp is not None and mp != self.manifest.path and mp.is_file():
            self.manifest = Manifest.load(mp)
        self._fill_table()

    def _fill_table(self) -> None:
        self._loading = True
        m = self.manifest
        self.table.setColumnCount(len(m.columns))
        self.table.setHorizontalHeaderLabels(m.columns)
        for j, c in enumerate(m.columns):
            self.table.horizontalHeaderItem(j).setToolTip(self._column_help(c))
        self.table.setRowCount(len(m.rows))
        for i, r in enumerate(m.rows):
            for j, c in enumerate(m.columns):
                self.table.setItem(i, j, QTableWidgetItem(r.get(c, "")))
        self.table.resizeColumnsToContents()
        self._loading = False
        for b in (self.b_dry, self.b_run, self.b_launcher):
            b.setEnabled(bool(m.rows))
        if m.rows:
            self.stack.show_content()
        else:
            self.stack.show_empty(
                "No runs yet",
                "Add the subjects to simulate (one run each), or import the targets recorded in "
                "Localite with <i>Import from Localite…</i> above.",
                "Add subjects…", self.add_subjects)
        self._update_summary()

    def _column_help(self, column: str) -> str:
        if column.lower() in RESERVED_HELP:
            return f"<b>{column}</b><br>{RESERVED_HELP[column.lower()]}"
        node = self.model.schema.get(re.sub(r"\(\d+\)", "", column))
        if node is None:
            return f"<b>{column}</b><br>Not a key of config_default.yaml."
        from .fields import display_text
        from .tooltip import html_for
        return html_for(node, column, display_text(node.default), note="Empty cell = the study config's value.")

    def _update_summary(self) -> None:
        m = self.manifest
        affix = self.model.get("io.output_affix") or ""
        issues = m.issues(self.model.schema, affix if isinstance(affix, str) else "")
        n_inc = sum(1 for r in m.rows if r.get("include", "").strip().lower() not in {"0", "false", "no", "n"})
        conds = [c or "(none)" for c in m.conditions()]
        if not m.rows:
            self.summary.setText("")
            return
        # counts on the first line, the conditions on their own, then the check result
        text = (f"{plural(len(m.rows), 'row')}, {n_inc} included, {plural(len(m.subjects()), 'subject')}"
                f"<br>Conditions: {', '.join(conds) or '—'}")
        if not issues:
            text += f"<br><span style='color:{style.active()['green']}'>No problems found</span>"
        for iss in issues:
            colour = style.active()["red"] if iss.level == "error" else style.active()["amber"]
            text += f"<br><span style='color:{colour}'>{iss.level}: {iss.where()}: {iss.message}</span>"
        self.summary.setText(text)

    def _cell_edited(self, item: QTableWidgetItem) -> None:
        if self._loading:
            return
        col = self.manifest.columns[item.column()]
        with self._edit("Edit run"):
            self.manifest.rows[item.row()][col] = item.text().strip()
        self._update_summary()

    def _path_edited(self) -> None:
        text = self.path.text().strip()
        if text == (self.model.get("study.manifest") or ""):
            return
        with self._edit("Change runs table file"):        # the path and the table it loads
            if text:
                self.model.set("study.manifest", text)
            else:
                self.model.reset("study.manifest")
            mp = self.manifest_path()
            if mp is not None and mp.is_file():
                self.manifest = Manifest.load(mp)
            else:
                # A new file: the runs entered so far are saved there
                self.manifest.path = mp
        self.changed.emit()
        self._fill_table()

    def _set_manifest_path(self, path: Path) -> None:
        cfg = self.model.path
        try:
            value = str(path.relative_to(cfg.parent)) if cfg is not None else str(path)
        except ValueError:
            value = str(path)
        self.path.setText(value)
        self._path_edited()

    def _browse(self) -> None:
        start = str(self.manifest_path() or self.default_path() or Path.home())
        chosen, _ = QFileDialog.getSaveFileName(self, "Study runs table", start,
                                                "Manifest (*.tsv *.csv)",
                                                options=QFileDialog.Option.DontConfirmOverwrite)
        if chosen:
            self._set_manifest_path(Path(chosen))

    # ── editing ─────────────────────────────────────────────────────────
    def add_row(self, values: dict[str, str] | None = None) -> None:
        with self._edit("Add run"):
            self.manifest.add_row(values)
        self._fill_table()

    def add_subjects(self, text: str | None = None) -> None:
        """One run per subject ID ('1-12, 15'); the condition column is left empty to fill in."""
        if text is None:
            text, ok = QInputDialog.getText(self, "Add subjects",
                                            "Subject IDs (ranges, commas, or spaces), e.g. 1-12, 15, sub-020:")
            if not ok:
                return
        try:
            ids = parse_subject_ids(text)
        except ValueError as e:
            QMessageBox.warning(self, "Add subjects", str(e))
            return
        first = not self.manifest.rows
        with self._edit("Add subjects"):
            for sid in ids:
                self.manifest.add_row({"subject_id": sid})
        self._fill_table()
        if ids:
            self.window.statusBar().showMessage(f"Added {plural(len(ids), 'run')}. " + (FIRST_ADD if first else ""), 8000)

    def remove_rows(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        with self._edit("Remove runs"):
            for r in rows:
                del self.manifest.rows[r]
        self._fill_table()

    def add_column(self, key: str | None = None) -> None:
        if key is None:
            keys = sorted(n.path for n in self.model.schema.leaves()
                          if not n.path.startswith(("study.", "path.")) and n.path != "subject_id")
            key, ok = QInputDialog.getItem(
                self, "Add setting column",
                "Setting to vary per run (a config key; type to search). Index a second transducer as "
                "transducer(2).serial.\nOr a run column: include (skip runs), config (per-run YAML), notes.",
                [k for k in ("condition", "include", "config", "notes") if k not in self.manifest.columns] + keys,
                0, True)
            if not ok or not key:
                return
        with self._edit("Add column"):
            self.manifest.add_column(key.strip())
        self._fill_table()

    def remove_column(self) -> None:
        cols = sorted({i.column() for i in self.table.selectedIndexes()})
        with self._edit("Remove column"):
            for c in [self.manifest.columns[j] for j in cols]:
                if c != "subject_id":
                    self.manifest.remove_column(c)
        self._fill_table()

    def import_neuronav(self, summary: str | Path | None = None) -> None:
        """Runs from the Localite ingestion outputs (neuronav_gui / neuronav_ingest_batch).

        Without an argument this opens the import dialog; with a summary CSV it
        adds one run per ingested target directly."""
        if summary is not None:
            imported = from_neuronav_summary(summary)
            replace = False
        else:
            from .localite_import import LocaliteImportDialog
            dlg = LocaliteImportDialog(self.window, self)
            if dlg.exec() != dlg.DialogCode.Accepted or dlg.result_manifest is None:
                return
            imported, replace = dlg.result_manifest, dlg.replace.isChecked()
            if dlg.folder.text() and dlg.folder.text() != (self.model.get("path.localite_post") or ""):
                self.model.set("path.localite_post", dlg.folder.text())
                self.changed.emit()
        if not imported.rows:
            QMessageBox.information(self, "Import", "No ingested targets to import.")
            return
        with self._edit("Import runs from Localite"):
            if replace:
                self.manifest.rows = []
            for r in imported.rows:
                self.manifest.add_row(r)
        self._fill_table()
        self.window.statusBar().showMessage(f"Imported {plural(len(imported.rows), 'run')} from the Localite ingestion", 6000)

    def default_path(self) -> Path | None:
        """<config name>_runs.tsv next to the study config (None while the config is unsaved)."""
        cfg = self.model.path
        return cfg.parent / f"{cfg.stem}_runs.tsv" if cfg is not None else None

    def save_as(self) -> Path | None:
        self._browse()
        path = self.manifest_path()
        return self.save_manifest(path) if path is not None else None

    def save_manifest(self, path: Path | None = None) -> Path | None:
        path = path or self.manifest_path() or self.manifest.path
        if path is None:
            path = self.default_path()
            if path is None or path.exists():
                self._browse()
                path = self.manifest_path()
                if path is None:
                    return None
            else:
                self._set_manifest_path(path)
        saved = self.manifest.save(path)
        self._update_summary()
        self.window.statusBar().showMessage(f"Saved {saved}", 5000)
        return saved

    # ── MATLAB ──────────────────────────────────────────────────────────
    def _prepare(self) -> bool:
        if not self.manifest.rows:
            QMessageBox.warning(self, "Study", "The table has no runs yet: add subjects first.")
            return False
        errors = [i for i in self.manifest.issues(self.model.schema) if i.level == "error"]
        if errors:
            QMessageBox.warning(self, "Study", "Fix the table's errors first:\n" +
                                "\n".join(f"{i.where()}: {i.message}" for i in errors[:8]))
            return False
        return self.save_manifest() is not None

    def _run_page(self):
        from .main_window import RUN
        return self.window.show_special(RUN)

    def dry_run(self) -> None:
        if not self._prepare():
            return
        run = self._run_page()
        if run._ready():
            run._start("Study dry run", launcher.study_dry_run_command(
                run.matlab.text(), self.model.path, self.window.prestus_root), self.model.path.parent)

    def write_launcher(self) -> Path | None:
        if not self._prepare() or not self._run_page()._ready(need_matlab=False):
            return None
        path = launcher.write_study_launcher(self.model.path, self.window.prestus_root)
        self.window.statusBar().showMessage(f"Study launcher written: {path}", 8000)
        return path

    def run_study(self) -> None:
        path = self.write_launcher()
        if path is None:
            return
        run = self._run_page()
        if run._ready():
            run._start("Study run", launcher.run_command(run.matlab.text(), path), path.parent)

    def group_report_path(self) -> Path | None:
        """The study's group report; for a config without a study (no study.name nor
        manifest), the report over the subjects in path.sim."""
        sim = self.model.get("path.sim") or None
        name, manifest = self.model.get("study.name") or None, self.manifest_path()
        if name or manifest is not None:
            return study_report_path(sim, name, manifest)
        return subject_list_report_path(sim, self.model.get("simulation.medium") or None,
                                        self.model.get("io.output_affix") or None)

    def open_group_report(self) -> None:
        path = self.group_report_path()
        if path is not None and path.is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
        else:
            QMessageBox.information(self, "Group report", "No group report yet" +
                                    (f" ({path})" if path is not None else " (set path.sim)") +
                                    ". Build it once the runs have finished.")

    def group_report(self) -> None:
        run = self._run_page()
        if self.manifest.rows and self.save_manifest() is None:
            return
        if run._ready():
            run._start("Group report", launcher.group_report_command(
                run.matlab.text(), self.model.path, self.window.prestus_root), self.model.path.parent)
