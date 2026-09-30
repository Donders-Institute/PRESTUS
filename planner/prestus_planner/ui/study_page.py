"""Study page: edit the study manifest (one row per run) and start study runs / the group report."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QFileDialog, QFrame, QHBoxLayout, QInputDialog, QLabel,
                             QLineEdit, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem,
                             QToolButton, QVBoxLayout, QWidget)

from .. import launcher
from ..study import Manifest, from_neuronav_summary, resolve_manifest

INTRO = ("A study runs PRESTUS for many subjects, each with its own settings. The manifest has one "
         "row per run: <tt>subject_id</tt>, an optional <tt>condition</tt> (e.g. the target; it becomes "
         "the run's output affix and the grouping in the group report), <tt>include</tt>, an optional "
         "per-row <tt>config</tt> YAML, and any config key as a column (e.g. "
         "<tt>placement.localite.target</tt>, <tt>transducer.target_isppa_wcm2</tt>). Empty cells keep "
         "the value of this config. MATLAB: <tt>prestus_study_start(config)</tt>.")


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
        title = QLabel("Study manifest")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        intro = QLabel(INTRO)
        intro.setObjectName("CardSub")
        intro.setWordWrap(True)
        v.addWidget(intro)

        row = QHBoxLayout()
        row.addWidget(QLabel("Manifest file"))
        self.path = QLineEdit()
        self.path.setPlaceholderText("e.g. study.tsv (relative to the study config)")
        self.path.editingFinished.connect(self._path_edited)
        browse = QToolButton(text="…")
        browse.clicked.connect(self._browse)
        row.addWidget(self.path, 1)
        row.addWidget(browse)
        v.addLayout(row)

        edit = QHBoxLayout()
        for text, slot in [("Add row", lambda: self.add_row()), ("Remove rows", self.remove_rows),
                           ("Add column…", lambda: self.add_column()),
                           ("Remove column", self.remove_column),
                           ("Import from Localite…", lambda: self.import_neuronav()),
                           ("Open neuronav_gui", lambda: self.window.open_neuronav_gui()),
                           ("Save manifest", lambda: self.save_manifest())]:
            b = QPushButton(text)
            b.clicked.connect(slot)
            edit.addWidget(b)
        edit.addStretch(1)
        v.addLayout(edit)

        self.table = QTableWidget()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.table.itemChanged.connect(self._cell_edited)
        self.table.setMinimumHeight(220)
        v.addWidget(self.table, 1)

        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self.summary)

        run = QHBoxLayout()
        self.b_dry = QPushButton("Check all rows in MATLAB (dry run)")
        self.b_dry.clicked.connect(self.dry_run)
        self.b_launcher = QPushButton("Write study launcher")
        self.b_launcher.clicked.connect(self.write_launcher)
        self.b_run = QPushButton("Run study locally")
        self.b_run.clicked.connect(self.run_study)
        self.b_report = QPushButton("Build group report")
        self.b_report.setToolTip("Group MNI maps + group report from the existing outputs")
        self.b_report.clicked.connect(self.group_report)
        for b in (self.b_dry, self.b_launcher, self.b_run, self.b_report):
            run.addWidget(b)
        run.addStretch(1)
        v.addLayout(run)

    # ── state ───────────────────────────────────────────────────────────
    @property
    def model(self):
        return self.window.model

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
        self.table.setRowCount(len(m.rows))
        for i, r in enumerate(m.rows):
            for j, c in enumerate(m.columns):
                self.table.setItem(i, j, QTableWidgetItem(r.get(c, "")))
        self.table.resizeColumnsToContents()
        self._loading = False
        self._update_summary()

    def _update_summary(self) -> None:
        m = self.manifest
        affix = self.model.get("io.output_affix") or ""
        issues = m.issues(self.model.schema, affix if isinstance(affix, str) else "")
        n_inc = sum(1 for r in m.rows if r.get("include", "").strip().lower() not in {"0", "false", "no", "n"})
        conds = [c or "(none)" for c in m.conditions()]
        text = (f"<b>{len(m.rows)}</b> row(s), <b>{n_inc}</b> included · "
                f"{len(m.subjects())} subject(s) · conditions: {', '.join(conds) or '—'}")
        if m.path:
            text += f"<br><b>File:</b> {m.path}"
        for iss in issues:
            colour = "#c9302c" if iss.level == "error" else "#b8860b"
            text += f"<br><span style='color:{colour}'>{iss.level}: {iss.where()}: {iss.message}</span>"
        self.summary.setText(text)

    def _cell_edited(self, item: QTableWidgetItem) -> None:
        if self._loading:
            return
        col = self.manifest.columns[item.column()]
        self.manifest.rows[item.row()][col] = item.text().strip()
        self._update_summary()

    def _path_edited(self) -> None:
        text = self.path.text().strip()
        if text == (self.model.get("study.manifest") or ""):
            return
        if text:
            self.model.set("study.manifest", text)
        else:
            self.model.reset("study.manifest")
        self.changed.emit()
        mp = self.manifest_path()
        self.manifest = Manifest.load(mp) if mp is not None and mp.is_file() else Manifest(path=mp)
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
        start = str(self.manifest_path() or (self.model.path.parent if self.model.path else Path.home()))
        chosen, _ = QFileDialog.getSaveFileName(self, "Study manifest", start,
                                                "Manifest (*.tsv *.csv)",
                                                options=QFileDialog.Option.DontConfirmOverwrite)
        if chosen:
            self._set_manifest_path(Path(chosen))

    # ── editing ─────────────────────────────────────────────────────────
    def add_row(self, values: dict[str, str] | None = None) -> None:
        self.manifest.add_row(values)
        self._fill_table()

    def remove_rows(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            del self.manifest.rows[r]
        self._fill_table()

    def add_column(self, key: str | None = None) -> None:
        if key is None:
            keys = sorted(n.path for n in self.model.schema.leaves()
                          if not n.path.startswith(("study.", "path.")) and n.path != "subject_id")
            key, ok = QInputDialog.getItem(self, "Add column", "Config key (dotted; index transducers "
                                           "as transducer(2).serial):", ["condition", "include", "config"] + keys,
                                           0, True)
            if not ok or not key:
                return
        self.manifest.add_column(key.strip())
        self._fill_table()

    def remove_column(self) -> None:
        cols = sorted({i.column() for i in self.table.selectedIndexes()})
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
        if replace:
            self.manifest.rows = []
        for r in imported.rows:
            self.manifest.add_row(r)
        self._fill_table()
        self.window.statusBar().showMessage(f"Imported {len(imported.rows)} run(s) from the Localite ingestion", 6000)

    def save_manifest(self, path: Path | None = None) -> Path | None:
        path = path or self.manifest_path() or self.manifest.path
        if path is None:
            self._browse()
            path = self.manifest_path()
            if path is None:
                return None
        saved = self.manifest.save(path)
        self._update_summary()
        self.window.statusBar().showMessage(f"Saved {saved}", 5000)
        return saved

    # ── MATLAB ──────────────────────────────────────────────────────────
    def _prepare(self) -> bool:
        if not self.manifest.rows:
            QMessageBox.warning(self, "Study", "The manifest has no rows.")
            return False
        errors = [i for i in self.manifest.issues(self.model.schema) if i.level == "error"]
        if errors:
            QMessageBox.warning(self, "Study", "Fix the manifest errors first:\n" +
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

    def group_report(self) -> None:
        run = self._run_page()
        if self.manifest.rows and self.save_manifest() is None:
            return
        if run._ready():
            run._start("Group report", launcher.group_report_command(
                run.matlab.text(), self.model.path, self.window.prestus_root), self.model.path.parent)
