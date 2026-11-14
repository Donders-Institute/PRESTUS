"""Dialog: study manifest rows from the Localite ingestion outputs of neuronav_gui."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog,
                             QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem, QPushButton,
                             QRadioButton, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)

from .. import neuronav_import as ni
from ..study import Manifest

HELP = ("Reads what the Localite page (or <tt>neuronav_gui</tt>) writes under <tt>path.localite_post</tt>: "
        "<tt>neuronav_ingest_summary.csv</tt> (batch window) and the per-session "
        "<tt>InstrumentMarker_*.json</tt> files (“Write outputs”). Each run uses "
        "<tt>placement.mode: localite</tt>, so PRESTUS reads the ingested markers itself when it runs. "
        "Ingest (or re-ingest) sessions on the <b>Localite page</b>, then <b>Scan</b>.")


class LocaliteImportDialog(QDialog):
    def __init__(self, window, parent: QWidget | None = None):
        super().__init__(parent)
        self.window = window
        self.model = window.model
        self.entries: list[ni.Entry] = []
        self.result_manifest: Manifest | None = None
        self.setWindowTitle("Import runs from Localite ingestion")
        self.resize(820, 620)
        v = QVBoxLayout(self)
        intro = QLabel(HELP)
        intro.setWordWrap(True)
        v.addWidget(intro)

        row = QHBoxLayout()
        row.addWidget(QLabel("Ingested markers folder"))
        self.folder = QLineEdit(str(self.model.get("path.localite_post") or ""))
        row.addWidget(self.folder, 1)
        b = QToolButton(text="…")
        b.clicked.connect(self._browse)
        row.addWidget(b)
        scan = QPushButton("Scan")
        scan.clicked.connect(lambda: self.scan())
        row.addWidget(scan)
        gui = QPushButton("Localite page")
        gui.setToolTip("Detect, name and write sessions (Placement tab ▸ Localite)")
        gui.clicked.connect(self._localite_page)
        row.addWidget(gui)
        v.addLayout(row)

        filt = QHBoxLayout()
        self.targets = QListWidget()
        self.sessions = QListWidget()
        for title, lw in (("Targets", self.targets), ("Sessions", self.sessions)):
            box = QVBoxLayout()
            box.addWidget(QLabel(f"<b>{title}</b>"))
            lw.itemChanged.connect(lambda _i: self._preview())
            box.addWidget(lw)
            filt.addLayout(box)
        opts = QVBoxLayout()
        opts.addWidget(QLabel("<b>Runs</b>"))
        self.g_target = QRadioButton("One run per target (single transducer)")
        self.g_session = QRadioButton("One run per session: targets → transducers 1, 2, …")
        n_tx = self.model.n_transducers
        (self.g_session if n_tx > 1 else self.g_target).setChecked(True)
        grp = QButtonGroup(self)
        for rb in (self.g_target, self.g_session):
            grp.addButton(rb)
            rb.toggled.connect(lambda _c: self._preview())
            opts.addWidget(rb)
        opts.addWidget(QLabel("Transducer order"))
        self.order = QComboBox()
        self.order.addItems(["Left (L…) first, then right (R…)", "Alphabetical", "As listed below"])
        self.order.currentIndexChanged.connect(lambda _i: self._preview())
        opts.addWidget(self.order)
        self.custom = QLineEdit()
        self.custom.setPlaceholderText("e.g. L_M1, R_M1")
        self.custom.textChanged.connect(lambda _t: self._preview())
        opts.addWidget(self.custom)
        self.check_n = QCheckBox(f"Require {n_tx} target(s) per session (the configured transducers)")
        self.check_n.setChecked(True)
        self.check_n.toggled.connect(lambda _c: self._preview())
        opts.addWidget(self.check_n)
        self.replace = QCheckBox("Replace the current manifest rows")
        opts.addWidget(self.replace)
        opts.addStretch(1)
        filt.addLayout(opts)
        v.addLayout(filt, 1)

        v.addWidget(QLabel("<b>Ingested series (QC)</b>"))
        self.qc_table = QTableWidget()
        v.addWidget(self.qc_table, 1)

        v.addWidget(QLabel("<b>Preview</b>"))
        self.table = QTableWidget()
        v.addWidget(self.table, 2)
        self.notes = QLabel()
        self.notes.setWordWrap(True)
        self.notes.setObjectName("CardSub")
        v.addWidget(self.notes)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.button(QDialogButtonBox.StandardButton.Ok).setText("Add runs")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        if self.folder.text():
            self.scan()

    def _localite_page(self) -> None:
        from .main_window import LOCALITE
        self.reject()
        self.window.show_special(LOCALITE)

    def _browse(self) -> None:
        d = QFileDialog.getExistingDirectory(self, "Ingested Localite markers", self.folder.text())
        if d:
            self.folder.setText(d)
            self.scan()

    def scan(self) -> None:
        folder = Path(self.folder.text())
        self.entries = ni.scan(folder) if folder.is_dir() else []
        for lw, values in ((self.targets, sorted({e.target for e in self.entries})),
                           (self.sessions, sorted({e.ses for e in self.entries}))):
            lw.blockSignals(True)
            lw.clear()
            for val in values:
                it = QListWidgetItem(val)
                it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                it.setCheckState(Qt.CheckState.Checked)
                lw.addItem(it)
            lw.blockSignals(False)
        if not self.custom.text():
            self.custom.setText(", ".join(ni.order_targets([e.target for e in self.entries], "lr")))
        self._preview()

    def _checked(self, lw: QListWidget) -> set[str]:
        return {lw.item(i).text() for i in range(lw.count())
                if lw.item(i).checkState() == Qt.CheckState.Checked}

    def _filtered_entries(self) -> list[ni.Entry]:
        tg, ss = self._checked(self.targets), self._checked(self.sessions)
        return [e for e in self.entries if e.target in tg and e.ses in ss]

    def build(self) -> tuple[Manifest, list[str]]:
        entries = self._filtered_entries()
        order: str | list[str] = ["lr", "alphabetical", "custom"][self.order.currentIndex()]
        if order == "custom":
            order = [t.strip() for t in self.custom.text().split(",") if t.strip()]
        grouping = "session" if self.g_session.isChecked() else "target"
        n = self.model.n_transducers if (grouping == "session" and self.check_n.isChecked()) else None
        return ni.build_manifest(entries, grouping, order, n)

    _QC_COLS = ["Subject", "Session", "Target", "Pos dev (mm)", "Rot dev (rad)",
               "Off axis (mm)", "Depth (mm)", "Pulses", "Source"]

    def _fill_qc_table(self) -> None:
        entries = self._filtered_entries()
        self.qc_table.setColumnCount(len(self._QC_COLS))
        self.qc_table.setHorizontalHeaderLabels(self._QC_COLS)
        self.qc_table.setRowCount(len(entries))

        def cell(v: float | None, fmt: str = "{:.2f}") -> str:
            return "" if v is None else fmt.format(v)

        for i, e in enumerate(entries):
            row = [e.sub, e.ses, e.target, cell(e.pos_dev_mm), cell(e.rot_dev_rad, "{:.3f}"),
                  cell(e.axis_dist_mm, "{:.1f}"), cell(e.depth_mm, "{:.1f}"),
                  cell(e.n_pulses, "{:.0f}"), e.source]
            for j, v in enumerate(row):
                self.qc_table.setItem(i, j, QTableWidgetItem(v))
        self.qc_table.resizeColumnsToContents()

    def _preview(self) -> None:
        self.custom.setEnabled(self.order.currentIndex() == 2)
        self._fill_qc_table()
        m, notes = self.build()
        self.table.setColumnCount(len(m.columns))
        self.table.setHorizontalHeaderLabels(m.columns)
        self.table.setRowCount(len(m.rows))
        for i, r in enumerate(m.rows):
            for j, c in enumerate(m.columns):
                self.table.setItem(i, j, QTableWidgetItem(r.get(c, "")))
        self.table.resizeColumnsToContents()
        n_json = sum(1 for e in self.entries if e.source == "json")
        head = (f"{len(self.entries)} ingested target(s) found"
                + (f" ({n_json} from InstrumentMarker JSON files)" if n_json else "")
                + f" → {len(m.rows)} run(s).") if self.entries else "Nothing found in this folder yet."
        self.notes.setText("<br>".join([head] + notes))
        self.result_manifest = m
