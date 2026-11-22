"""Help ▸ Check for updates…: compare the PRESTUS checkout with the latest release and,
when asked, move to it (update.py does the work; here it runs off the GUI thread)."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QLabel, QMessageBox, QPlainTextEdit,
                             QVBoxLayout)

from .. import update
from . import style


class _Worker(QThread):
    """Runs one update.py call; `done` carries its return value."""
    done = pyqtSignal(object)
    line = pyqtSignal(str)

    def __init__(self, fn, parent=None):
        super().__init__(parent)
        self._fn = fn

    def run(self) -> None:
        try:
            self.done.emit(self._fn(self.line.emit))
        except Exception as err:               # shown in the dialog, never lost in the thread
            self.done.emit(err)


class UpdateDialog(QDialog):
    def __init__(self, prestus_root: Path, parent=None, remote: str = update.REPO_URL):
        super().__init__(parent)
        self.setWindowTitle("PRESTUS updates")
        self.root, self.remote = Path(prestus_root), remote
        self.status: update.Status | None = None
        self.result_: update.Result | None = None
        self._worker: _Worker | None = None

        v = QVBoxLayout(self)
        self.headline = QLabel("<b>Checking for a new release…</b>")
        self.headline.setTextFormat(Qt.TextFormat.RichText)
        self.details = QLabel()
        self.details.setWordWrap(True)
        self.details.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.technical = QLabel()               # the git state: small print for maintainers
        self.technical.setWordWrap(True)
        self.technical.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.technical.setStyleSheet(f"color: {style.active()['ink3']}; font-size: 11px; "
                                     "padding-top: 8px;")
        self.log = QPlainTextEdit(readOnly=True)
        self.log.setVisible(False)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        self.update_btn = self.buttons.addButton("Update…", QDialogButtonBox.ButtonRole.ActionRole)
        self.update_btn.setEnabled(False)
        self.update_btn.clicked.connect(self.start_update)
        self.buttons.rejected.connect(self.reject)
        for w in (self.headline, self.details, self.log, self.technical, self.buttons):
            v.addWidget(w)
        self.setMinimumWidth(520)

    # ── checking ─────────────────────────────────────────────────────────
    def start_check(self) -> None:
        self._run(lambda _log: update.check(self.root, self.remote), self.show_status)

    def show_status(self, st) -> None:
        if isinstance(st, Exception):
            self.headline.setText("<b>The check failed</b>")
            self.details.setText(str(st))
            return
        self.status = st
        self.headline.setText(f"<b>{st.headline()}</b>")
        self.details.setText("\n\n".join(st.lines()[1:]))
        self.technical.setText("\n".join(t for t in (f"Folder: {st.root}", st.technical()) if t))
        self.update_btn.setText(f"Update to {st.latest}…" if st.latest else "Update…")
        self.update_btn.setEnabled(st.can_update)

    # ── updating ─────────────────────────────────────────────────────────
    def confirm(self) -> bool:
        """Ask before the checkout moves (running simulations would switch code)."""
        box = QMessageBox(self)
        box.setWindowTitle("PRESTUS updates")
        box.setIcon(QMessageBox.Icon.Question)
        box.setText(f"Update PRESTUS to {self.status.latest}?")
        box.setInformativeText(f"{update.RUNNING_JOBS_NOTE}\n\nFolder: {self.root}")
        go = box.addButton("Update", QMessageBox.ButtonRole.AcceptRole)
        box.setDefaultButton(box.addButton(QMessageBox.StandardButton.Cancel))
        box.exec()
        return box.clickedButton() is go

    def start_update(self) -> None:
        if self.status is None or not self.status.can_update or not self.confirm():
            return
        self.update_btn.setEnabled(False)
        self.headline.setText(f"<b>Updating to {self.status.latest}…</b>")
        self.log.setVisible(True)
        target = self.status.latest
        self._run(lambda log: update.apply(self.root, target, self.remote, log=log),
                  self.show_result)

    def show_result(self, res) -> None:
        if isinstance(res, Exception):
            self.headline.setText("<b>The update failed</b>")
            self.log.appendPlainText(str(res))
            return
        self.result_ = res
        if not res.ok:
            self.headline.setText("<b>PRESTUS was not updated</b>")
        elif res.changed:
            self.headline.setText(f"<b>PRESTUS is now on {res.target}: restart the Nexus</b>")
            self.details.setVisible(False)                 # the check is out of date
            self.technical.setText(f"Folder: {self.root}")
        else:
            self.headline.setText("<b>PRESTUS is up to date</b>")

    # ── plumbing ─────────────────────────────────────────────────────────
    def _run(self, fn, then) -> None:
        self._worker = _Worker(fn, self)
        self._worker.line.connect(self.log.appendPlainText)
        self._worker.done.connect(then)
        self._worker.start()

    def reject(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.wait()                # git finishes its step; never killed half-way
        super().reject()


def show_update_dialog(prestus_root: Path, parent=None) -> None:
    dialog = UpdateDialog(prestus_root, parent)
    dialog.start_check()
    dialog.exec()
