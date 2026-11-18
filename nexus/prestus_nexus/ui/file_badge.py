"""A file icon that names the YAML files behind a selection on hover.

Selectors show what was chosen (a device, a driving system); the library files it resolves to
are details most users do not need in the text. A FileBadge sits next to the selector and lists
them in its tooltip (role, file name, folder), so file names are found the same way everywhere.
"""

from __future__ import annotations

from html import escape
from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QComboBox, QToolButton

from . import icons, style
from .tooltip import FAINT


class FileBadge(QToolButton):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAutoRaise(True)
        self.setIcon(icons.icon("file", style.active()["ink3"]))
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.files: list[tuple[str, str]] = []
        self.hide()

    def set_files(self, files) -> None:
        """`files`: (role, path or note) pairs, e.g. ("Device", Path(...)); empty paths are
        skipped. A str that is not a path (e.g. 'none: looked for X.yaml') is shown as is."""
        self.files = [(role, str(p)) for role, p in files if p]
        lines = []
        for role, p in self.files:
            path = Path(p)
            if path.suffix in (".yaml", ".yml") and path.parent != Path("."):
                lines.append(f"<b>{escape(role)}</b>: {escape(path.name)}"
                             f"<br><span style='color:{FAINT}'>{escape(str(path.parent))}</span>")
            else:
                lines.append(f"<b>{escape(role)}</b>: {escape(p)}")
        self.setToolTip("<br>".join(lines))
        self.setVisible(bool(lines))


def fit_combo_to_items(combo: QComboBox, min_width: int = 0) -> None:
    """Widen `combo` and its drop-down list to the longest item, so names are not cut."""
    fm = combo.fontMetrics()
    text = max((fm.horizontalAdvance(combo.itemText(i)) for i in range(combo.count())), default=0)
    combo.setMinimumWidth(max(min_width, text + 48))         # + arrow and frame
    combo.view().setMinimumWidth(text + 32)
