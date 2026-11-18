"""Find a setting by its label, config key, choice labels or help text (Edit > Find setting,
Ctrl+F) and jump to it with MainWindow.goto. Only keys the window can show are listed;
a key hidden by the current settings (visible_if) is listed but cannot be opened."""

from __future__ import annotations

from dataclasses import dataclass

from PyQt6.QtCore import QEvent, Qt
from PyQt6.QtWidgets import QDialog, QLabel, QLineEdit, QListWidget, QListWidgetItem, QVBoxLayout

from .. import labels
from ..schema import Node

LIMIT = 50                  # results listed; the count says when there are more


@dataclass
class Entry:
    node: Node
    label: str
    where: str              # 'Tab › Card'
    text: str               # lower-case search text: label, key, choices, help


def entries(schema, section_of, can_open) -> list[Entry]:
    """One entry per setting the window can open. section_of(path) is the card a key is shown
    on (MainWindow.display_section); can_open(path) whether goto reaches it."""
    tab_label = {t["key"]: t.get("label", t["key"]) for t in schema.tabs}
    out = []
    for node in schema.leaves():
        if node.derived or not can_open(node.path):
            continue
        section = section_of(node.path)
        tab = tab_label.get(schema.tab_of(section) or "", "")
        card = schema.section_label(section)
        where = f"{tab} › {card}" if tab and tab != card else card
        name = labels.label(node)
        choices = " ".join(f"{c} {labels.choice(node, c)}" for c in node.choices or [])
        text = " ".join([name, node.path, choices, node.help or "", where]).lower()
        out.append(Entry(node, name, where, text))
    return out


def search(items: list[Entry], query: str) -> list[Entry]:
    """Entries containing every word of `query`, best first: words in the label count
    most, then in the key, then anywhere (choices, help, card)."""
    words = query.lower().split()
    if not words:
        return []
    scored = []
    for e in items:
        if not all(w in e.text for w in words):
            continue
        name, path = e.label.lower(), e.node.path.lower()
        score = sum(3 if w in name else 2 if w in path else 0 for w in words)
        q = query.strip().lower()
        score += 5 * (name == q or path == q or path.endswith("." + q)) + 2 * name.startswith(q)
        scored.append((-score, len(e.label), e.label, e))
    return [e for *_, e in sorted(scored, key=lambda s: s[:3])]


class FieldSearch(QDialog):
    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.setWindowTitle("Find setting")
        self.resize(560, 420)
        lay = QVBoxLayout(self)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText("Setting name, config key, or words from its help")
        self.edit.setClearButtonEnabled(True)
        self.edit.textChanged.connect(self._update)
        self.edit.installEventFilter(self)
        self.list = QListWidget()
        self.list.itemActivated.connect(self._open)
        self.count = QLabel()
        self.count.setObjectName("CardSub")
        lay.addWidget(self.edit)
        lay.addWidget(self.list, 1)
        lay.addWidget(self.count)
        self._items: list[Entry] = []

    def open_search(self) -> None:
        w = self.window
        self._items = entries(w.model.schema, w.display_section, w.can_goto)
        self.edit.selectAll()
        self._update(self.edit.text())
        self.show()
        self.raise_()
        self.activateWindow()
        self.edit.setFocus()

    def _update(self, text: str) -> None:
        self.list.clear()
        hits = search(self._items, text)
        model = self.window.model
        for e in hits[:LIMIT]:
            shown = model.is_visible(e.node.path)
            line2 = f"{e.where} · {e.node.path}" + ("" if shown else " · not used with the current settings")
            item = QListWidgetItem(f"{e.label}\n{line2}")
            item.setData(Qt.ItemDataRole.UserRole, e.node.path)
            item.setToolTip(e.node.help or "")
            if not shown:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.list.addItem(item)
        first = next((i for i in range(self.list.count()) if self.list.item(i).flags() & Qt.ItemFlag.ItemIsEnabled), -1)
        self.list.setCurrentRow(first)
        if not text.strip():
            self.count.setText("Type to search all settings.")
        else:
            more = f" (showing {LIMIT})" if len(hits) > LIMIT else ""
            self.count.setText(f"{len(hits)} setting{'s' * (len(hits) != 1)} found{more}.")

    def _open(self, item: QListWidgetItem | None) -> None:
        if item is None or not item.flags() & Qt.ItemFlag.ItemIsEnabled:
            return
        self.hide()
        self.window.goto(item.data(Qt.ItemDataRole.UserRole))

    def eventFilter(self, obj, event) -> bool:
        """Up/Down move through the results while typing; Enter opens the selected one."""
        if obj is self.edit and event.type() == QEvent.Type.KeyPress:
            key = event.key()
            if key in (Qt.Key.Key_Down, Qt.Key.Key_Up, Qt.Key.Key_PageDown, Qt.Key.Key_PageUp):
                self.list.keyPressEvent(event)
                return True
            if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                self._open(self.list.currentItem())
                return True
        return super().eventFilter(obj, event)
