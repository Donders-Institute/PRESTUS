"""Acoustic properties card body: medium_properties as one tissue × property table.

The table shows the values this config uses (toolbox default, project defaults, study), one
column per tissue the medium models: water alone for a 'water' medium, every tissue for
'layered' and 'phantom' (a phantom models only the tissues its layers assign). It is greyed and
read-only while the card is locked; unlocked, a cell edit sets that key, and an emptied cell
resets it. Rows that matter for the thermal simulation only (unc.THERMAL) are shaded;
absorption follows attenuation (unc.PROPERTIES order). Absorption shares the attenuation power law (y = alpha_power): a tissue sets either
absorption_fraction or absorption_coeff, and the table shows the other one as implied
(thermal_absorption_properties.m). Hovering a cell shows its key, its default and the
uncertainty-variant values of that property (uncertainty.py).

For a layered medium with a segmentation, each column header's tooltip gives the tissue's volume
in it, and tissues the segmentation gives no voxels to under the current layers (check_layers.m, as
mirrored in tissue_census.py) are greyed but stay editable. The labels are counted in the
background, once per file.

The bounds plot of the uncertainty variants is on the Uncertainty sub-tab (uncertainty_view.py).
The header links to the ITRUSST tissue-property overview (unc.ITRUSST_PROPERTIES_URL).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from PyQt6.QtCore import QSize, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QBrush, QColor, QFont
from PyQt6.QtWidgets import (QAbstractItemView, QHBoxLayout, QHeaderView, QLabel, QPushButton,
                             QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .. import segmentation as sg
from .. import tissue_census as tc
from .. import uncertainty as unc
from ..checks import check_value
from ..merge import is_empty
from .tooltip import FAINT
from . import choice_icons, icons, style
from .fields import display_text
from .uncertainty_view import itrusst_button

SECTION = "medium_properties"
PROP_LABEL = {k: (label, unit) for k, label, unit in unc.PROPERTIES}
PROP_LABEL["absorption_coeff"] = ("Absorption coefficient", "dB/(cm·MHz^y), y of attenuation")
PROP_LABEL["absorption_fraction"] = ("Absorption fraction of attenuation", "")
PATH_ROLE = Qt.ItemDataRole.UserRole
_COUNT = ThreadPoolExecutor(max_workers=1, thread_name_prefix="tissue-census")


def _parse(text: str):
    """Cell text -> value: '' = reset, a number, else None (invalid)."""
    text = text.strip()
    if text in ("", "[]", "—"):
        return []
    try:
        x = float(text)
    except ValueError:
        return None
    return int(x) if x.is_integer() and "." not in text and "e" not in text.lower() else x


def _stat_key(f: Path) -> tuple | None:
    try:
        st = f.stat()
    except OSError:
        return None
    return (str(f), st.st_mtime, st.st_size)


def _num(v) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


class MediumPanel(QWidget):
    changed = pyqtSignal()
    goto = pyqtSignal(str)              # a key path to open (the note's links)

    def __init__(self, model, prestus_root: str | Path, parent=None):
        super().__init__(parent)
        self.model = model
        self.prestus_root = Path(prestus_root)
        self.locked = True
        self._filling = False
        self._sets_for: tuple | None = None
        self.sets: unc.PropertySets | None = None
        # segmentation census: which tissues get voxels (tissue_census.py)
        self.census: dict[str, tc.TissueUse] | None = None
        self.census_state = "off"       # 'off' | 'counting' | 'done' | 'failed'
        self.census_file: Path | None = None
        self._count = None              # pending (file, label_volumes future)
        self._failed_key = None
        self._census_error = ""
        self._count_timer = QTimer(self, interval=100)
        self._count_timer.timeout.connect(self._poll_count)

        node = model.schema.get(SECTION)
        self.tissues = [t for t in (node.children if node else []) if t.is_group]
        self.props: list[str] = []
        for t in self.tissues:
            for p in t.children:
                if not p.is_group and p.name not in self.props:
                    self.props.append(p.name)
        # absorption right after attenuation, as in unc.PROPERTIES (config order puts it last)
        order = [k for k, _, _ in unc.PROPERTIES]
        self.props.sort(key=lambda p: order.index(p) if p in order else len(order))
        self.thermal_rows = set(self.props) & set(unc.THERMAL)

        v = QVBoxLayout(self)
        v.setContentsMargins(0, 4, 0, 0)
        head = QHBoxLayout()
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        self.note.setTextFormat(Qt.TextFormat.RichText)
        self.note.linkActivated.connect(self.goto.emit)
        head.addStretch(1)
        self.reset_all = QPushButton("Reset all to defaults")
        self.reset_all.setProperty("compact", True)
        self.reset_all.clicked.connect(self._reset_all)
        self.open_itrusst = itrusst_button()
        head.addWidget(self.open_itrusst, 0, Qt.AlignmentFlag.AlignTop)
        head.addWidget(self.reset_all, 0, Qt.AlignmentFlag.AlignTop)
        v.addLayout(head)

        self.table = QTableWidget(len(self.props), len(self.tissues))
        self.table.setHorizontalHeaderLabels([unc.TISSUE_LABEL.get(t.name, t.name) for t in self.tissues])
        self.table.horizontalHeader().setIconSize(QSize(22, 22))
        self.table.setVerticalHeaderLabels([self._row_label(p) for p in self.props])
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.itemChanged.connect(self._item_changed)
        v.addWidget(self.table)
        v.addWidget(self.note)                  # how to read and edit the table, below it

        self.refresh()

    # ── helpers ─────────────────────────────────────────────────────────
    @staticmethod
    def _row_label(prop: str) -> str:
        label, unit = PROP_LABEL.get(prop, (prop, ""))
        return f"{label} [{unit}]" if unit else label

    def _path(self, row: int, col: int) -> str:
        return f"{SECTION}.{self.tissues[col].name}.{self.props[row]}"

    def medium(self) -> str:
        return str(self.model.get("simulation.medium") or "layered")

    def shown_tissues(self) -> list[str]:
        """Tissues of the current medium: water alone for 'water', all otherwise."""
        names = [t.name for t in self.tissues]
        return ["water"] if self.medium() == "water" and "water" in names else names

    def _header_icons(self) -> None:
        """A head-layer glyph above each tissue column, its own layer accented; with a
        segmentation census, unused tissues greyed with the reason under their name (the
        volume is in the tooltip)."""
        c = style.active()
        for col, t in enumerate(self.tissues):
            glyph = choice_icons.tissue_glyph(t.name)
            item = self.table.horizontalHeaderItem(col)
            if item is None:
                continue
            use = self.census.get(t.name) if self.census else None
            name = unc.TISSUE_LABEL.get(t.name, t.name)
            used = use is None or use.modelled
            if glyph:
                item.setIcon(icons.icon(glyph, c["ink3"] if used else c["ink4"], 20,
                                        accent=c["accent"] if used else c["ink4"]))
            if use is None:
                item.setText(name)
                item.setToolTip("")
                item.setForeground(QBrush())
                continue
            # the volume is in the tooltip; the header names only why a tissue is not modelled
            sub = {"modelled": "", "remainder": "rest", "absent": "no voxels",
                   "merged": "→ skull (pCT)", "no_layer": "no layer"}[use.state]
            # 'Skull (single layer)' and long names ('Skull trabecular') wrap to fit equal columns
            head = name.replace(" (", "\n(")
            first, _, rest = head.partition("\n")
            if len(first) > 12 and " " in first:
                first = first.replace(" ", "\n", 1)
            item.setText(first + ("\n" + rest if rest else "") + ("\n" + sub if sub else ""))
            item.setForeground(QColor(c["ink"] if used else c["ink4"]))
            state = ("Modelled" if used else "<b>Not modelled</b> for this subject")
            item.setToolTip(f"<b>{name}</b>: {state}. {use.note}"
                            + (f"<br>{tc.volume_text(use.volume_ml)} in the whole-head segmentation; the "
                               "simulation grid is a crop of it." if use.state == "modelled" else ""))

    def _unused(self, tissue: str) -> tc.TissueUse | None:
        use = self.census.get(tissue) if self.census else None
        return use if use is not None and not use.modelled else None

    def _fit_height(self) -> None:
        header = self.table.horizontalHeader()
        # the size hint follows header text changed this refresh (height() only after a layout pass)
        h = max(header.height(), header.sizeHint().height()) + 2 * self.table.frameWidth()
        h += sum(self.table.rowHeight(r) for r in range(self.table.rowCount()))
        self.table.setFixedHeight(h + 2)

    def _bounds_tip(self, tissue: str, prop: str) -> str:
        if self.sets is None:
            return ""
        lines = []
        for variant in unc.VARIANTS:
            if variant not in self.sets.sets:
                continue
            if variant in unc.PAIRS["trans"] and prop in unc.THERMAL:
                continue
            if self.sets.inherited(variant, tissue, prop):
                continue
            x = self.sets.value(variant, tissue, prop)
            if isinstance(x, (int, float)):
                lines.append(f"{unc.VARIANT_SHORT[variant]}: {x:.4g}")
        return ("<br>Uncertainty variants: " + ", ".join(lines)) if lines else ""

    def _implied(self, tissue: str, prop: str) -> float | None:
        """The absorption field a tissue does not set, from the one it sets (same power law)."""
        base = f"{SECTION}.{tissue}."
        alpha = _num(self.model.get(base + "alpha_coeff"))
        coeff = self.model.get(base + "absorption_coeff")
        if prop == "absorption_coeff" and is_empty(coeff):
            frac = _num(self.model.get(base + "absorption_fraction"))
            frac = 1.0 if frac is None else frac
            return frac * alpha if alpha is not None else None
        if prop == "absorption_fraction" and not is_empty(coeff):
            c = _num(coeff)
            return c / alpha if c is not None and alpha else None
        return None

    # ── state ───────────────────────────────────────────────────────────
    def set_locked(self, locked: bool) -> None:
        if locked != self.locked:
            self.locked = locked
            self.refresh()

    def refresh(self) -> None:
        self._refresh_sets()
        self._refresh_census()
        self._header_icons()                 # follows the theme
        c = style.active()
        shown = self.shown_tissues()
        thermal_bg = QColor(c["thermal"])
        thermal_bg.setAlpha(28)
        self._filling = True
        n_set = 0
        for col, tissue in enumerate(self.tissues):
            self.table.setColumnHidden(col, tissue.name not in shown)
            own = {p.name: p for p in tissue.children if not p.is_group}
            for row, prop in enumerate(self.props):
                node = own.get(prop)
                item = self.table.item(row, col)
                if item is None:
                    item = QTableWidgetItem()
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                    self.table.setItem(row, col, item)
                item.setBackground(QBrush(thermal_bg) if prop in self.thermal_rows else QBrush())
                if node is None:
                    item.setText("")
                    item.setFlags(Qt.ItemFlag.NoItemFlags)
                    continue
                path = node.path
                item.setData(PATH_ROLE, path)
                value = self.model.get(path)
                implied = self._implied(tissue.name, prop)
                font = QFont(item.font())
                modified = self.model.is_modified(path)
                n_set += self.model.is_overridden(path)
                font.setBold(modified)
                font.setItalic(implied is not None)
                item.setFont(font)
                if implied is not None:
                    item.setText(f"({implied:.4g})")
                else:
                    item.setText(display_text(value) if not is_empty(value) else "—")
                err = check_value(node, value)
                unused = self._unused(tissue.name)
                if self.locked or unused is not None:
                    colour = c["ink4"]                  # greyed: read-only, or not modelled
                elif err:
                    colour = c["red"]
                elif implied is not None:
                    colour = c["ink3"]
                else:
                    colour = c["accent"] if modified else c["ink"]
                item.setForeground(QColor(colour))
                flags = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                if not self.locked:
                    flags |= Qt.ItemFlag.ItemIsEditable
                item.setFlags(flags)
                extra = ""
                if implied is not None:
                    extra = ("<br><i>Implied by the absorption coefficient: coefficient / attenuation "
                             "coefficient.</i>" if prop == "absorption_fraction" else
                             "<br><i>Implied by the absorption fraction: fraction × attenuation "
                             "coefficient, with the attenuation power law y.</i>")
                    extra += "<br>A tissue sets one of the two; typing a value here switches to it."
                elif prop == "absorption_coeff":
                    extra = "<br>Absorption follows the attenuation power law (y = alpha_power)."
                if prop in self.thermal_rows:
                    extra += "<br>Thermal simulation only."
                if unused is not None:
                    extra += f"<br><i>Not modelled for this subject: {unused.note}</i>"
                tip = (f"<b>{unc.TISSUE_LABEL.get(tissue.name, tissue.name)}: {self._row_label(prop)}</b>"
                       f"<br>Default: {display_text(self.model.baseline(path)) or '(empty)'}"
                       + (f"<br><span style='color:{c['red']}'>{err}</span>" if err else "")
                       + extra + self._bounds_tip(tissue.name, prop)
                       + f"<br><span style='color:{FAINT}'>Key: {path}</span>")
                item.setToolTip(tip)
        self._filling = False
        self._fit_height()
        self.reset_all.setVisible(n_set > 0)
        self.reset_all.setEnabled(not self.locked)
        legend = (f"<span style='background:{thermal_bg.name(QColor.NameFormat.HexArgb)}'>"
                  "&nbsp;Shaded&nbsp;</span>: thermal only · (brackets): implied"
                  + (" · bold: changed" if n_set else ""))
        lock = ("Locked: show advanced options and unlock to edit." if self.locked else
                "Edit a cell to change it; empty it to reset.")
        phantom = (" Phantom: only tissues its layers assign are modelled."
                   if self.medium() == "phantom" else "")
        bounds = ""
        if self.medium() != "water" and self.model.get("simulation.uncertainty") is not True:
            bounds = " <a href='simulation.uncertainty'>Uncertainty analysis</a> plots the ITRUSST ranges."
        self.note.setText(lock + " " + legend + "." + phantom + bounds + self._census_note())

    # ── segmentation census ─────────────────────────────────────────────
    def _census_source(self) -> Path | None:
        """final_tissues.nii.gz of a layered medium, if it is there."""
        if self.medium() != "layered":
            return None
        f = tc.seg_file(sg.m2m_path(self.model))
        return f if f is not None and f.is_file() else None

    def _refresh_census(self) -> None:
        """Count the segmentation's labels in the background (once per file), then map them
        through the current layers (cheap, on every refresh)."""
        f = self._census_source()
        self.census_file = f
        if f is None:
            self.census, self.census_state = None, "off"
            return
        volumes = tc.cached_volumes(f)
        if volumes is None:
            self.census = None
            if self.census_state == "failed" and self._failed_key == _stat_key(f):
                return                          # not again until the file changes
            if self._count is None:
                self.census_state = "counting"
                self._count = (f, _COUNT.submit(tc.label_volumes, f))
                self._count_timer.start()
            return
        eff = self.model.effective()
        props, layers = eff.get("medium_properties"), eff.get("layers")
        names = list(props) if isinstance(props, dict) else [t.name for t in self.tissues]
        pct = self.model.get("pct.enabled") in (1, True)
        self.census = tc.tissue_use(layers if isinstance(layers, dict) else {}, names, volumes, pct)
        self.census_state = "done"

    def _poll_count(self) -> None:
        if self._count is None or not self._count[1].done():
            return
        self._count_timer.stop()
        f, future = self._count
        self._count = None
        try:
            future.result()
        except Exception as e:                  # unreadable file: say so, keep the table as is
            self.census, self.census_state = None, "failed"
            self._census_error = f"{type(e).__name__}: {e}"
            self._failed_key = _stat_key(f)
            self.note.setText(self.note.text() + self._census_note())
            return
        self.refresh()

    def wait_census(self, timeout: float = 30.0) -> None:
        """Finish a pending count now (tests and offscreen use)."""
        if self._count is not None:
            self._count[1].result(timeout=timeout)
            self._poll_count()

    def _census_note(self) -> str:
        c = style.active()
        name = self.census_file.parent.name if self.census_file else ""
        if self.census_state == "counting":
            return f"<br>Counting the tissue labels of {name}…"
        if self.census_state == "failed":
            return (f"<br><span style='color:{c['amber']}'>Could not read {name}/{tc.SEG_FILE} "
                    f"({self._census_error}).</span>")
        if self.census_state != "done" or not self.census:
            return ""
        unused = [unc.TISSUE_LABEL.get(t, t) for t, u in self.census.items() if not u.modelled]
        text = f"<br>{name}: whole-head volumes in the column tooltips."
        if unused:
            text += f" Greyed, not modelled: {', '.join(unused)}."
        return text

    def _sets_shown(self) -> bool:
        """Variant values in the cell tooltips: uncertainty analysis on and a medium with brain and skull."""
        return self.medium() != "water" and self.model.get("simulation.uncertainty") is True

    def _refresh_sets(self) -> None:
        if not self._sets_shown():
            self.sets = None
            return
        sig = unc.signature(self.model, self.prestus_root)
        if sig == self._sets_for and self.sets is not None:
            return
        self._sets_for = sig
        self.sets = unc.property_sets(self.model, self.prestus_root)

    # ── editing ─────────────────────────────────────────────────────────
    def _item_changed(self, item: QTableWidgetItem) -> None:
        if self._filling or self.locked:
            return
        path = item.data(PATH_ROLE)
        if not path:
            return
        value = _parse(item.text().strip("()"))
        if value is None:                       # not a number: show the stored value again
            self.refresh()
            return
        if value == []:
            self.model.reset(path)
        else:
            self.model.set(path, value)
            self._one_absorption_field(path)
        self.refresh()
        self.changed.emit()

    def _one_absorption_field(self, path: str) -> None:
        """A tissue sets absorption_fraction or absorption_coeff, never both
        (thermal_absorption_model.m): setting one clears the other in this config. A fraction
        inherited from the defaults is cleared by PRESTUS (medium_absorption_layer.m)."""
        tissue_path, prop = path.rsplit(".", 1)
        if prop == "absorption_fraction":
            other = f"{tissue_path}.absorption_coeff"
            if not is_empty(self.model.get(other)):
                self.model.set(other, [])
        elif prop == "absorption_coeff":
            other = f"{tissue_path}.absorption_fraction"
            if self.model.is_overridden(other):
                self.model.reset(other)

    def _reset_all(self) -> None:
        if self.locked:
            return
        for col in range(len(self.tissues)):
            for row in range(len(self.props)):
                item = self.table.item(row, col)
                path = item.data(PATH_ROLE) if item is not None else None
                if path and self.model.is_overridden(path):
                    self.model.reset(path)
        self.refresh()
        self.changed.emit()

    def focus(self, path: str) -> bool:
        for col in range(len(self.tissues)):
            for row in range(len(self.props)):
                if self._path(row, col) == path:
                    self.table.setCurrentCell(row, col)
                    self.table.setFocus(Qt.FocusReason.OtherFocusReason)
                    return True
        return False
