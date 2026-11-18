"""Medium › Uncertainty card, shown while simulation.uncertainty is on: the ITRUSST bounds plot,
then a closed 'Values per variant' fold with the default medium properties next to the four
variants (Temp-H/L bound heating, Trans-H/L bound intensity). The table is locked like the
acoustic property table; unlocked, it edits simulation.uncertainty_overrides, and edited values
are bold red with a warning above the plot. draw_ranges plots the range of each pair for brain and skull
(one column per property with its own scale, as in the ITRUSST tissue-property overview), with a
line at the value this config uses. The Trans pair runs acoustic only, so thermal properties
show the Temp pair alone."""

from __future__ import annotations

import copy
from pathlib import Path

from PyQt6.QtCore import Qt, QUrl, pyqtSignal
from PyQt6.QtGui import QColor, QDesktopServices, QFont
from PyQt6.QtWidgets import (QAbstractItemView, QFrame, QHBoxLayout, QHeaderView, QLabel, QMessageBox,
                             QPushButton, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout,
                             QWidget)

from .. import uncertainty as unc
from . import icons, style
from .tooltip import FAINT

try:
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    HAVE_MPL = True
except ImportError:          # pragma: no cover
    HAVE_MPL = False

SKULL_TICK = {"skull": "S", "skull_trabecular": "T", "skull_cortical": "C"}


def _fmt(v) -> str:
    if v is None or v == [] or v == "":
        return "—"
    try:
        return f"{float(v):.4g}"
    except (TypeError, ValueError):
        return str(v)


def _colours(c: dict) -> dict:
    # report convention: thermal content warm, acoustic content cool
    return {"default": c["ink"], "temp": c["thermal"], "trans": c["accent"]}


PAIR_DX = {"temp": -0.17, "trans": 0.17}        # x offset of each pair within a tissue column
BOUND_MARKER = {"high": "^", "low": "v"}


def pair_values(sets: unc.PropertySets, pair: str, tissue: str, key: str) -> dict[str, tuple[float, bool]]:
    """{'high'|'low': (value, inherited)} for a pair; empty for the trans pair on a thermal
    property (acoustic only, so the inherited default is no bound) or when neither file sets it."""
    if pair == "trans" and key in unc.THERMAL:
        return {}
    out = {}
    for bound, variant in zip(("high", "low"), unc.PAIRS[pair]):
        x = sets.value(variant, tissue, key) if variant in sets.sets else None
        if isinstance(x, (int, float)):
            out[bound] = (float(x), sets.inherited(variant, tissue, key))
    if all(inh for _, inh in out.values()):
        return {}
    return out


def any_inherited(sets: unc.PropertySets, tissues) -> bool:
    """A plotted variant value is the inherited default (drawn hollow)."""
    return any(inh for key in unc.PLOTTED for t in tissues for pair in unc.PAIRS
               for _, inh in pair_values(sets, pair, t, key).values())


def draw_ranges(fig, sets: unc.PropertySets, tissues: tuple[str, ...], c: dict,
                legend: bool = True, hollow: bool | None = None) -> None:
    """One column per property, each on its own scale. Per tissue, the Temp pair (left) and the
    Trans pair (right) each draw a bar between their two values, with a marker per variant
    (hollow where inherited from the default); a tick marks the default. Thermal properties
    show the Temp pair only. `legend`: draw the legend (the first of the stacked plots only);
    `hollow`: list the hollow marker (default: when this plot draws one)."""
    fig.clear()
    fig.set_facecolor(c["surface"])
    col = _colours(c)
    props = [p for p in unc.PROPERTIES if p[0] in unc.PLOTTED]
    axes = fig.subplots(1, len(props))
    any_hollow = False
    for ax, (key, label, unit) in zip(axes, props):
        ax.set_facecolor(c["surface"])
        values = []
        for gi, tissue in enumerate(tissues):
            for pair in unc.PAIRS:
                v = pair_values(sets, pair, tissue, key)
                if not v:
                    continue
                x = gi + PAIR_DX[pair]
                values += [val for val, _ in v.values()]
                if len(v) == 2:
                    ax.vlines(x, min(val for val, _ in v.values()), max(val for val, _ in v.values()),
                              color=col[pair], linewidth=2.5, alpha=0.45, zorder=1)
                for bound, (val, inh) in v.items():
                    any_hollow = any_hollow or inh
                    ax.plot(x, val, BOUND_MARKER[bound], color=col[pair], markersize=5, zorder=3,
                            markerfacecolor=c["surface"] if inh else col[pair])
            d = sets.value("default", tissue, key)
            if isinstance(d, (int, float)):
                values.append(float(d))
                ax.hlines(float(d), gi - 0.34, gi + 0.34, color=col["default"], linewidth=1.2, zorder=2)
        if values:
            lo, hi = min(values), max(values)
            pad = (hi - lo) * 0.12 or abs(hi) * 0.1 or 1
            ax.set_ylim(lo - pad, hi + pad)
        if key in ("attenuation", "absorption"):
            unit = f"{unit} @ {sets.freq_mhz:g} MHz"
        ax.set_title(f"{label}\n[{unit}]" if unit else label, fontsize=style.PLOT_PT["normal"], color=c["ink2"])
        ax.set_xlim(-0.6, len(tissues) - 0.4)
        if len(tissues) > 1:
            ax.set_xticks(range(len(tissues)), [SKULL_TICK.get(t, t) for t in tissues], fontsize=style.PLOT_PT["small"])
        else:
            ax.set_xticks([])
        ax.tick_params(labelsize=style.PLOT_PT["small"], colors=c["ink3"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(c["line2"])
    # one entry per run: the default and the four variants (colour = pair, marker = H/L);
    # the hollow-marker entry only when a hollow marker is drawn
    def marker(pair, bound, label):
        return Line2D([], [], color=col[pair], marker=BOUND_MARKER[bound], linestyle="none",
                      markersize=6, label=label)
    handles = [Line2D([], [], color=col["default"], linewidth=1.2, label="Default"),
               marker("temp", "high", "Temp-H"), marker("temp", "low", "Temp-L"),
               marker("trans", "high", "Trans-H"), marker("trans", "low", "Trans-L")]
    if not legend:
        fig.subplots_adjust(left=0.05, right=0.96, top=0.82, bottom=0.12, wspace=0.9)
        return
    if any_hollow if hollow is None else hollow:
        handles.append(Line2D([], [], color=c["ink3"], marker="o", markerfacecolor=c["surface"],
                              linestyle="none", markersize=5, label="hollow: default value used"))
    leg = fig.legend(handles=handles, loc="upper center", ncol=len(handles), fontsize=style.PLOT_PT["small"], frameon=False,
                     bbox_to_anchor=(0.5, 1.0))
    for text in leg.get_texts():
        text.set_color(c["ink2"])
    fig.subplots_adjust(left=0.05, right=0.96, top=0.68, bottom=0.12, wspace=0.9)


def itrusst_button() -> QPushButton:
    """Opens the ITRUSST tissue-property overview (reference values and bounds) in the browser."""
    b = QPushButton("Open ITRUSST tissue properties ↗")
    b.setProperty("compact", True)
    b.setToolTip(f"Reference acoustic and thermal properties and their bounds: {unc.ITRUSST_PROPERTIES_URL}")
    b.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(unc.ITRUSST_PROPERTIES_URL)))
    return b


CELL_ROLE = Qt.ItemDataRole.UserRole          # (variant, tissue, prop) of an editable cell


class UncertaintyCard(QFrame):
    """Bounds plot of the Temp-H/L, Trans-H/L medium properties; the values in a fold below.
    The table is locked (greyed, read-only) like the acoustic property table; unlocked, a variant
    cell edit sets simulation.uncertainty_overrides.<variant>.<tissue>.<prop>, and an emptied cell
    (or the shipped bound typed back) removes the edit. Edited values are bold red, and a warning
    above the plot counts them. Hidden unless simulation.uncertainty."""

    changed = pyqtSignal()

    def __init__(self, model, prestus_root: str | Path, parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.prestus_root = Path(prestus_root)
        self.locked = True
        self._filling = False
        self._shown_for: tuple | None = None
        self.sets: unc.PropertySets | None = None

        v = QVBoxLayout(self)
        v.setContentsMargins(14, 10, 14, 12)
        title = QLabel("Uncertainty variants: medium properties")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        self.errors = QLabel()                  # unreadable bounds files: shown with the fold closed
        self.errors.setObjectName("CardSub")
        self.errors.setWordWrap(True)
        self.errors.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.errors)
        self.warning = QLabel()                 # edited bounds: tinted box, shown with the fold closed
        self.warning.setWordWrap(True)
        self.warning.setTextFormat(Qt.TextFormat.RichText)
        self.warning.linkActivated.connect(self._link)
        v.addWidget(self.warning)

        self.canvases: dict[str, object] = {}
        if HAVE_MPL:
            from .mpl_theme import ThemedCanvas      # drawn in the Nexus theme
            self.caption = QLabel()
            self.caption.setObjectName("CardSub")
            self.caption.setWordWrap(True)
            v.addWidget(self.caption)
            for name, heading in (("brain", "Brain"),
                                  ("skull", "Skull: S single layer · T trabecular · C cortical")):
                h = QLabel(heading)
                h.setObjectName("SubGroupTitle")
                v.addWidget(h)
                canvas = ThemedCanvas(Figure(figsize=(9, 2.6)), export_name=f"property_bounds_{name}")
                canvas.setMinimumHeight(230)
                v.addWidget(canvas)
                self.canvases[name] = canvas
        # below the plot, folded: the values behind it
        self.b_values, self.values = self._fold("Values per variant")
        v.addWidget(self.b_values)
        v.addWidget(self.values)
        vv = self.values.layout()
        head = QHBoxLayout()
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        self.note.setTextFormat(Qt.TextFormat.RichText)
        head.addWidget(self.note, 1)
        self.reset_all = QPushButton("Reset edits to ITRUSST bounds")
        self.reset_all.setProperty("compact", True)
        self.reset_all.clicked.connect(self._reset_all)
        head.addWidget(self.reset_all, 0, Qt.AlignmentFlag.AlignTop)
        self.lock_button = QToolButton()
        self.lock_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.lock_button.clicked.connect(self._toggle_lock)
        head.addWidget(self.lock_button, 0, Qt.AlignmentFlag.AlignTop)
        vv.addLayout(head)
        self.table = QTableWidget(0, 4 + len(unc.VARIANTS))
        self.table.setHorizontalHeaderLabels(["Tissue", "Property", "Unit", "Default",
                                              *(unc.VARIANT_SHORT[v] for v in unc.VARIANTS)])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        # all rows at once: the page scrolls, not the table
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.table.itemChanged.connect(self._item_changed)
        vv.addWidget(self.table)
        # the Uncertainty analysis card's row of the overrides key: kept out of sight (this table
        # edits it); goto on the key opens the fold
        self._hidden = QWidget()
        self._hidden.hide()
        from .section import SectionPage
        self._hidden_grid = SectionPage._new_grid()
        self._hidden.setLayout(self._hidden_grid)
        v.addWidget(self._hidden)
        self.refresh()

    def attach_overrides(self, page) -> None:
        """Take the overrides row of the Uncertainty analysis card (`page`) out of sight."""
        page.move_rows((unc.OVERRIDES_KEY,), self._hidden_grid, reveal=lambda: self.b_values.setChecked(True))

    @staticmethod
    def _fold(text: str) -> tuple[QToolButton, QWidget]:
        """A closed fold: a toggle button and the widget (with a QVBoxLayout) it shows."""
        button = QToolButton(text=text)
        button.setCheckable(True)
        button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        button.setArrowType(Qt.ArrowType.RightArrow)
        button.setAutoRaise(True)
        body = QWidget()
        lay = QVBoxLayout(body)
        lay.setContentsMargins(0, 4, 0, 0)
        body.hide()

        def toggle(on: bool) -> None:
            body.setVisible(on)
            button.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        button.toggled.connect(toggle)
        return button, body

    # ── state ───────────────────────────────────────────────────────────
    def set_locked(self, locked: bool) -> None:
        if locked != self.locked:
            self.locked = locked
            self._shown_for = None
            self.refresh()

    def _toggle_lock(self) -> None:
        if self.locked:
            answer = QMessageBox.warning(
                self, "Unlock the uncertainty bounds",
                "The variants take the lower and upper bounds of the ITRUSST tissue-property table. "
                "An edited value means the variant no longer brackets these bounds. Edit only for a "
                "documented reason and report the values you used.\n\nUnlock for editing?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel)
            if answer != QMessageBox.StandardButton.Yes:
                return
        self.set_locked(not self.locked)

    def refresh(self) -> None:
        on = self.model.get("simulation.uncertainty") is True
        self.setVisible(on)
        if not on:
            return
        sig = (unc.signature(self.model, self.prestus_root), self.locked, style.active()["surface"])
        if sig == self._shown_for:
            return
        self._shown_for = sig
        self.sets = sets = unc.property_sets(self.model, self.prestus_root)
        c = style.active()
        self._fill_note(sets, c)
        self._fill_table(sets, c)
        hollow = any_inherited(sets, (*unc.BRAIN, *unc.SKULL))
        for name, tissues in (("brain", unc.BRAIN), ("skull", unc.SKULL)):
            if name in self.canvases:      # one legend, on the first plot
                draw_ranges(self.canvases[name].figure, sets, tissues, c,
                            legend=name == "brain", hollow=hollow)
                self.canvases[name].draw_idle()

    def _fill_note(self, sets: unc.PropertySets, c: dict) -> None:
        errors = [f"<b>{unc.VARIANT_LABEL[v]}</b>: <span style='color:{c['red']}'>{why}</span>"
                  for v, why in sets.errors.items()]
        self.errors.setText("<br>".join(errors))
        self.errors.setVisible(bool(errors))
        n = sets.n_edits()
        if n:
            tint = QColor(c["red"])
            tint.setAlpha(26)
            which = ", ".join(unc.VARIANT_SHORT[v] for v in sets.edits)
            self.warning.setStyleSheet(f"QLabel {{ background: {tint.name(QColor.NameFormat.HexArgb)}; "
                                       f"color: {c['ink']}; border-radius: 6px; padding: 8px 10px; }}")
            self.warning.setText(
                f"<b style='color:{c['red']}'>⚠ {n} bound value{'s' if n != 1 else ''} edited ({which}).</b> "
                "These variants no longer bracket the ITRUSST tissue-property bounds; report the values "
                "you used. <a href='show'>Show values</a>")
        self.warning.setVisible(bool(n))
        lock = ("Locked (greyed): unlock to edit a variant value. " if self.locked else
                "Edit a variant cell to change its bound; empty it to restore the ITRUSST bound. ")
        self.note.setText(lock + "Greyed values are not in the variant's bounds and use the default. "
                          "Trans-H/L run acoustic only, so their thermal properties are not used. "
                          f"Edited values are <b style='color:{c['red']}'>bold red</b>.")
        self.lock_button.setText("Unlock" if self.locked else "Lock")
        self.lock_button.setIcon(icons.icon("lock" if self.locked else "unlock"))
        self.lock_button.setToolTip("Read-only: the ITRUSST bounds. Click to edit." if self.locked
                                    else "Make the table read-only again.")
        self.reset_all.setVisible(bool(n))
        self.reset_all.setEnabled(not self.locked)
        if hasattr(self, "caption"):
            freq = (f"Attenuation is compared at the transducer frequency ({sets.freq_mhz:g} MHz)."
                    if sets.freq_from_transducer else
                    f"No transducer frequency set: attenuation is compared at {sets.freq_mhz:g} MHz.")
            self.caption.setText("ITRUSST bounds per property and the value this config uses (line). "
                                 "Temp-H/L bound brain heating, Trans-H/L brain intensity (acoustic "
                                 "properties only). " + freq)

    def _fill_table(self, sets: unc.PropertySets, c: dict) -> None:
        self._filling = True
        present = list(sets.sets.get("default", {}))
        first = [t for t in (*unc.BRAIN, *unc.SKULL) if t in present]
        tissues = first + [t for t in present if t not in first]
        rows = [(t, p) for t in tissues for p in unc.PROPERTIES]
        self.table.clearSpans()
        self.table.setRowCount(len(rows))
        grey, red = QColor(c["ink4"]), QColor(c["red"])
        fixed = Qt.ItemFlag.ItemIsEnabled
        for r, (tissue, (key, label, unit)) in enumerate(rows):
            if key in ("attenuation", "absorption"):
                unit = f"{unit} @ {sets.freq_mhz:g} MHz"
            cells = [unc.TISSUE_LABEL.get(tissue, tissue), label, unit,
                     _fmt(sets.value("default", tissue, key))]
            for i, text in enumerate(cells):
                item = QTableWidgetItem(text)
                item.setFlags(fixed)
                self.table.setItem(r, i, item)
            for i, variant in enumerate(unc.VARIANTS, start=4):
                if variant not in sets.sets:
                    item = QTableWidgetItem("—")
                    item.setFlags(fixed)
                elif variant in unc.PAIRS["trans"] and key in unc.THERMAL:
                    item = QTableWidgetItem("—")
                    item.setFlags(fixed)
                    item.setForeground(grey)
                    item.setToolTip("Trans variants run acoustic only: thermal properties are not used")
                else:
                    item = QTableWidgetItem(_fmt(sets.value(variant, tissue, key)))
                    edited = sets.edited(variant, tissue, key)
                    tips = []
                    if edited:
                        font = QFont(item.font())
                        font.setBold(True)
                        item.setFont(font)
                        item.setForeground(red)
                        tips.append(f"<b style='color:{c['red']}'>Edited</b>: ITRUSST bound "
                                    f"{_fmt(sets.bound(variant, tissue, key))}")
                    elif self.locked or sets.inherited(variant, tissue, key):
                        item.setForeground(grey)
                    if not edited and sets.inherited(variant, tissue, key):
                        tips.append("Not in the variant's bounds: the default value is used")
                    flags = fixed
                    if unc.editable(variant, key):
                        item.setData(CELL_ROLE, (variant, tissue, key))
                        if not self.locked:
                            flags |= Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEditable
                        tips.append(f"<span style='color:{FAINT}'>Key: {unc.OVERRIDES_KEY}."
                                    f"{variant}.{tissue}.{key}</span>")
                    item.setFlags(flags)
                    item.setToolTip("<br>".join(tips))
                self.table.setItem(r, i, item)
        # grouped by tissue: its name once, in a cell spanning its property rows
        n_props = len(unc.PROPERTIES)
        for k in range(len(tissues)):
            if n_props > 1:
                self.table.setSpan(k * n_props, 0, n_props, 1)
            cell = self.table.item(k * n_props, 0)
            if cell is not None:
                cell.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
                font = QFont(cell.font())
                font.setBold(True)
                cell.setFont(font)
        self._filling = False
        self._fit_table_height()

    def _fit_table_height(self) -> None:
        t = self.table
        h = t.horizontalHeader().sizeHint().height() + 2 * t.frameWidth()
        h += sum(t.rowHeight(r) for r in range(t.rowCount()))
        if t.horizontalScrollBar().isVisible():
            h += t.horizontalScrollBar().sizeHint().height()
        t.setFixedHeight(h + 2)

    # ── editing ─────────────────────────────────────────────────────────
    def _link(self, href: str) -> None:
        if href == "show":
            self.b_values.setChecked(True)

    def _item_changed(self, item: QTableWidgetItem) -> None:
        if self._filling or self.locked:
            return
        cell = item.data(CELL_ROLE)
        if not cell:
            return
        from .medium_table import _parse
        value = _parse(item.text())
        if value is None:                       # not a number: show the stored value again
            self._shown_for = None
            self.refresh()
            return
        self.set_edit(*cell, value)

    def set_edit(self, variant: str, tissue: str, prop: str, value) -> None:
        """Set (or, with [] or the shipped bound, remove) one edit of the overrides key. An
        absorption field replaces an edit of the other one (uncertainty_medium_layer.m clears the
        bound's other field)."""
        ov = copy.deepcopy(unc.overrides(self.model))
        props = ov.setdefault(variant, {}).setdefault(tissue, {})
        bound = self.sets.bound(variant, tissue, prop) if self.sets else None
        if value == [] or (bound is not None and value == bound):
            props.pop(prop, None)
        else:
            props[prop] = value
            if prop in unc.ABSORPTION:
                props.pop(next(k for k in unc.ABSORPTION if k != prop), None)
        if not props:
            ov[variant].pop(tissue)
        if not ov[variant]:
            ov.pop(variant)
        self._store(ov)

    def _reset_all(self) -> None:
        if not self.locked:
            self._store({})

    def _store(self, ov: dict) -> None:
        if ov:
            self.model.set(unc.OVERRIDES_KEY, ov)
        else:
            self.model.reset(unc.OVERRIDES_KEY)
        self._shown_for = None
        self.refresh()
        self.changed.emit()
