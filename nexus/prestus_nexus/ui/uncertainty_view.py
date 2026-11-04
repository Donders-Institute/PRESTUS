"""Simulation tab card, shown while simulation.uncertainty is on: the default medium properties
and the four variants (Temp-H/L bound heating, Trans-H/L bound intensity) side by side, and the
range of each pair for brain and skull (one column per property with its own scale, as in the
ITRUSST tissue-property overview). The Trans pair runs acoustic only, so thermal properties show
the Temp pair alone."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (QAbstractItemView, QFrame, QHeaderView, QLabel, QTableWidget,
                             QTableWidgetItem, QVBoxLayout)

from .. import uncertainty as unc
from . import style

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
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


def draw_ranges(fig, sets: unc.PropertySets, tissues: tuple[str, ...], c: dict) -> None:
    """One column per property, each on its own scale. Per tissue, the Temp pair (left) and the
    Trans pair (right) each draw a bar between their two values, with a marker per variant
    (hollow where inherited from the default); a tick marks the default. Thermal properties
    show the Temp pair only."""
    fig.clear()
    fig.set_facecolor(c["surface"])
    col = _colours(c)
    props = [p for p in unc.PROPERTIES if p[0] in unc.PLOTTED]
    axes = fig.subplots(1, len(props))
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
        ax.set_title(f"{label}\n{unit}" if unit else label, fontsize=8, color=c["ink2"])
        ax.set_xlim(-0.6, len(tissues) - 0.4)
        if len(tissues) > 1:
            ax.set_xticks(range(len(tissues)), [SKULL_TICK.get(t, t) for t in tissues], fontsize=7)
        else:
            ax.set_xticks([])
        ax.tick_params(labelsize=7, colors=c["ink3"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(c["line2"])
    handles = [Line2D([], [], color=col["default"], linewidth=1.2, label="Default"),
               Line2D([], [], color=col["temp"], linewidth=2.5, alpha=0.6, label="Temp-H/L (heating bounds)"),
               Line2D([], [], color=col["trans"], linewidth=2.5, alpha=0.6,
                      label="Trans-H/L (intensity bounds, acoustic only)"),
               Line2D([], [], color=c["ink3"], marker="^", linestyle="none", markersize=5, label="H variant"),
               Line2D([], [], color=c["ink3"], marker="v", linestyle="none", markersize=5, label="L variant"),
               Line2D([], [], color=c["ink3"], marker="o", markerfacecolor=c["surface"], linestyle="none",
                      markersize=5, label="hollow: not in the variant file, default used")]
    leg = fig.legend(handles=handles, loc="upper center", ncol=3, fontsize=7, frameon=False,
                     bbox_to_anchor=(0.5, 1.0))
    for text in leg.get_texts():
        text.set_color(c["ink2"])
    fig.subplots_adjust(left=0.05, right=0.96, top=0.68, bottom=0.12, wspace=0.9)


class UncertaintyCard(QFrame):
    """Default and Temp-H/L, Trans-H/L medium properties; hidden unless simulation.uncertainty."""

    def __init__(self, model, prestus_root: str | Path, parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.model = model
        self.prestus_root = Path(prestus_root)
        self._shown_for: tuple | None = None
        self.sets: unc.PropertySets | None = None

        v = QVBoxLayout(self)
        v.setContentsMargins(14, 10, 14, 12)
        title = QLabel("Uncertainty variants: medium properties")
        title.setObjectName("CardTitle")
        v.addWidget(title)
        self.note = QLabel()
        self.note.setObjectName("CardSub")
        self.note.setWordWrap(True)
        self.note.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(self.note)

        self.canvases: dict[str, object] = {}
        if HAVE_MPL:
            for name, heading in (("brain", "Brain"), ("skull", "Skull: S single layer · T trabecular · C cortical")):
                h = QLabel(heading)
                h.setObjectName("SubGroupTitle")
                v.addWidget(h)
                fig = Figure(figsize=(9, 2.6))
                canvas = FigureCanvasQTAgg(fig)
                canvas.setMinimumHeight(230)
                v.addWidget(canvas)
                self.canvases[name] = canvas

        self.table = QTableWidget(0, 4 + len(unc.VARIANTS))
        self.table.setHorizontalHeaderLabels(["Tissue", "Property", "Unit", "Default",
                                              *(unc.VARIANT_SHORT[v] for v in unc.VARIANTS)])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setMinimumHeight(320)
        v.addWidget(self.table)
        self.refresh()

    def _signature(self) -> tuple:
        files = []
        for variant in unc.VARIANTS:
            f = unc.variant_file(self.model, self.prestus_root, variant)
            try:
                files.append((str(f), f.stat().st_mtime))
            except OSError:
                files.append((str(f), None))
        base = repr(sorted((t, sorted(p.items())) for t, p in unc.default_properties(self.model).items()))
        return (tuple(files), base, unc.transducer_mhz(self.model), style.active()["surface"])

    def refresh(self) -> None:
        on = self.model.get("simulation.uncertainty") is True
        self.setVisible(on)
        if not on:
            return
        sig = self._signature()
        if sig == self._shown_for:
            return
        self._shown_for = sig
        self.sets = sets = unc.property_sets(self.model, self.prestus_root)
        c = style.active()
        self._fill_note(sets, c)
        for name, tissues in (("brain", unc.BRAIN), ("skull", unc.SKULL)):
            canvas = self.canvases.get(name)
            if canvas is not None:
                draw_ranges(canvas.figure, sets, tissues, c)
                canvas.draw_idle()
        self._fill_table(sets, c)

    def _fill_note(self, sets: unc.PropertySets, c: dict) -> None:
        lines = []
        for variant in unc.VARIANTS:
            label = unc.VARIANT_LABEL[variant]
            if variant in sets.errors:
                lines.append(f"<b>{label}</b>: <span style='color:{c['red']}'>{sets.errors[variant]}</span>")
            else:
                lines.append(f"<b>{label}</b>: {sets.files[variant]}")
        freq = (f"Attenuation is compared at the transducer frequency ({sets.freq_mhz:g} MHz)."
                if sets.freq_from_transducer else
                f"No transducer frequency set: attenuation is compared at {sets.freq_mhz:g} MHz.")
        lines.append("Each variant file is merged onto the default medium properties; tissues and "
                     "properties a file does not list (greyed below) keep the default. Trans-H/L run "
                     "acoustic only, so their thermal properties are not used. " + freq)
        self.note.setText("<br>".join(lines))

    def _fill_table(self, sets: unc.PropertySets, c: dict) -> None:
        present = list(sets.sets.get("default", {}))
        first = [t for t in (*unc.BRAIN, *unc.SKULL) if t in present]
        tissues = first + [t for t in present if t not in first]
        rows = [(t, p) for t in tissues for p in unc.PROPERTIES]
        self.table.setRowCount(len(rows))
        grey = QColor(c["ink4"])
        for r, (tissue, (key, label, unit)) in enumerate(rows):
            if key in ("attenuation", "absorption"):
                unit = f"{unit} @ {sets.freq_mhz:g} MHz"
            cells = [unc.TISSUE_LABEL.get(tissue, tissue), label, unit,
                     _fmt(sets.value("default", tissue, key))]
            for i, text in enumerate(cells):
                self.table.setItem(r, i, QTableWidgetItem(text))
            for i, variant in enumerate(unc.VARIANTS, start=4):
                if variant not in sets.sets:
                    item = QTableWidgetItem("—")
                elif variant in unc.PAIRS["trans"] and key in unc.THERMAL:
                    item = QTableWidgetItem("—")
                    item.setForeground(grey)
                    item.setToolTip("Trans variants run acoustic only: thermal properties are not used")
                else:
                    item = QTableWidgetItem(_fmt(sets.value(variant, tissue, key)))
                    if sets.inherited(variant, tissue, key):
                        item.setForeground(grey)
                        item.setToolTip("Not in the variant file: the default value is used")
                self.table.setItem(r, i, item)
