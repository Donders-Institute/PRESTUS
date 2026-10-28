"""Simulation tab card, shown while simulation.uncertainty is on: the default, conservative and
liberal medium properties side by side, and their ranges for brain and skull (one column per
property with its own scale, as in the ITRUSST tissue-property overview)."""

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
    if v is None:
        return "—"
    try:
        return f"{float(v):.4g}"
    except (TypeError, ValueError):
        return str(v)


def _colours(c: dict) -> dict:
    return {"default": c["ink"], "conservative": c["info"], "liberal": c["red"]}


def draw_ranges(fig, sets: unc.PropertySets, tissues: tuple[str, ...], c: dict) -> None:
    """One column per property, each on its own scale: a bar spans the conservative–liberal
    range, ticks mark the two variants (dotted where inherited from the default), a dot the
    default."""
    fig.clear()
    fig.set_facecolor(c["surface"])
    col = _colours(c)
    props = [p for p in unc.PROPERTIES if p[0] in unc.PLOTTED]
    axes = fig.subplots(1, len(props))
    for ax, (key, label, unit) in zip(axes, props):
        ax.set_facecolor(c["surface"])
        values = []
        for gi, tissue in enumerate(tissues):
            v = {k: sets.value(k, tissue, key) for k in ("default", *unc.VARIANTS) if k in sets.sets}
            v = {k: float(x) for k, x in v.items() if isinstance(x, (int, float))}
            values += list(v.values())
            bounds = [v[k] for k in unc.VARIANTS if k in v]
            if len(bounds) == 2:
                lo, hi = min(bounds), max(bounds)
                ax.bar(gi, hi - lo, bottom=lo, width=0.42, color=c["line2"], zorder=1)
            for k in unc.VARIANTS:
                if k in v:
                    ax.hlines(v[k], gi - 0.34, gi + 0.34, color=col[k], linewidth=2, zorder=2,
                              linestyles=":" if sets.inherited(k, tissue, key) else "-")
            if "default" in v:
                ax.plot(gi, v["default"], "o", color=col["default"], markersize=4, zorder=3)
        if values:
            lo, hi = min(values), max(values)
            pad = (hi - lo) * 0.12 or abs(hi) * 0.1 or 1
            ax.set_ylim(lo - pad, hi + pad)
        if key == "attenuation":
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
    handles = [Line2D([], [], color=col["default"], marker="o", linestyle="none", markersize=4, label="Default"),
               Line2D([], [], color=col["conservative"], linewidth=2, label="Conservative (lower bound)"),
               Line2D([], [], color=col["liberal"], linewidth=2, label="Liberal (upper bound)"),
               Line2D([], [], color=c["ink3"], linewidth=2, linestyle=":", label="not in the variant file: default used")]
    leg = fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=7, frameon=False,
                     bbox_to_anchor=(0.5, 1.0))
    for text in leg.get_texts():
        text.set_color(c["ink2"])
    fig.subplots_adjust(left=0.05, right=0.96, top=0.72, bottom=0.12, wspace=0.9)


class UncertaintyCard(QFrame):
    """Default / conservative / liberal medium properties; hidden unless simulation.uncertainty."""

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

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Tissue", "Property", "Unit", "Default",
                                              "Conservative", "Liberal"])
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
                     "properties a file does not list (greyed below) keep the default. " + freq)
        self.note.setText("<br>".join(lines))

    def _fill_table(self, sets: unc.PropertySets, c: dict) -> None:
        present = list(sets.sets.get("default", {}))
        first = [t for t in (*unc.BRAIN, *unc.SKULL) if t in present]
        tissues = first + [t for t in present if t not in first]
        rows = [(t, p) for t in tissues for p in unc.PROPERTIES]
        self.table.setRowCount(len(rows))
        grey = QColor(c["ink4"])
        for r, (tissue, (key, label, unit)) in enumerate(rows):
            if key == "attenuation":
                unit = f"{unit} @ {sets.freq_mhz:g} MHz"
            cells = [unc.TISSUE_LABEL.get(tissue, tissue), label, unit,
                     _fmt(sets.value("default", tissue, key))]
            for i, text in enumerate(cells):
                self.table.setItem(r, i, QTableWidgetItem(text))
            for i, variant in enumerate(unc.VARIANTS, start=4):
                if variant not in sets.sets:
                    item = QTableWidgetItem("—")
                else:
                    item = QTableWidgetItem(_fmt(sets.value(variant, tissue, key)))
                    if sets.inherited(variant, tissue, key):
                        item.setForeground(grey)
                        item.setToolTip("Not in the variant file: the default value is used")
                self.table.setItem(r, i, item)
