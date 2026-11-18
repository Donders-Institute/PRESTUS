"""The Nexus tab bar.

- Bold labels for chosen tabs (e.g. the sub-tab of the placement method in use). QTabBar has
  no per-tab font, so the labels are painted here; tab shapes are the style's own.
- Status colour (set_tab's `color`): the label (or the compact icon) is drawn in it, e.g. red
  for a tab with errors; the style sheet fixes the text colour, so these labels are painted here.
- Compact mode (the main tabs, set_compactable): when the window is too narrow for every full
  label, the tabs show their icon (and a short text such as an issue count) instead, and the
  full label moves to the tooltip. Labels are never shortened with '…'. A tab set without a
  short icon (QuickStart) keeps its label.
"""

from __future__ import annotations

from PyQt6.QtCore import QRect, QSize, Qt
from PyQt6.QtGui import QColor, QFont, QFontMetrics, QIcon
from PyQt6.QtWidgets import QStyle, QStyleOptionTab, QStylePainter, QTabBar, QTabWidget


class EmphasisTabBar(QTabBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._bold: set[int] = set()
        # compact mode (main tabs only, set_compactable): full and short (text, icon, tooltip)
        self._compactable = False
        self._compact = False
        self._full: dict[int, tuple[str, QIcon, str]] = {}
        self._short: dict[int, tuple[str, QIcon]] = {}
        self._color: dict[int, QColor] = {}

    # ── bold ────────────────────────────────────────────────────────────
    def bold_tabs(self) -> set[int]:
        return set(self._bold)

    def set_bold(self, indices) -> None:
        indices = {i for i in indices if 0 <= i < self.count()}
        if indices == self._bold:
            return
        self._bold = indices
        for i in range(self.count()):           # tab widths follow the text: re-lay out
            self.setTabText(i, self.tabText(i))
        self.update()

    def _bold_font(self):
        f = self.font()
        f.setBold(True)
        return f

    def _label_font(self):
        """The weight the style sheet gives tab labels (QTabBar::tab font-weight: 600). A
        coloured label is painted here with the widget font, which does not carry it."""
        f = self.font()
        f.setWeight(QFont.Weight.DemiBold)
        return f

    # ── compact mode ────────────────────────────────────────────────────
    # The texts and icons are set on the tabs themselves (the native style draws them), so
    # owners set them through set_tab, never with setTabText / setTabIcon.
    def set_compactable(self, on: bool = True) -> None:
        self._compactable = on

    def set_tab(self, index: int, text: str, icon: QIcon | None = None, tooltip: str = "",
                short_icon: QIcon | None = None, short_text: str = "",
                color: QColor | str | None = None) -> None:
        """Tab `index`: its label, icon and tooltip, what it shows in compact mode (`short_icon`,
        `short_text`) and the colour its label is drawn in (None: the style's). Without a
        short_icon the tab keeps its label (QuickStart)."""
        self._full[index] = (text, icon or QIcon(), tooltip)
        if color:
            self._color[index] = QColor(color)
        else:
            self._color.pop(index, None)
        if short_icon is not None:
            self._short[index] = (short_text, short_icon)
        else:
            self._short.pop(index, None)
        self._apply(index)
        self._update_mode()

    def is_compact(self) -> bool:
        return self._compact

    def label(self, index: int) -> str:
        """The tab's full label, also while compact mode shows only its icon."""
        return self._full[index][0] if index in self._full else self.tabText(index)

    def label_color(self, index: int) -> QColor | None:
        return self._color.get(index)

    def _apply(self, index: int) -> None:
        if index not in self._full:
            return
        text, icon, tip = self._full[index]
        if self._compact and index in self._short:
            short_text, short_icon = self._short[index]
            self.setTabText(index, short_text)
            self.setTabIcon(index, short_icon)
            self.setTabToolTip(index, f"<b>{text.strip()}</b>" + (f"<br>{tip}" if tip else ""))
        else:
            self.setTabText(index, text)
            self.setTabIcon(index, icon)
            self.setTabToolTip(index, tip)

    def _width_as(self, index: int, text: str, icon: QIcon) -> int:
        """The tab's width with `text` and `icon`, sized by the style as QTabBar::tabSizeHint."""
        opt = QStyleOptionTab()
        self.initStyleOption(opt, index)
        opt.text, opt.icon = text, icon
        style = self.style()
        fm = QFontMetrics(self._bold_font() if index in self._bold else self._label_font())
        hspace = style.pixelMetric(QStyle.PixelMetric.PM_TabBarTabHSpace, opt, self)
        vspace = style.pixelMetric(QStyle.PixelMetric.PM_TabBarTabVSpace, opt, self)
        icon_w = 0 if icon.isNull() else self.iconSize().width() + (4 if text else 0)
        contents = QSize(fm.size(Qt.TextFlag.TextShowMnemonic, text).width() + hspace + icon_w,
                         max(fm.height(), self.iconSize().height()) + vspace)
        return style.sizeFromContents(QStyle.ContentsType.CT_TabBarTab, opt, contents, self).width()

    def _total(self, compact: bool) -> int:
        total = 0
        for i in range(self.count()):
            if not self.isTabVisible(i):
                continue
            if i not in self._full:
                total += super().tabSizeHint(i).width()
            elif compact and i in self._short:
                total += self._width_as(i, *self._short[i])
            else:
                total += self._width_as(i, *self._full[i][:2])
        return total

    def full_labels_width(self) -> int:
        return self._total(False)

    def _update_mode(self) -> None:
        if not self._compactable:
            return
        compact = self._total(False) > self.width()
        if compact != self._compact:
            self._compact = compact
            for i in range(self.count()):
                self._apply(i)
            self.updateGeometry()

    def resizeEvent(self, event) -> None:              # noqa: N802
        super().resizeEvent(event)
        self._update_mode()

    def minimumSizeHint(self) -> QSize:                # noqa: N802
        """Compactable: narrow enough for the icons, so the window can shrink below the labels."""
        size = super().minimumSizeHint()
        if self._compactable:
            size.setWidth(self._total(True))
        return size

    def tabSizeHint(self, index: int) -> QSize:       # noqa: N802
        size = super().tabSizeHint(index)
        if index in self._bold:
            text = self.tabText(index)
            extra = QFontMetrics(self._bold_font()).horizontalAdvance(text) - \
                QFontMetrics(self.font()).horizontalAdvance(text)
            size.setWidth(size.width() + max(0, extra))
        return size

    def label_rect(self, index: int) -> QRect:
        """Where tab `index` draws its icon and label: the tab's natural width centred in the
        tab, so a tab stretched to share the bar (setExpanding) keeps its icon beside the
        label rather than at the tab's left edge."""
        rect = self.tabRect(index)
        if self.tabIcon(index).isNull():
            return rect
        width = min(rect.width(), self._width_as(index, self.tabText(index), self.tabIcon(index)))
        return QRect(rect.center().x() - width // 2 + 1, rect.y(), width, rect.height())

    def separators(self) -> list[int]:
        """x of the divider between each pair of neighbouring visible tabs."""
        shown = [i for i in range(self.count()) if self.isTabVisible(i)]
        return [self.tabRect(i).right() + 1 for i in shown[:-1]]

    def paintEvent(self, event) -> None:               # noqa: N802
        painter = QStylePainter(self)
        order = [i for i in range(self.count()) if i != self.currentIndex()] + \
            ([self.currentIndex()] if self.currentIndex() >= 0 else [])     # current on top
        for i in order:
            if not self.isTabVisible(i):
                continue
            opt = QStyleOptionTab()
            self.initStyleOption(opt, i)
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabShape, opt)
            opt.rect = self.label_rect(i)
            painter.setFont(self._bold_font() if i in self._bold else self.font())
            colour = self._color.get(i)
            if colour is None or not opt.text:
                painter.drawControl(QStyle.ControlElement.CE_TabBarTabLabel, opt)
                continue
            text, opt.text = opt.text, ""
            painter.drawControl(QStyle.ControlElement.CE_TabBarTabLabel, opt)     # icon only
            rect = self.style().subElementRect(QStyle.SubElement.SE_TabBarTabText, opt, self)
            painter.setFont(self._bold_font() if i in self._bold else self._label_font())
            painter.setPen(colour)
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextShowMnemonic, text)
        # light dividers between the tabs, a third of the bar's height short at each end
        from . import style
        painter.setPen(QColor(style.active()["line2"]))
        inset = self.height() // 3
        for x in self.separators():
            painter.drawLine(x, inset, x, self.height() - inset)


def tab_widget() -> QTabWidget:
    """A document-mode QTabWidget with an EmphasisTabBar (set before any tab is added),
    fitted to the window width."""
    w = QTabWidget()
    bar = EmphasisTabBar(w)
    w.setTabBar(bar)
    w.setDocumentMode(True)
    fit_to_width(bar)
    return w


def fit_to_width(bar: QTabBar) -> None:
    """Tabs that share the window width and show their full labels, never shortened (…) and
    without scroll arrows; an EmphasisTabBar switches to icons when the labels do not fit."""
    bar.setExpanding(True)
    bar.setElideMode(Qt.TextElideMode.ElideNone)
    bar.setUsesScrollButtons(False)
