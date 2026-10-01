"""Look of the Nexus: the palette and logo of the PRESTUS HTML reports."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QIcon, QPainter, QPainterPath, QPalette, QPixmap
from PyQt6.QtWidgets import QApplication

# Palette of the report skin (functions/report/css_styles_base.m, light theme).
# "accent" is the acoustic-content colour (acoustic = cool/blue, thermal = warm) and stays
# reserved for acoustic data (transducer_view, timing_view) — it is not the brand colour.
# "brand" (the logo's teal) is for chrome: the header, active tab, primary button, focus ring.
#
# C is also what the header bar is always styled from (see _stylesheet()'s #HeaderBar rules),
# regardless of the active theme, so the logo sits on its familiar light background in dark
# mode too. Anything else should read from the active palette (style.active()), not C.
C = {
    "bg": "#e9edf1", "bg2": "#dfe5ec", "surface": "#ffffff", "surface2": "#f5f7fa", "surface3": "#eef2f6",
    "ink": "#10161d", "ink2": "#3d4956", "ink3": "#6b7886", "ink4": "#9aa7b4",
    "line": "#dde3ea", "line2": "#cdd6df",
    "accent": "#1f6fb2", "accent_soft": "#e7f1f9", "thermal": "#d9482a", "thermal_soft": "#fceee9",
    "brand": "#097b9a", "brand_soft": "#e2f0f3",
    "green": "#1f9d57", "amber": "#d98a1f", "red": "#d63b3b", "info": "#2f7fd1",
}

# Dark theme, same roles as C, tuned to keep the same brand/green/amber/red/info hues legible
# on a dark surface.
C_DARK = {
    "bg": "#181d22", "bg2": "#12161a", "surface": "#20262c", "surface2": "#272e35", "surface3": "#2f373f",
    "ink": "#e8edf2", "ink2": "#c2cbd4", "ink3": "#8fa0ac", "ink4": "#647680",
    "line": "#39424a", "line2": "#49535c",
    "accent": "#5aa8e6", "accent_soft": "#1d3040", "thermal": "#e8815e", "thermal_soft": "#3a2820",
    "brand": "#2bb6d6", "brand_soft": "#173239",
    "green": "#3fcb7c", "amber": "#e3a63e", "red": "#e2685f", "info": "#5b9fe6",
}

_dark = False


def active() -> dict:
    """The palette currently applied to everything except the header bar."""
    return C_DARK if _dark else C


def is_dark() -> bool:
    return _dark


def logo_path(prestus_root: Path | None) -> Path | None:
    for root in [prestus_root, Path(__file__).resolve().parents[3]]:
        if root is None:
            continue
        p = Path(root) / "functions" / "report" / "assets" / "logo_PRESTUS.png"
        if p.is_file():
            return p
    return None


def logo_pixmap(prestus_root: Path | None, height: int = 34) -> QPixmap | None:
    p = logo_path(prestus_root)
    if p is None:
        return None
    pix = QPixmap(str(p))
    if pix.isNull():
        return None
    return pix.scaledToHeight(height, Qt.TransformationMode.SmoothTransformation)


def app_icon(prestus_root: Path | None) -> QIcon | None:
    """App/Dock icon: the logo centered on a white canvas with macOS-style rounded corners.

    The source PNG is a non-square opaque rectangle, so handing it to QIcon
    directly leaves the square icon slot's margins transparent (visible as
    see-through in the macOS Dock/Finder); pad it onto white instead. This
    runs as a plain script rather than a bundled .app, so macOS won't apply
    its own squircle mask — round the corners ourselves (macOS's ~22.37%
    corner-radius convention) so the icon reads as native in Dock/Finder.
    """
    p = logo_path(prestus_root)
    if p is None:
        return None
    pix = QPixmap(str(p))
    if pix.isNull():
        return None
    side = max(pix.width(), pix.height())
    square = QPixmap(side, side)
    square.fill(Qt.GlobalColor.transparent)
    painter = QPainter(square)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    radius = side * 0.2237
    path = QPainterPath()
    path.addRoundedRect(0, 0, side, side, radius, radius)
    painter.setClipPath(path)
    painter.fillRect(square.rect(), QColor(C["surface"]))
    painter.drawPixmap((side - pix.width()) // 2, (side - pix.height()) // 2, pix)
    painter.end()
    return QIcon(square)


def _stylesheet(c: dict) -> str:
    """The app stylesheet for palette `c`, except the header bar, which is always styled
    from the light palette `C` (see module docstring on `C`) so it keeps its usual look —
    logo included — no matter the active theme."""
    return f"""
QMainWindow, QDialog {{ background: {c['bg']}; }}
QWidget {{ color: {c['ink']}; font-size: 13px; }}
QToolTip {{ background: #10161d; color: #ffffff; border: none; padding: 6px 8px; border-radius: 4px; }}

#HeaderBar {{ background: {C['surface']}; border-bottom: 1px solid {C['line']}; }}
#HeaderTitle {{ font-size: 16px; font-weight: 600; color: {C['ink']}; }}
#HeaderSub {{ color: {C['ink3']}; }}
#HeaderBar QPushButton {{ background: {C['surface']}; color: {C['ink']}; border: 1px solid {C['line2']}; }}
#HeaderBar QPushButton:hover {{ border-color: {C['brand']}; }}
#HeaderBar QPushButton[primary="true"] {{ background: {C['brand']}; color: white; border-color: {C['brand']}; }}
.Badge, QLabel[badge="true"] {{ border-radius: 9px; padding: 2px 8px; font-weight: 600; }}

QTabWidget::pane {{ border: none; background: {c['bg']}; }}
QTabBar::tab {{ background: transparent; color: {c['ink3']}; padding: 9px 18px; margin-right: 2px;
               border: none; border-bottom: 3px solid transparent; font-weight: 600; }}
QTabBar::tab:selected {{ color: {c['brand']}; border-bottom: 3px solid {c['brand']}; }}
QTabBar::tab:hover {{ color: {c['ink']}; }}

#TabIntro {{ color: {c['ink2']}; }}
#Card {{ background: {c['surface']}; border: 1px solid {c['line']}; border-radius: 8px; }}
#CardTitle {{ font-size: 14px; font-weight: 600; color: {c['ink']}; }}
#CardSub {{ color: {c['ink3']}; }}
#SubGroup {{ background: {c['surface2']}; border: 1px solid {c['line']}; border-radius: 6px; }}
#SubGroupTitle {{ font-weight: 600; color: {c['ink2']}; }}
#FieldLabel {{ color: {c['ink']}; }}
#FieldKey {{ color: {c['ink4']}; font-size: 11px; }}
#AdvancedBadge {{ color: {c['ink4']}; font-size: 10px; }}

QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QPlainTextEdit, QTextEdit {{
    background: {c['surface']}; border: 1px solid {c['line2']}; border-radius: 4px; padding: 4px 6px; }}
QLineEdit:focus, QComboBox:focus, QPlainTextEdit:focus {{ border: 1px solid {c['brand']}; }}
QLineEdit:disabled {{ background: {c['surface3']}; color: {c['ink3']}; }}
QPushButton {{ background: {c['surface']}; border: 1px solid {c['line2']}; border-radius: 5px; padding: 5px 12px; }}
QPushButton:hover {{ border-color: {c['brand']}; }}
QPushButton[primary="true"] {{ background: {c['brand']}; color: white; border-color: {c['brand']}; font-weight: 600; }}
QPushButton:disabled {{ color: {c['ink4']}; }}
QToolButton {{ border: none; padding: 2px 4px; border-radius: 5px; }}
QToolButton:hover {{ background: {c['surface3']}; }}
QCheckBox {{ spacing: 6px; }}
QTableWidget, QTreeWidget, QListWidget {{ background: {c['surface']}; border: 1px solid {c['line']};
    border-radius: 6px; gridline-color: {c['line']}; }}
QHeaderView::section {{ background: {c['surface2']}; border: none; border-bottom: 1px solid {c['line']};
    padding: 4px 6px; font-weight: 600; color: {c['ink2']}; }}
QScrollArea {{ border: none; background: transparent; }}
QScrollArea > QWidget > QWidget {{ background: transparent; }}
QDockWidget {{ color: {c['ink']}; }}
QDockWidget::title {{ background: {c['surface2']}; padding: 4px 8px; }}
QMenuBar {{ background: {c['surface']}; color: {c['ink']}; }}
QMenuBar::item:selected {{ background: {c['surface3']}; }}
QMenu {{ background: {c['surface']}; color: {c['ink']}; border: 1px solid {c['line']}; }}
QMenu::item:selected {{ background: {c['brand_soft']}; color: {c['ink']}; }}
QStatusBar {{ background: {c['surface']}; border-top: 1px solid {c['line']}; color: {c['ink3']}; }}
"""


def apply(app: QApplication, prestus_root: Path | None = None, dark: bool = False) -> None:
    global _dark
    _dark = dark
    c = active()
    app.setStyle("Fusion")
    pal = app.palette()
    for role, key in [(QPalette.ColorRole.Window, "bg"), (QPalette.ColorRole.Base, "surface"),
                      (QPalette.ColorRole.AlternateBase, "surface2"), (QPalette.ColorRole.Text, "ink"),
                      (QPalette.ColorRole.WindowText, "ink"), (QPalette.ColorRole.Button, "surface"),
                      (QPalette.ColorRole.ButtonText, "ink"), (QPalette.ColorRole.Highlight, "brand"),
                      (QPalette.ColorRole.ToolTipBase, "ink"), (QPalette.ColorRole.PlaceholderText, "ink4")]:
        pal.setColor(role, QColor(c[key]))
    pal.setColor(QPalette.ColorRole.HighlightedText, QColor("#ffffff"))
    pal.setColor(QPalette.ColorRole.ToolTipText, QColor("#ffffff"))
    app.setPalette(pal)
    app.setStyleSheet(_stylesheet(c))
    icon = app_icon(prestus_root)
    if icon is not None:
        app.setWindowIcon(icon)
