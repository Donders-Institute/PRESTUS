"""Colour and typography constants mirroring the MATLAB sty() palette."""

# Hex colours
ACCENT      = "#0b7b99"   # teal
ACCENT_LT   = "#d9eef2"   # teal tint
ACCENT_DK   = "#075266"   # teal dark
MANDATORY   = "#a91d25"   # deep red
SUCCESS     = "#2da13d"   # green
BG_DARK     = "#1f1f20"   # near-black header
BG_PANEL    = "#ffffff"
BG_FIG      = "#f1f1f2"
BG_EDIT     = "#f6f6f8"
TEXT        = "#151516"
TEXT_SUB    = "#666670"
BORDER      = "#d1d1d4"

FONT        = "Helvetica Neue, Segoe UI, Arial, sans-serif"
FS          = 13
FS_SM       = 11
FS_LG       = 15
FS_XL       = 20

# Qt stylesheet fragments
QSS_BASE = f"""
QWidget {{
    font-family: "Helvetica Neue", "Segoe UI", Arial, sans-serif;
    font-size: {FS}px;
    color: {TEXT};
}}
QTabWidget::pane {{
    border: none;
    background: {BG_PANEL};
}}
QTabBar::tab {{
    background: {BG_FIG};
    color: {TEXT_SUB};
    padding: 8px 16px;
    border: none;
    border-right: 2px solid {BORDER};
    min-width: 140px;
    text-align: left;
}}
QTabBar::tab:selected {{
    background: {BG_PANEL};
    color: {ACCENT};
    border-right: 3px solid {ACCENT};
    font-weight: bold;
}}
QTabBar::tab:hover {{
    background: {ACCENT_LT};
    color: {ACCENT_DK};
}}
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background: {BG_EDIT};
    border: 1px solid {BORDER};
    border-radius: 4px;
    padding: 3px 6px;
}}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{
    border: 1px solid {ACCENT};
}}
QPushButton {{
    background: {ACCENT};
    color: white;
    border: none;
    border-radius: 4px;
    padding: 5px 12px;
}}
QPushButton:hover {{
    background: {ACCENT_DK};
}}
QPushButton:pressed {{
    background: {ACCENT_DK};
}}
QPushButton.browse {{
    background: #dde3f0;
    color: {TEXT};
    padding: 4px 8px;
    min-width: 36px;
}}
QPushButton.browse:hover {{
    background: {ACCENT_LT};
}}
QPushButton.success {{
    background: {SUCCESS};
    font-weight: bold;
    font-size: {FS_LG}px;
    padding: 10px 16px;
}}
QPushButton.success:hover {{
    background: #228c31;
}}
QPushButton.danger {{
    background: {MANDATORY};
}}
QPushButton.secondary {{
    background: #dde3f0;
    color: {TEXT};
}}
QPushButton.secondary:hover {{
    background: {ACCENT_LT};
    color: {ACCENT_DK};
}}
QScrollArea {{
    border: none;
    background: {BG_PANEL};
}}
QLabel.section {{
    background: {ACCENT_LT};
    color: {ACCENT_DK};
    font-weight: bold;
    font-size: {FS_SM}px;
    padding: 3px 8px;
}}
QLabel.note {{
    color: {TEXT_SUB};
    font-size: {FS_SM}px;
}}
QLabel.mandatory-star {{
    color: {MANDATORY};
    font-weight: bold;
}}
QPlainTextEdit.log {{
    background: #191f18;
    color: #ccf0c0;
    font-family: "Courier New", monospace;
    font-size: 12px;
    border: none;
}}
QTableWidget {{
    gridline-color: {BORDER};
    background: {BG_EDIT};
}}
QTableWidget::item:selected {{
    background: {ACCENT_LT};
    color: {TEXT};
}}
QHeaderView::section {{
    background: {BG_FIG};
    border: 1px solid {BORDER};
    padding: 4px;
    font-size: {FS_SM}px;
}}
QListWidget {{
    background: {BG_EDIT};
    border: 1px solid {BORDER};
    border-radius: 4px;
}}
QListWidget::item:selected {{
    background: {ACCENT};
    color: white;
}}
QTextEdit.yaml-viewer {{
    font-family: "Courier New", monospace;
    font-size: 12px;
    background: {BG_EDIT};
    border: 1px solid {BORDER};
}}
QSplitter::handle {{
    background: {BORDER};
}}
"""
