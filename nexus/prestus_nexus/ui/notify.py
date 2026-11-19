"""Tell the user that a run or a cluster job ended while the Nexus is in the background:
a system notification where the platform offers one (QSystemTrayIcon), and the Dock or
taskbar alert (QApplication.alert) in any case. Nothing while the Nexus is the active window."""

from __future__ import annotations

from PyQt6.QtWidgets import QApplication, QSystemTrayIcon, QWidget

_tray: QSystemTrayIcon | None = None


def notify(window: QWidget, title: str, text: str) -> bool:
    """Notify unless window is active; returns whether it did."""
    global _tray
    if window is None or window.isActiveWindow():
        return False
    try:
        if QSystemTrayIcon.isSystemTrayAvailable() and QSystemTrayIcon.supportsMessages():
            if _tray is None:
                _tray = QSystemTrayIcon(window.windowIcon(), QApplication.instance())
                _tray.setToolTip("PRESTUS Nexus")
            _tray.show()
            _tray.showMessage(title, text, QSystemTrayIcon.MessageIcon.Information, 10_000)
        QApplication.alert(window)
    except Exception:                    # a notification must never break the Nexus
        return False
    return True
