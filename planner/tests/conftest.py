import sys
from pathlib import Path

import pytest

PLANNER = Path(__file__).resolve().parents[1]
PRESTUS_ROOT = PLANNER.parent
sys.path.insert(0, str(PLANNER))

from prestus_planner.schema import Schema  # noqa: E402


@pytest.fixture(scope="session")
def prestus_root() -> Path:
    return PRESTUS_ROOT


@pytest.fixture(scope="session")
def schema(prestus_root) -> Schema:
    return Schema(prestus_root / "config" / "config_default.yaml")


@pytest.fixture(autouse=True)
def _no_modal_dialogs(monkeypatch):
    """Offscreen, a modal dialog would block forever: fail instead.

    Tests that expect a dialog patch it themselves (their monkeypatch wins).
    """
    try:
        from PyQt6.QtWidgets import QMessageBox
    except ImportError:
        return

    def unexpected(*args, **kwargs):
        pytest.fail(f"unexpected dialog: {args[1:3]}")
    for name in ("information", "warning", "critical", "question"):
        monkeypatch.setattr(QMessageBox, name, unexpected)
