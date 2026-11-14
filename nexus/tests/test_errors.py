"""ui.errors: guarded actions and the exception hook show a dialog instead of crashing."""

import os
import sys

import pytest
import yaml

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from prestus_nexus.ui import errors  # noqa: E402


@pytest.fixture
def shown(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(errors, "show_error", lambda parent, title, msg, details="": seen.append((title, msg, details)))
    monkeypatch.setenv("PRESTUS_NEXUS_CACHE", str(tmp_path))
    return seen


def test_guarded_turns_expected_errors_into_a_dialog(shown, tmp_path):
    with errors.guarded(None, "Save study") as g:
        open(tmp_path / "missing" / "x.yaml", encoding="utf-8")
    assert g.failed and isinstance(g.error, OSError)
    title, msg, details = shown[0]
    assert title == "Save study failed" and "x.yaml" in msg and "Traceback" in details
    assert "Save study" in (tmp_path / errors.LOG_NAME).read_text(encoding="utf-8")
    with errors.guarded(None, "Open study") as g:
        yaml.safe_load("a: [")
    assert g.failed and len(shown) == 2
    with errors.guarded(None, "Reload") as g:
        pass
    assert not g.failed


def test_guarded_lets_other_errors_through(shown):
    with pytest.raises(KeyError):
        with errors.guarded(None, "Save"):
            raise KeyError("x")
    assert not shown


def test_guarded_as_decorator(shown):
    @errors.guarded(None, "New study")
    def act():
        raise ValueError("bad value")
    assert act() is None and shown[0][1].endswith("ValueError: bad value")


def test_excepthook_logs_and_shows(shown, monkeypatch):
    monkeypatch.setattr(sys, "excepthook", sys.excepthook)
    import threading
    monkeypatch.setattr(threading, "excepthook", threading.excepthook)
    errors.install_excepthook()
    try:
        raise RuntimeError("boom")
    except RuntimeError:
        sys.excepthook(*sys.exc_info())
    assert shown and "RuntimeError: boom" in shown[0][1]
