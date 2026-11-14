"""Main-window wiring, offscreen: menus, undo, the Transducers sub-tab, guarded saves."""

import os
import stat

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui import errors  # noqa: E402
from prestus_nexus.ui.main_window import MainWindow  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path):
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "session.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def _menu(window, title):
    for action in window.menuBar().actions():
        if action.text().replace("&", "") == title:
            return action.menu()
    raise AssertionError(f"no {title} menu")


def test_menus_and_names(window):
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert titles[:2] == ["File", "Edit"] and "Help" in titles
    tools = [a.text() for a in _menu(window, "Tools").actions()]
    for name in ("Localite", "PlanTUS", "TUS entry", "Calibration", "Usage statistics…"):
        assert name in tools
    shortcuts = {a.shortcut().toString() for a in _menu(window, "Help").actions()}
    assert "F1" in shortcuts


def test_undo_reverts_a_change(window):
    assert not window.undo_action.isEnabled()
    window.model.set("subject_id", 7)
    assert window.undo_action.isEnabled()
    window.undo_action.trigger()
    assert window.model.get("subject_id") != 7
    assert window.redo_action.isEnabled()
    window.redo_action.trigger()
    assert window.model.get("subject_id") == 7


def test_transducers_subtab_follows_count(window):
    tab = window.tab_pages["placement"]
    index = tab.subtabs.indexOf(window.transducers_page)
    tab.refresh()
    assert not tab.subtabs.isTabVisible(index)
    window.model.add_transducer()
    tab.refresh()
    assert tab.subtabs.isTabVisible(index)


@pytest.mark.skipif(os.name == "nt", reason="POSIX permissions")
def test_failed_save_keeps_file_and_window(window, tmp_path, monkeypatch):
    shown = []
    monkeypatch.setattr(errors, "show_error", lambda *a, **k: shown.append(a))
    folder = tmp_path / "ro"
    folder.mkdir()
    study = folder / "config_x.yaml"
    study.write_text("subject_id: 1\n")
    folder.chmod(stat.S_IRUSR | stat.S_IXUSR)       # no new files: the atomic write must fail
    try:
        window.model.path = study
        window.model.set("subject_id", 3)
        window.dirty = True
        assert window.save() is False
        assert shown and window.dirty
        assert study.read_text() == "subject_id: 1\n"
    finally:
        folder.chmod(stat.S_IRWXU)

