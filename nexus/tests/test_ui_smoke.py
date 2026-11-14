"""End-to-end UI flow, offscreen: new study -> edit -> save -> reopen."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from prestus_nexus import checks, yamlio  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import REVIEW, MainWindow  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_new_study_edit_and_save(app, schema, prestus_root, tmp_path, monkeypatch):
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    anat = tmp_path / "anat"
    anat.mkdir()
    (anat / "sub-002_T1w.nii.gz").write_bytes(b"")
    (anat / "m2m_sub-002").mkdir()
    session = Session(tmp_path / "session.yaml")

    # Untitled window (no start dialog) -> File > New study…
    first = MainWindow(ConfigModel(schema), session, prestus_root)
    window = first.new_study(tmp_path / "config_pilot.yaml", copy_default=True)
    study = tmp_path / "config_pilot.yaml"
    assert study.is_file() and (tmp_path / "config_default.yaml").is_file()
    assert window.model.schema.default_path == tmp_path / "config_default.yaml"
    assert str(study) in session.recent()
    for key, value in {"subject_id": 2, "path.anat": str(anat), "path.sim": str(tmp_path / "sims")}.items():
        window.model.set(key, value)
    assert not [i for i in checks.subject_issues(window.model) if i.level == "error"]

    window.goto("placement.manual.coordinate_space")
    page = window.pages["placement"]
    row = next(r for r in page.rows if r.node.path == "placement.manual.coordinate_space")
    row.editor.combo.setCurrentText("ras")                       # user picks RAS
    assert not any(r.node.path == "transducer.trans_pos" for r in window._page("transducer").rows)
    row = next(r for r in page.rows if r.node.path == "transducer.trans_pos")   # on Placement
    row.editor.edit.setText("-62.1, 12.4, 48")
    row.editor._commit()
    assert window.dirty
    window._page(REVIEW).refresh()
    assert window.save()

    data = yamlio.load_file(study)
    assert data["subject_id"] == 2
    assert data["path"]["anat"] == str(anat)
    assert data["placement"]["manual"]["coordinate_space"] == "ras"
    assert data["transducer"]["trans_pos"] == [-62.1, 12.4, 48]
    assert str(study) in Session(tmp_path / "session.yaml").recent()

    # File > Open… of the saved study gives a window on the same file
    reopened = window.open_study(study)
    assert reopened.model.get("placement.manual.coordinate_space") == "ras"
    reopened.dirty = False
    reopened.close()
