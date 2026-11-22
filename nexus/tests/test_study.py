"""Study manifest model, Localite import, study launcher and the Study page."""

import os
import stat
import subprocess

import pytest

from prestus_nexus import launcher
from prestus_nexus.model import ConfigModel
from prestus_nexus.study import (Manifest, from_neuronav_summary, parse_subject_ids, resolve_manifest,
                                 sanitize_affix, study_report_path, subject_list_report_path)


def test_manifest_roundtrip_tsv_and_csv(tmp_path):
    m = Manifest()
    m.add_row({"subject_id": "1", "condition": "M1", "analysis.target_mni_mm": "[-38, -20, 60]"})
    m.add_row({"subject_id": "2", "condition": "M1"})
    for name in ("s.tsv", "s.csv"):
        p = m.save(tmp_path / name)
        back = Manifest.load(p)
        assert back.columns == ["subject_id", "condition", "analysis.target_mni_mm"]
        assert back.rows[0]["analysis.target_mni_mm"] == "[-38, -20, 60]"
        assert back.rows[1]["analysis.target_mni_mm"] == ""
    assert "# " not in (tmp_path / "s.tsv").read_text()


def test_parse_subject_ids():
    assert parse_subject_ids("1-3, 7 sub-020") == ["1", "2", "3", "7", "sub-020"]
    assert parse_subject_ids("008-010;3 3") == ["008", "009", "010", "3"]
    assert parse_subject_ids("  ") == []
    for bad in ("5-2", "abc", "1-x"):
        with pytest.raises(ValueError):
            parse_subject_ids(bad)


def test_manifest_comments_and_key_columns(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("# comment\nsubject_id\tcondition\tnotes\tplacement.localite.target\n1\tA\tx\tL_MD\n")
    m = Manifest.load(p)
    assert m.key_columns() == ["placement.localite.target"]
    assert m.subjects() == [1] and m.conditions() == ["A"]


def test_affixes_match_matlab_rule():
    m = Manifest()
    for r in ({"subject_id": "1", "condition": "M1"}, {"subject_id": "1", "condition": "PUL"},
              {"subject_id": "2", "io.output_affix": "custom"}, {"subject_id": "3"}):
        m.add_row(r)
    assert m.affixes("_desc-temphigh") == ["_desc-temphigh-M1", "_desc-temphigh-PUL", "_desc-custom", "_desc-temphigh"]
    assert m.affixes("") == ["_desc-M1", "_desc-PUL", "_desc-custom", ""]
    assert sanitize_affix("cond A") == "_desc-cond_A"


def test_issues(schema):
    m = Manifest(columns=["subject_id", "condition", "include", "placement.localite.targte", "bad key"])
    m.add_row({"subject_id": "1", "condition": "A"})
    m.add_row({"subject_id": "sub-001", "condition": "A"})         # same outputs as row 1
    m.add_row({"subject_id": "x"})
    m.add_row({"subject_id": "4", "include": "maybe"})
    m.add_row({"subject_id": "1", "condition": "A", "include": "0"})  # excluded: no clash
    by = {(i.row, i.column): i for i in m.issues(schema)}
    assert by[(None, "bad key")].level == "error"
    assert "did you mean placement.localite.target" in by[(None, "placement.localite.targte")].message
    assert by[(1, "condition")].level == "error"
    assert by[(2, "subject_id")].level == "error"
    assert by[(3, "include")].level == "error"
    assert (4, "condition") not in by


def test_import_neuronav_summary(tmp_path):
    p = tmp_path / "neuronav_ingest_summary.csv"
    p.write_text("sub,ses,target,series_index,status\n"
                 "sub-001,ses-01,L_MD,1,2 target(s) written\n"
                 "sub-001,ses-01,R_MD,2,2 target(s) written\n"
                 "sub-002,ses-01,L_MD,1,1 target(s) written\n"
                 "sub-002,ses-02,L_MD,1,1 target(s) written\n"
                 "sub-003,ses-01,,,no usable markers\n")
    m = from_neuronav_summary(p)
    assert [r["subject_id"] for r in m.rows] == ["1", "1", "2", "2"]
    assert [r["condition"] for r in m.rows] == ["L_MD", "R_MD", "L_MD-ses-01", "L_MD-ses-02"]
    assert all(r["placement.mode"] == "localite" for r in m.rows)
    assert m.rows[3]["placement.localite.session"] == "ses-02"
    assert not m.issues()


def test_resolve_manifest(tmp_path):
    cfg = tmp_path / "cfg" / "config_s.yaml"
    assert resolve_manifest(cfg, "study.tsv") == tmp_path / "cfg" / "study.tsv"
    assert resolve_manifest(cfg, str(tmp_path / "abs" / "m.tsv")) == tmp_path / "abs" / "m.tsv"
    assert resolve_manifest(cfg, "") is None


def test_study_launcher_and_commands(tmp_path, prestus_root):
    study = tmp_path / "config_grp.yaml"
    path = launcher.write_study_launcher(study, prestus_root)
    assert path.name == "run_study_grp.m"
    text = path.read_text()
    assert f"prestus_study_start('{study.resolve()}', 'Platform', 'matlab');" in text
    argv = launcher.study_dry_run_command("/x/matlab", study, prestus_root)
    assert argv[:2] == ["/x/matlab", "-batch"] and "'DryRun', true" in argv[2]
    assert launcher.OK_SENTINEL in argv[2]
    argv = launcher.group_report_command("/x/matlab", study, prestus_root)
    assert "prestus_group_report_start(" in argv[2]


# ── Study page (offscreen) ──────────────────────────────────────────────

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.mark.posix
def test_study_page_edit_save_and_dry_run(schema, prestus_root, tmp_path):
    pytest.importorskip("PyQt6")
    from PyQt6.QtCore import QEventLoop, QTimer
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, STUDY, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sim"))
    m.save(tmp_path / "config_g.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window.show_special(STUDY)
    page.refresh()

    page.path.setText("study.tsv")
    page._path_edited()
    assert window.model.get("study.manifest") == "study.tsv"
    page.add_column("transducer.target_isppa_wcm2")
    page.add_row({"subject_id": "1", "condition": "A", "transducer.target_isppa_wcm2": "30"})
    page.add_row({"subject_id": "2", "condition": "A"})
    page.table.item(1, 2).setText("20")                  # edit a cell in the table
    assert page.manifest.rows[1]["transducer.target_isppa_wcm2"] == "20"
    saved = page.save_manifest()
    assert saved == tmp_path / "study.tsv"
    assert Manifest.load(saved).rows[1]["transducer.target_isppa_wcm2"] == "20"
    assert "2 rows, 2 included" in page.summary.text()

    # dry run through the Run page with a fake MATLAB that records its arguments
    fake = tmp_path / "matlab"
    fake.write_text(f"#!/bin/sh\necho \"$@\" > {tmp_path / 'args.txt'}\necho '{launcher.OK_SENTINEL}'\nexit 0\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    run = window._page(RUN)
    run.matlab.setText(str(fake))
    results = []
    loop = QEventLoop()
    run.finished.connect(lambda task, ok: (results.append((task, ok)), loop.quit()))
    QTimer.singleShot(20000, loop.quit)
    page.dry_run()
    loop.exec()
    assert results == [("Study dry run", True)]
    assert "prestus_study_start" in (tmp_path / "args.txt").read_text()


def test_study_page_empty_state_add_subjects_default_save(schema, prestus_root, tmp_path):
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import STUDY, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m = ConfigModel(schema)
    m.save(tmp_path / "config_e.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window.show_special(STUDY)
    page.refresh()
    assert page.stack.is_empty() and page.summary.text() == ""
    page.add_subjects("1-3")
    assert not page.stack.is_empty()
    assert [r["subject_id"] for r in page.manifest.rows] == ["1", "2", "3"]
    assert "subject_id" in page.table.horizontalHeaderItem(0).toolTip()
    # no table file yet: saved as <config>_runs.tsv next to the config, rows kept
    saved = page.save_manifest()
    assert saved == tmp_path / "config_e_runs.tsv"
    assert window.model.get("study.manifest") == "config_e_runs.tsv"
    assert len(Manifest.load(saved).rows) == 3
    # pointing at a new file keeps the rows instead of clearing the table
    page.path.setText("other.tsv")
    page._path_edited()
    assert len(page.manifest.rows) == 3


MATLAB = launcher.find_matlab()


@pytest.mark.matlab
@pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")
def test_study_dry_run_in_matlab(schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sim"))
    m.set("study.manifest", "study.tsv")
    study = m.save(tmp_path / "config_d.yaml")
    man = Manifest(columns=["subject_id", "condition", "transducer.target_isppa_wcm2"])
    man.add_row({"subject_id": "1", "condition": "A", "transducer.target_isppa_wcm2": "30"})
    man.add_row({"subject_id": "2", "condition": "B"})
    man.save(tmp_path / "study.tsv")
    proc = subprocess.run(launcher.study_dry_run_command(MATLAB, study, prestus_root),
                          capture_output=True, text=True, timeout=600, cwd=tmp_path)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0 and launcher.OK_SENTINEL in out, out[-3000:]
    assert "Dry run: 2 run(s) ready, 0 with config errors" in out


@pytest.mark.parametrize("name,manifest,expected", [
    ("", "runs.tsv", "group/runs/group_runs_report.html"),           # manifest file name
    ("pilot study", "runs.tsv", "group/pilot_study/group_pilot_study_report.html"),
    (None, None, None),                                               # no name at all
])
def test_study_report_path(tmp_path, name, manifest, expected):
    # study_setup.m: <path.sim>/group/<name>/, generate_group_report.m: group_<name>_report.html
    mp = tmp_path / manifest if manifest else None
    got = study_report_path(str(tmp_path / "sims"), name, mp)
    assert got == (tmp_path / "sims" / expected if expected else None)
    assert study_report_path(None, "x", mp) is None


def test_study_page_open_group_report(schema, prestus_root, tmp_path, monkeypatch):
    pytest.importorskip("PyQt6")
    from PyQt6.QtGui import QDesktopServices
    from PyQt6.QtWidgets import QApplication, QMessageBox

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import STUDY, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sims"))
    m.set("study.manifest", "runs.tsv")
    m.save(tmp_path / "config_g.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window.show_special(STUDY)
    opened, told = [], []
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: told.append(a[2]))
    page.open_group_report()                                 # not built yet: says so
    assert not opened and told and "No group report yet" in told[0]
    report = tmp_path / "sims" / "group" / "runs" / "group_runs_report.html"
    report.parent.mkdir(parents=True)
    report.write_text("<html></html>")
    page.open_group_report()
    assert opened == [report.as_posix()]


@pytest.mark.parametrize("medium,affix,expected", [
    ("layered", None, "group_layered_report.html"),
    ("water", "_desc-pilot", "group_water_report_desc-pilot.html"),
    (None, None, None),                                               # no medium
])
def test_subject_list_report_path(tmp_path, medium, affix, expected):
    # generate_group_report.m without a study: <path.sim>/group_<medium>_report<affix>.html
    got = subject_list_report_path(str(tmp_path), medium, affix)
    assert got == (tmp_path / expected if expected else None)
    assert subject_list_report_path(None, "layered", None) is None


@pytest.mark.parametrize("study_config", [False, True])
def test_study_page_group_report_without_a_study(schema, prestus_root, tmp_path, monkeypatch, study_config):
    # a config without study.name / manifest opens the report over the subjects in path.sim;
    # a study config never falls back to it (negative control)
    pytest.importorskip("PyQt6")
    from PyQt6.QtGui import QDesktopServices
    from PyQt6.QtWidgets import QApplication, QMessageBox

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import STUDY, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sims"))
    m.set("simulation.medium", "layered")
    if study_config:
        m.set("study.name", "pilot")
    m.save(tmp_path / "config_g.yaml")
    report = tmp_path / "sims" / "group_layered_report.html"
    report.parent.mkdir(parents=True)
    report.write_text("<html></html>")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window.show_special(STUDY)
    opened, told = [], []
    monkeypatch.setattr(QDesktopServices, "openUrl", lambda url: opened.append(url.toLocalFile()))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: told.append(a[2]))
    page.open_group_report()
    if study_config:
        assert not opened and "No group report yet" in told[0]
    else:
        assert opened == [report.as_posix()] and not told
    window.dirty = False
    window.close()
