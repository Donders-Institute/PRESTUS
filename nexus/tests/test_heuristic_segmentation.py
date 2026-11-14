"""Heuristic page (choose a candidate) and the segmentation step (status, Run segmentation)."""

import pytest
from PyQt6.QtWidgets import QApplication

from prestus_nexus import launcher
from prestus_nexus import segmentation as sg
from prestus_nexus.model import ConfigModel
from prestus_nexus.outputs import dir_placement_run
from prestus_nexus.session import Session
from prestus_nexus.ui.main_window import HEURISTIC, MainWindow

TPARS = ("idx,trans_x,trans_y,trans_z,targ_x,targ_y,targ_z,dist_to_target,prop_intersect,"
         "mean_dist_skin,var_dist_skin,mean_dist_skull,var_dist_skull\n"
         "4,161,138,246,117,139,161,85.2,0.02,3.1,0.4,5.2,0.6\n"
         "7,160,137,245,117,139,161,84.9,0.31,2.9,0.5,5.0,0.7\n"
         "9,150,140,240,117,139,161,80.1,0.04,3.3,0.2,5.5,0.3\n")
CHOICE = "trans_x\t150\ntrans_y\t140\ntrans_z\t240\ntarg_x\t117\ntarg_y\t139\ntarg_z\t161\n"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def window(app, schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sim"))
    m.set("subject_id", 9)
    m.set("placement.heuristic.target_name", "M1")
    m.save(tmp_path / "config_t.yaml")
    w = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def _runs(monkeypatch):
    calls = []
    import prestus_nexus.ui.heuristic_page as hp
    monkeypatch.setattr(hp, "run_partial", lambda window, task, name, what, overrides, extra="", **kw:
                        calls.append((name, dict(overrides), extra)) or "job")
    return calls


def test_heuristic_page_lists_and_places_a_candidate(window, monkeypatch):
    calls = _runs(monkeypatch)
    page = window.pages[HEURISTIC]
    window.model.set("placement.mode", "heuristic")
    page.refresh()
    assert page.stack.is_empty()                                 # no candidates yet
    page.compute()
    name, ov, extra = calls[-1]
    assert name == "heuristic" and ov["io.overwrite_files"] == "always" and extra == launcher.PLACEMENT_ONLY
    run = dir_placement_run(window.model, "heuristic")
    run.mkdir(parents=True)
    (run / "tpars_sub-009_M1.csv").write_text(TPARS)
    (run / "sub-009_M1.txt").write_text(CHOICE)
    page.refresh()
    assert [c.idx for c in page.cands] == [4, 7, 9]
    assert page.choice == 9 and page.current == 9               # the heuristic's choice, shown first
    assert page.table.rowCount() == 3
    page.select(7)
    assert "above criterion_intersection" in page.readout.text()
    page.place()
    name, ov, _ = calls[-1]
    assert ov["io.overwrite_files"] == "never"                  # reuse the candidate list
    assert window.model.get("placement.heuristic.candidate") == 7
    page.use_heuristic()
    assert window.model.get("placement.heuristic.candidate") == []


def test_segmentation_status_and_run(window, monkeypatch, tmp_path):
    m = window.model
    m.set("path.seg", str(tmp_path / "seg"))
    assert sg.status(m).state == "missing"
    assert "the layered head model" in sg.needed_by(m)
    m2m = tmp_path / "seg" / "m2m_sub-009"
    (m2m / "toMNI").mkdir(parents=True)
    (m2m / "final_tissues.nii.gz").write_bytes(b"")
    assert sg.status(m).state == "incomplete"
    (m2m / "toMNI" / "MNI2Conform_nonl.nii.gz").write_bytes(b"")
    assert sg.status(m).state == "found"
    ov = dict(launcher.segmentation_overrides(rerun=True))
    assert ov["modules.segmentation_only"] == 1 and ov["simulation.medium"] == "layered"
    assert ov["io.overwrite_simnibs"] == 1
    script = launcher.write_segmentation_script(m.path, window.prestus_root).read_text()
    assert "'modules.segmentation_only', 1" in script and "PLACEMENT" not in script


def test_placement_banner_asks_for_segmentation(window, tmp_path, monkeypatch):
    m = window.model
    m.set("path.seg", str(tmp_path / "seg"))
    m.set("startup.simnibs_bin_path", str(tmp_path))
    m.set("placement.mode", "heuristic")
    panel = window.placement_seg_panel
    panel.refresh()
    assert not panel.isHidden() and panel.button.isEnabled()
    assert "heuristic placement" in panel.label.text()
    started = []
    import prestus_nexus.ui.segmentation_panel as sp
    monkeypatch.setattr(sp, "run_partial", lambda window, task, name, what, overrides, extra="", **kw:
                        started.append(dict(overrides)) or "job")
    panel.run()
    assert started and started[0]["modules.segmentation_only"] == 1 and started[0]["io.overwrite_simnibs"] == 0
    assert "Submitted" in panel.label.text()
    m.set("placement.mode", "manual")
    m.set("simulation.medium", "water")
    panel.refresh()
    assert panel.isHidden()                                      # nothing needs a segmentation
