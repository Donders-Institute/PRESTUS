"""results.py finds the runs in a subject's output folder (by output affix) and the Results
tab shows each one: here a folder as run_multitransducer_modes leaves it (async variants
with per-transducer sub-runs, coherent, sequential steps) plus a plain run."""

import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

from prestus_nexus import results as rs
from prestus_nexus.model import ConfigModel

TAG, MED = "sub-003", "layered"
RUNS = {  # affix: (writes a CSV, NIfTI maps)
    "": (True, ["Itp"]),
    "_desc-coherent": (True, ["Itp", "heating", "heatrise"]),
    "_desc-async_trans-01": (True, ["Itp"]),
    "_desc-async_trans-02": (True, ["intensity"]),           # older name of the ISPPA map
    "_desc-async_desc-asyncJoint": (True, ["Itp", "heatrise", "CEM43_end"]),
    "_desc-async_desc-asyncPreproc": (False, ["medium_masks"]),   # not a result
    "_desc-seq-t1": (True, ["Itp", "heating"]),
    "_desc-seq-t2": (False, ["heating"]),
}


def _folder(tmp_path, nifti=lambda p: p.write_bytes(b"")):
    base = tmp_path / "sims" / TAG
    (base / "nii").mkdir(parents=True)
    for affix, (has_csv, maps) in RUNS.items():
        if has_csv:
            (base / f"{TAG}_{MED}{affix}.csv").write_text(
                "subject_id,Isppa_brain,riseT_brain\n3,1.0,0.1\n3,%s,0.2\n" % (len(affix) + 1))
            (base / f"{TAG}_{MED}_report{affix}.html").write_text("<html></html>")
        for m in maps:
            nifti(base / "nii" / f"{TAG}_{MED}_T1w{affix}_{m}.nii.gz")
    (base / f"{TAG}_{MED}_desc-coherent_targeting.csv").write_text("a\n1\n")
    (base / f"{TAG}_{MED}_sequential_report.html").write_text("<html></html>")
    (base / f"{TAG}_placement_desc-async_desc-asyncPreproc.json").write_text(json.dumps(
        {"mode": "manual", "transducers": [{"trans_pos": [1, 2, 3], "focus_pos": [4, 5, 6]},
                                           {"trans_pos": [7, 8, 9], "focus_pos": [10, 11, 12]}]}))
    return base


def _model(schema, tmp_path, coupling="async") -> ConfigModel:
    m = ConfigModel(schema)
    m.set("subject_id", 3)
    m.set("path.sim", str(tmp_path / "sims"))
    m.set("simulation.medium", MED)
    m.set("simulation.transducer_coupling", coupling)
    m.set("io.output_affix", "async")
    m.add_transducer()
    return m


def test_discover_runs(tmp_path):
    runs = rs.discover_runs(_folder(tmp_path), TAG)
    by_affix = {r.affix: r for r in runs}
    assert set(by_affix) == {a for a in RUNS if a != "_desc-async_desc-asyncPreproc"}
    assert by_affix["_desc-async_trans-02"].maps["Intensity (ISPPA)"].name.endswith("_intensity.nii.gz")
    assert by_affix["_desc-seq-t2"].csv is None and "Max temperature" in by_affix["_desc-seq-t2"].maps
    assert by_affix["_desc-coherent"].report.name == f"{TAG}_{MED}_report_desc-coherent.html"


def test_describe_affix():
    assert rs.describe_affix("_desc-async_trans-02") == ("async · transducer 2 alone", "transducer", 2)
    assert rs.describe_affix("_desc-async-pfov_desc-asyncJoint")[1:] == ("joint", None)
    assert "sweep 3" in rs.describe_affix("_desc-asyncJointS03")[0]
    assert rs.describe_affix("")[0] == "(no affix)"


def test_default_run_is_async_joint(schema, tmp_path):
    _folder(tmp_path)
    m = _model(schema, tmp_path)
    runs = rs.runs_for(m)
    assert runs[rs.default_run(runs, m)].affix == "_desc-async_desc-asyncJoint"
    m.set("simulation.transducer_coupling", "coherent")
    m.set("io.output_affix", "coherent")
    assert runs[rs.default_run(runs, m)].affix == "_desc-coherent"


def test_related_runs_and_positions(tmp_path):
    runs = rs.discover_runs(_folder(tmp_path), TAG)
    joint = next(r for r in runs if r.kind == "joint")
    parts = rs.related_runs(runs, joint)
    assert [r.transducer for r in parts] == [1, 2]
    assert rs.placement_positions(joint) == [([1, 2, 3], [4, 5, 6]), ([7, 8, 9], [10, 11, 12])]
    assert rs.placement_positions(parts[1]) == [([7, 8, 9], [10, 11, 12])]
    assert rs.placement_positions(next(r for r in runs if r.affix == "_desc-seq-t1")) == []


def test_read_metrics_last_row(tmp_path):
    runs = rs.discover_runs(_folder(tmp_path), TAG)
    run = next(r for r in runs if r.affix == "_desc-coherent")
    assert rs.read_metrics(run.csv)["riseT_brain"] == "0.2"


def test_results_page_lists_runs(schema, tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    nib = pytest.importorskip("nibabel")
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.results_page import ResultsPage
    app = QApplication.instance() or QApplication([])  # noqa: F841  keep it alive

    def nifti(path):
        vol = np.zeros((8, 8, 8), np.float32)
        vol[2, 3, 4] = 5
        nib.save(nib.Nifti1Image(vol, np.eye(4)), str(path))

    _folder(tmp_path, nifti)
    anat = tmp_path / "anat"
    anat.mkdir()
    nib.save(nib.Nifti1Image(np.ones((8, 8, 8), np.float32), np.eye(4)), str(anat / "sub-003_T1w.nii.gz"))
    m = _model(schema, tmp_path)
    m.set("path.anat", str(anat))
    page = ResultsPage(SimpleNamespace(model=m))
    page.refresh()
    assert page.subtabs.currentWidget() is page.overview      # the reporting overview opens first
    assert page.overview.report is not None and page.acoustic.viewer.anat is None
    page.show_view(page.acoustic)                             # maps: under Spatial
    assert page.subtabs.currentWidget() is page.spatial
    page.acoustic.viewer.wait()                               # volumes load in the background

    labels = [page.run_combo.itemText(i) for i in range(page.run_combo.count())]
    assert "async · both transducers (async)" in labels
    assert not any("alone" in lbl for lbl in labels)          # sub-runs are parts, not runs
    assert page.run_combo.currentText() == "async · both transducers (async)"
    assert page.parts.count() == 3                            # joint, transducer 1, transducer 2
    assert len(page.acoustic.viewer.markers) == 4             # 2 transducers + 2 foci
    assert np.allclose(page.acoustic.viewer.cursor, [2, 3, 4])   # centred on the map's peak
    viewer = page.acoustic.viewer                             # ISPPA: 0..peak, in W/cm², read out
    viewer._draw()
    assert viewer._over_lim == (0.0, 5.0) and viewer.threshold.value() == 10
    assert "≥ 0.5 W/cm²" in viewer.threshold_label.text()
    assert "overlay 5 W/cm²" in viewer.readout.text()
    assert viewer._colorbar is not None and viewer._colorbar.ax.get_visible()
    page.parts.setCurrentIndex(2)
    page.acoustic.viewer.wait()
    assert page.acoustic.viewer.overlay.path.name.endswith("_trans-02_intensity.nii.gz")
    assert len(page.acoustic.viewer.markers) == 2
    assert page.compare.rowCount() == len(labels)
    assert page.report_combo.findText(f"{TAG}_{MED}_sequential_report.html") >= 0

    # Sub-tabs other than the shown one load when opened.
    assert page.thermal.viewer.overlay is None and page.thermal.viewer.anat is None
    page.show_view(page.thermal)
    page.thermal.viewer.wait()
    assert page.thermal.viewer.overlay is not None

    # Showing the tab again with nothing changed reloads nothing; a changed file does.
    overlay = page.thermal.viewer.overlay
    page.refresh()
    assert page.thermal.combo.currentText() and page.thermal.viewer.overlay is overlay
    csv_path = tmp_path / "sims" / TAG / f"{TAG}_{MED}_desc-coherent.csv"
    os.utime(csv_path, ns=(csv_path.stat().st_atime_ns, csv_path.stat().st_mtime_ns + 10**9))
    page.compare.setRowCount(0)
    page.refresh()
    assert page.compare.rowCount() == len(labels)
    assert page.run_combo.currentText() == "async · both transducers (async)"   # kept


def _heat_file(path, n=5):
    sio = pytest.importorskip("scipy.io")
    t = np.arange(n + 1, dtype=float)
    layer = lambda k: {"brain": np.linspace(37, 37 + k, n), "skull": np.linspace(37, 38 + k, n)}
    seq = np.array([(0.0, "off", 1), (1.0, "on", 1), (2.0, "off", 1), (3.0, "on", 1), (4.0, "off", 1)],
                   dtype=[("time", "O"), ("status", "O"), ("recorded", "O")])
    path.parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(str(path), {
        "timeseries": {"T": layer(1), "Tdiff": layer(0), "CEM43": layer(0), "CEM43_iso": {}},
        "HEAT": {"recordedtime": t, "temp": {"Brain": np.full(n + 1, 37.0), "TARGET": t + 37}},
        "time_status_seq": seq})


def test_read_heat(tmp_path):
    path = tmp_path / "cache" / f"{TAG}_{MED}_HEAT.mat"
    _heat_file(path)
    heat = rs.read_heat(path)
    ts = heat.metrics["Temperature (tissue max)"]
    assert set(ts.tissues) == {"brain", "skull"}
    assert np.allclose(ts.time, [1, 2, 3, 4, 5])               # recordedtime without t=0
    assert heat.metrics["Temperature (focal axis)"].time.size == 6
    assert "CEM43 ISO (tissue max)" not in heat.metrics        # empty in this file
    assert heat.on == [(0.0, 1.0), (2.0, 3.0)]
    with pytest.raises(rs.HeatReadError):
        rs.read_heat(tmp_path / "missing.mat")


def test_heating_tab(schema, tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    pytest.importorskip("matplotlib")
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.results_page import ResultsPage
    app = QApplication.instance() or QApplication([])  # noqa: F841  keep it alive

    base = _folder(tmp_path)
    _heat_file(base / "cache" / f"{TAG}_{MED}_HEAT_desc-async_desc-asyncJoint.mat")
    page = ResultsPage(SimpleNamespace(model=_model(schema, tmp_path)))
    page.refresh()
    assert page.heat.data is None                              # not opened yet
    page.subtabs.setCurrentWidget(page.heat)
    assert page.heat.metric.currentText() == "Temperature (tissue max)"
    assert page.heat.shown_tissues() == ["brain", "skull"]
    page.heat.tissues.item(1).setCheckState(Qt.CheckState.Unchecked)
    assert page.heat.shown_tissues() == ["brain"]
    assert len(page.heat.ax.get_lines()) == 1
    page.heat.metric.setCurrentText("Temperature rise (tissue max)")
    assert page.heat.shown_tissues() == ["brain"]              # choice kept across metrics
    assert page.heat.shade.isEnabled()
    page.parts.setCurrentIndex(1)                              # transducer 1: the joint run's heat
    assert page.heat.data is not None
    page.run_combo.setCurrentText("coherent")
    assert page.heat.data is None and "not found" in page.heat.note.text()


LAYERED_CSV = {"subject_id": 3, "freq_Hz": 500000, "Isppa": 40.0, "Ispta": 4.0, "Isppa_brain": 12.5,
               "Isppa_skull": 30.0, "Isppa_skin": 5.0, "Psptp_brain": 0.6e6, "Psptp_skull": 1.1e6,
               "Psptp_skin": 0.3e6, "Psptp": 2.4e6, "MI_brain": 0.85, "MI_skull": 1.5, "MI_skin": 0.4,
               "MI_tc": 0.9, "real_focal_distance_ep_mm": 41.2, "focal_flhm_center_ep_mm": 43.0,
               "focal_vol_3dB_mm3": 120.0, "focal_vol_6dB_mm3": 260.0, "freefield_isppa_wcm2": 30.0,
               "TIC": 0.4, "riseT_brain": 0.3, "riseT_skull": 2.6, "riseT_skin": 0.1, "maxT_brain": 37.3,
               "maxT_skull": 39.6, "maxT_skin": 37.1, "CEM43_brain": 0.01, "CEM43_skull": 0.5,
               "CEM43_skin": 0.001, "CEM43iso_brain": 0.02, "CEM43iso_skull": 9.0, "CEM43iso_skin": 0.002}


def _write_csv(path, cols):
    path.write_text(",".join(cols) + "\n" + ",".join(str(v) for v in cols.values()) + "\n")


def test_itrusst_overview(schema, tmp_path):
    from prestus_nexus import itrusst_report as ir
    base = tmp_path / "sims" / TAG
    base.mkdir(parents=True)
    _write_csv(base / f"{TAG}_layered.csv", LAYERED_CSV)
    _write_csv(base / f"{TAG}_water.csv", {"subject_id": 3, "freq_Hz": 500000, "Isppa": 30.0, "Ispta": 3.0,
                                         "Psptp": 1.0e6, "MI": 1.4, "focal_vol_6dB_mm3": 300.0})
    _write_csv(base / f"{TAG}_layered_targeting.csv", {"placement_mode": "plantus", "focus_error_mm": 3.5,
                                                       "focus_error_axial_mm": 3.0, "focus_error_lateral_mm": 1.8})
    (base / f"{TAG}_placement.json").write_text(json.dumps(
        {"mode": "plantus", "transducers": {"trans_pos": [1, 2, 3], "focus_pos": [4, 5, 6]}}))
    m = ConfigModel(schema)
    for k, v in {"subject_id": 3, "path.sim": str(tmp_path / "sims"), "simulation.medium": "layered",
                 "timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 80, "timing.ptri": 80, "timing.ptrd": 80,
                 "transducer.freq_hz": 250000}.items():
        m.set(k, v)
    runs = rs.runs_for(m)
    run = next(r for r in runs if r.medium == "layered")
    rep_ = ir.build(m, run, runs)
    rows = {(sec.title, r.label): r for sec in rep_.sections for r in sec.rows}
    safety = {r.label: r for sec in rep_.sections if sec.title.startswith("ITRUSST") for r in sec.rows}
    # limits and colour rule of the HTML report: above > limit, near >= 50 %, within below
    assert safety["Max pressure"].status == "above" and safety["Max pressure"].value == "2.4 MPa"
    assert safety["MI (transcranial)"].status == "within" and safety["MI (skull)"].status == "near"
    assert safety["Temp rise (skull)"].status == "above" and safety["Max temp (skull)"].status == "above"
    assert safety["MI (free water, from the free-water run)"].value == "1.4"
    assert rep_.verdict == "review" and set(rep_.above) == {"Max pressure", "Temp rise (skull)", "Max temp (skull)"}
    # the frequency the run used (CSV), not the config's
    assert rows[("Transducer", "Fundamental frequency")].value == "500 kHz"
    assert rows[("Transducer", "Fundamental frequency")].source == "run CSV"
    # pulse timing derived from the config
    assert rows[("Pulse timing", "Pulse repetition frequency (PRF = 1/PRI)")].value == "5 Hz"
    assert rows[("Pulse timing", "Pulse duty cycle (PD/PRI)")].value == "10 %"
    assert rows[("Pulse timing", "Total sonication duration")].value == "80 s"
    # free water from the water run, in situ and targeting from their files
    fw = "Free-water acoustic parameters (simulated)"
    assert rows[(fw, "ISPPA")].value == "30 W/cm²"
    assert rows[(fw, "Focal width, axial (−3 dB FWHM)")].value == ir.NA      # no CSV column, no map
    assert rows[("Placement and targeting", "Placement method")].value == "plantus"
    assert rows[("Placement and targeting", "Planned focus (target)")].value == "[4, 5, 6] T1 voxel"
    assert rows[("Placement and targeting", "Focus error (focal centroid to planned focus)")].value == "3.5 mm"
    th = "Thermal (simulated)"
    assert rows[(th, "CEM43 (skull, k-Wave form)")].value == "0.5 min"
    m.set("thermal.cem43_iso", 1)                         # the ISO columns, as the report reads them
    rows = {(sec.title, r.label): r for sec in ir.build(m, run, runs).sections for r in sec.rows}
    assert rows[(th, "CEM43 (skull, ISO form)")].value == "9 min"
    assert rows[(th, "CEM43 (skull, ISO form)")].status == "near"
    md = ir.to_markdown(rep_)
    assert "| Max pressure | 2.4 MPa | ≤ 2 MPa | above | run CSV |" in md


def test_itrusst_overview_missing_outputs(schema, tmp_path):
    from prestus_nexus import itrusst_report as ir
    _folder(tmp_path)                                         # CSVs: Isppa_brain and riseT_brain only
    m = _model(schema, tmp_path)
    runs = rs.runs_for(m)
    rep_ = ir.build(m, runs[0], runs)
    safety = [r for r in rep_.sections[0].rows]
    assert [r.label for r in safety] == ["Temp rise (brain)"] and rep_.verdict == "within"
    rows = {r.label: r for sec in rep_.sections for r in sec.rows}
    assert rows["MI, transcranial (MItc)"].value == ir.NA and rows["MI, transcranial (MItc)"].status == ""
    assert rows["Manufacturer"].value == ir.NA
    assert any("free-water run" in sec.note for sec in rep_.sections)
