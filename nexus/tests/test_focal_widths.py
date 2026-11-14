"""Focal widths (focal_widths.py) and the per-tissue ISPTA of the Results overview."""

import json
import os

import numpy as np
import pytest

from prestus_nexus import focal_widths as fwm
from prestus_nexus import itrusst_report as ir
from prestus_nexus import results as rs
from prestus_nexus.model import ConfigModel

K = 2 * np.sqrt(2 * np.log(2))          # FWHM of exp(-x²/2σ²) is K σ


def gaussian(shape, centre, sigmas, rows=np.eye(3)):
    """exp(-Σ d_i²/2σ_i²), d_i along the rows of ``rows`` (voxel units)."""
    grid = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1) - np.asarray(centre)
    d = grid @ np.asarray(rows, dtype=float).T
    return np.exp(-(d ** 2 / (2 * np.asarray(sigmas, dtype=float) ** 2)).sum(axis=-1))


def test_crossing_width_interpolates_between_samples():
    # samples 1.0 at 0, 0.6 at ±1, 0.2 at ±2: half crossed at ±(1 + 0.1/0.4) = ±1.25
    prof = np.array([0.2, 0.6, 1.0, 0.6, 0.2])
    assert fwm.crossing_width(prof, 2, 1.0, 0.5) == pytest.approx(2.5)
    assert fwm.crossing_width(prof, 2, 0.5, 0.5) == pytest.approx(1.25)       # step in mm
    # one side never drops below half (or leaves the map): no width
    assert fwm.crossing_width(np.array([0.9, 1.0, 0.2]), 1, 1.0, 0.5) is None
    assert fwm.crossing_width(np.array([np.nan, 0.9, 1.0, 0.2]), 2, 1.0, 0.5) is None


def test_widths_axis_aligned_and_anisotropic():
    data = gaussian((41, 41, 81), (20, 20, 40), (2, 4, 8))
    ax, lo, hi = fwm.widths_at(data, np.eye(4), (20, 20, 40), (0, 0, 1))
    assert ax == pytest.approx(8 * K, rel=0.02)
    assert lo == pytest.approx(2 * K, rel=0.02) and hi == pytest.approx(4 * K, rel=0.02)
    # the axis direction's sign does not matter
    assert fwm.widths_at(data, np.eye(4), (20, 20, 40), (0, 0, -3))[0] == pytest.approx(ax)


def test_widths_oblique_axis_and_voxel_size():
    u = np.array([1.0, 0, 1]) / np.sqrt(2)
    rows = [np.cross([0, 1, 0], u), [0, 1, 0], u]
    data = gaussian((61, 61, 61), (30, 30, 30), (3, 3, 9), rows)
    ax, lo, hi = fwm.widths_at(data, np.eye(4), (30, 30, 30), u)
    assert ax == pytest.approx(9 * K, rel=0.02)
    assert lo == pytest.approx(3 * K, rel=0.03) and hi == pytest.approx(3 * K, rel=0.03)
    # 0.5 mm voxels: the same focus in voxels is half as wide in mm
    aff = np.diag([0.5, 0.5, 0.5, 1.0])
    ax2, lo2, _ = fwm.widths_at(data, aff, (30, 30, 30), u)
    assert ax2 == pytest.approx(ax / 2, rel=1e-6) and lo2 == pytest.approx(lo / 2, rel=1e-6)


def test_widths_truncated_lobe():
    data = gaussian((31, 31, 31), (15, 15, 1), (3, 3, 20))
    ax, lo, hi = fwm.widths_at(data, np.eye(4), (15, 15, 1), (0, 0, 1))
    assert ax is None and lo == pytest.approx(3 * K, rel=0.02)


def test_focal_peak_inside_mask():
    data = gaussian((21, 21, 21), (10, 10, 14), (2, 2, 4))
    data[10, 10, 2] = 5.0                                   # near-field hot spot outside the mask
    mask = np.zeros_like(data)
    mask[5:16, 5:16, 8:21] = 1
    assert fwm.focal_peak(data, mask) == (10, 10, 14)
    assert fwm.focal_peak(data, None) == (10, 10, 2)


def test_measure_on_nifti_is_cached(tmp_path):
    nib = pytest.importorskip("nibabel")
    data = gaussian((41, 41, 61), (20, 20, 35), (3, 3, 7)).astype(np.float32)
    mask = (data > 0.25).astype(np.uint8)
    data[20, 20, 3] = 4.0                                   # outside the focal mask
    aff = np.diag([0.8, 0.8, 0.8, 1.0])
    nib.save(nib.Nifti1Image(data, aff), str(tmp_path / "itp.nii.gz"))
    nib.save(nib.Nifti1Image(mask, aff), str(tmp_path / "mask.nii.gz"))
    trans, focus = [21, 21, 1], [21, 21, 41]                # 1-based T1 voxels
    fw = fwm.measure(tmp_path / "itp.nii.gz", tmp_path / "mask.nii.gz", trans, focus)
    assert fw.axial_mm == pytest.approx(7 * K * 0.8, rel=0.02)
    assert fw.lateral_min_mm == pytest.approx(3 * K * 0.8, rel=0.02)
    assert fw.source == fwm.SOURCE_MAP and fw.note == ""
    files = list((fwm._cache_file("x").parent).glob("*.json"))
    assert len(files) == 1 and json.loads(files[0].read_text())["axial_mm"] == pytest.approx(fw.axial_mm)
    fwm._memory.clear()                                     # a new session reads the disk cache
    axis = tuple(float(a) for a in (*trans, *focus))
    assert fwm.cached(tmp_path / "itp.nii.gz", tmp_path / "mask.nii.gz", axis) == fw
    # without the mask the near-field spot is the peak; said in the note
    fw2 = fwm.measure(tmp_path / "itp.nii.gz", None, trans, focus)
    assert "no focal mask" in fw2.note


def test_from_csv():
    assert fwm.from_csv({"Isppa": "1"}) is None
    fw = fwm.from_csv({"focal_fwhm_axial_mm": "30.5", "focal_fwhm_lateral_min_mm": "4",
                       "focal_fwhm_lateral_max_mm": "NaN"})
    assert (fw.axial_mm, fw.lateral_min_mm, fw.lateral_max_mm, fw.source) == (30.5, 4.0, None, "run CSV")


# ── per-tissue ISPTA ───────────────────────────────────────────────────────

def test_tissue_ispta_from_run_duty_cycle(schema):
    m = ConfigModel(schema)
    m.set("timing.pd", 0.05)
    m.set("timing.pri", 0.1)                                 # config says 50 %; the run used 10 %
    metrics = {"Isppa": "40", "Ispta": "4", "Isppa_brain": "12.5", "Isppa_skull": "30"}
    v, src = ir.tissue_ispta(metrics, "brain", m)
    assert v == pytest.approx(1.25) and "run CSV" in src
    assert ir.tissue_ispta(metrics, "skull", m)[0] == pytest.approx(3.0)
    assert ir.tissue_ispta(metrics, "skin", m)[0] is None    # no ISPPA for it


def test_tissue_ispta_falls_back_to_config(schema):
    m = ConfigModel(schema)
    m.set("timing.pd", 0.02)
    m.set("timing.pri", 0.2)
    v, src = ir.tissue_ispta({"Isppa_brain": "12.5", "Ispta": "NaN"}, "brain", m)
    assert v == pytest.approx(1.25) and "config" in src
    m.set("timing.pri", 0)
    assert ir.tissue_ispta({"Isppa_brain": "12.5"}, "brain", m)[0] is None


# ── overview rows ──────────────────────────────────────────────────────────

def _overview_folder(schema, tmp_path, extra=None):
    base = tmp_path / "sims" / "sub-003"
    (base / "nii").mkdir(parents=True)
    cols = {"subject_id": 3, "Isppa": 40.0, "Ispta": 4.0, "Isppa_brain": 12.5, "Isppa_skull": 30.0,
            "Isppa_skin": 5.0, **(extra or {})}
    (base / "sub-003_layered.csv").write_text(",".join(cols) + "\n" + ",".join(map(str, cols.values())) + "\n")
    m = ConfigModel(schema)
    m.set("subject_id", 3)
    m.set("path.sim", str(tmp_path / "sims"))
    m.set("simulation.medium", "layered")
    return base, m


def _rows(rep):
    return {r.label: r for sec in rep.sections if sec.title.startswith("Estimated") for r in sec.rows}


def test_overview_ispta_and_csv_widths(schema, tmp_path):
    base, m = _overview_folder(schema, tmp_path, {"focal_fwhm_axial_mm": 31.0,
                                                  "focal_fwhm_lateral_min_mm": 4.5,
                                                  "focal_fwhm_lateral_max_mm": 5.0})
    runs = rs.runs_for(m)
    rows = _rows(ir.build(m, runs[0], runs))
    assert rows["ISPTA (brain, pulse-train average, PD/PRI)"].value == "1.25 W/cm²"
    assert rows["ISPTA (skull, pulse-train average, PD/PRI)"].value == "3 W/cm²"
    assert rows["ISPTA (whole grid, pulse-train average, PD/PRI)"].value == "4 W/cm²"
    assert rows["Focal width, axial (−3 dB FWHM)"].value == "31 mm"
    assert rows["Focal width, axial (−3 dB FWHM)"].source == "run CSV"
    assert rows["Focal width, lateral (−3 dB FWHM), narrowest / widest"].value == "4.5 / 5 mm"
    assert ir.width_runs(runs[0], runs) == []                # nothing to measure


def test_overview_widths_measured_from_map(schema, tmp_path):
    nib = pytest.importorskip("nibabel")
    base, m = _overview_folder(schema, tmp_path)
    data = gaussian((31, 31, 41), (15, 15, 25), (2, 2, 5)).astype(np.float32)
    nib.save(nib.Nifti1Image(data, np.eye(4)), str(base / "nii" / "sub-003_layered_T1w_Itp.nii.gz"))
    runs = rs.runs_for(m)
    run = runs[0]
    assert ir.width_runs(run, runs) == []                    # no placement record: no axis
    (base / "sub-003_placement.json").write_text(json.dumps(
        {"transducers": [{"trans_pos": [16, 16, 1], "focus_pos": [16, 16, 26]}]}))
    assert ir.width_runs(run, runs) == [run]
    rows = _rows(ir.build(m, run, runs, {run.stem: ir.PENDING}))
    assert rows["Focal width, axial (−3 dB FWHM)"].value == ir.PENDING
    fw = fwm.measure(*fwm.inputs_for(run))
    rows = _rows(ir.build(m, run, runs, {run.stem: fw}))
    assert float(rows["Focal width, axial (−3 dB FWHM)"].value.split()[0]) == pytest.approx(5 * K, rel=0.02)
    assert rows["Focal width, axial (−3 dB FWHM)"].source.startswith(fwm.SOURCE_MAP)
    assert "no focal mask" in rows["Focal width, axial (−3 dB FWHM)"].source


def test_overview_panel_measures_off_the_gui_thread(schema, tmp_path):
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    nib = pytest.importorskip("nibabel")
    from PyQt6.QtWidgets import QApplication
    from prestus_nexus.ui.results_page import OverviewPanel
    app = QApplication.instance() or QApplication([])
    base, m = _overview_folder(schema, tmp_path)
    data = gaussian((31, 31, 41), (15, 15, 25), (2, 2, 5)).astype(np.float32)
    nib.save(nib.Nifti1Image(data, np.eye(4)), str(base / "nii" / "sub-003_layered_T1w_Itp.nii.gz"))
    (base / "sub-003_placement.json").write_text(json.dumps(
        {"transducers": [{"trans_pos": [16, 16, 1], "focus_pos": [16, 16, 26]}]}))
    runs = rs.runs_for(m)
    panel = OverviewPanel()
    panel.load(m, runs[0], runs)
    label = "Focal width, axial (−3 dB FWHM)"
    assert _rows(panel.report)[label].value == ir.PENDING
    import time
    deadline = time.monotonic() + 20
    while _rows(panel.report)[label].value == ir.PENDING and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.02)
    shown = _rows(panel.report)[label].value
    assert float(shown.split()[0]) == pytest.approx(5 * K, rel=0.02)
    panel.load(m, runs[0], runs)                             # cached: shown at once
    assert _rows(panel.report)[label].value == shown
