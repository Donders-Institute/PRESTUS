"""Beam axis tab on a synthetic head, offscreen: section, medium, free-water profiles."""

import math
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")
pytest.importorskip("matplotlib")
nib = pytest.importorskip("nibabel")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.ui.beam_view import BeamPage, profile_from_payload  # noqa: E402
from prestus_nexus.ui.placement_viewer import PlacementPage  # noqa: E402

GRID = (121, 121, 121)
AFFINE = np.diag([1.0, 1, 1, 1])
AFFINE[:3, 3] = -60


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _radius():
    g = np.indices(GRID).reshape(3, -1).T - 60
    return np.linalg.norm(g, axis=1).reshape(GRID)


@pytest.fixture
def subject(tmp_path):
    """1 mm T1 with a 50 mm 'head' sphere: brain, skull 40-45 mm, scalp 45-50 mm."""
    r = _radius()
    anat = tmp_path / "anat"
    (anat / "m2m_sub-004").mkdir(parents=True)
    nib.save(nib.Nifti1Image((r <= 50).astype(np.float32) * 100, AFFINE), str(anat / "sub-004_T1w.nii.gz"))
    seg = np.where(r <= 40, 2, np.where(r <= 45, 4, np.where(r <= 50, 5, 0))).astype(np.uint8)
    nib.save(nib.Nifti1Image(seg, AFFINE), str(anat / "m2m_sub-004" / "final_tissues.nii.gz"))
    return anat


def _pages(schema, prestus_root, anat, placed=True):
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for k, v in {"subject_id": 4, "path.anat": str(anat), "transducer.serial": "CTX_500_026",
                 "transducer.combo.ds_serial": "SC_203_035", "transducer.focal_distance_ep": 60,
                 "placement.manual.coordinate_space": "ras"}.items():
        m.set(k, v)
    viewer = PlacementPage(m)
    viewer.refresh()
    if placed:
        viewer._write("focus_pos", np.array([0.0, 0.0, 10.0]))
        viewer.direction = np.array([0.0, 0.0, 1.0])
        viewer._snap()                                  # exit plane 5 mm outside the scalp
    beam = BeamPage(viewer)
    beam.refresh()
    return m, viewer, beam


class _FakeBridge:
    """Synchronous stand-in for MatlabBridge: records calls, answers nexus_calibration."""

    def __init__(self, state="ready"):
        self.state, self.calls = state, []

    def call(self, fn, args=None, on_done=None, on_error=None, on_progress=None):
        self.calls.append(args["action"])
        if args["action"] == "open":
            on_done({"combo": args["serial"]})
            return
        z = [float(x) for x in range(0, 141)]
        imod = [math.exp(-((x - 60) / 15) ** 2) * 4.0 for x in z]
        on_done({"results": {"geometric": {"z": z, "Imod": [imod], "focus": 61.0},
                             "library": {"z": z, "Imod": [imod], "focus": 60.0}}})


def _click(beam, ax):
    beam._clicked(type("E", (), {"inaxes": ax, "button": 1, "xdata": 50, "ydata": 0.5})())


def test_empty_states_name_what_is_missing(app, schema, prestus_root, subject, tmp_path):
    m, viewer, beam = _pages(schema, prestus_root, subject, placed=False)
    assert beam.stack.is_empty() and beam.info.text() == ""
    m.set("path.anat", str(tmp_path / "nowhere"))
    viewer.refresh()
    beam.refresh()
    assert beam.stack.is_empty()


def test_section_medium_and_exit_plane_gap(app, schema, prestus_root, subject, tmp_path):
    """Exit plane -> medium, from the SimNIBS labels via layers, then from a run's medium."""
    m, viewer, beam = _pages(schema, prestus_root, subject)
    assert not beam.stack.is_empty()
    text = beam.info.text()
    assert "SimNIBS labels" in text and "Exit plane → medium" in text
    assert "+5." in text or "+4.5" in text          # 0.5 mm sampling of a 1 mm segmentation
    sim = tmp_path / "sim"
    (sim / "sub-004" / "nii").mkdir(parents=True)
    med = np.where(_radius() <= 53, 3, 1).astype(np.uint8)   # skin (index 3) swollen by 3 mm
    nib.save(nib.Nifti1Image(med, AFFINE), str(sim / "sub-004" / "nii" / "sub-004_layered_T1w_medium_masks.nii.gz"))
    m.set("path.sim", str(sim))
    beam.refresh()
    text = beam.info.text()
    assert "last run" in text and ("+2." in text or "+1.5" in text)
    beam._scrolled(type("E", (), {"inaxes": beam.beam_ax, "button": "up"})())
    beam.turn_right.click()
    assert beam.beam_angle == 30


def test_without_segmentation_the_section_still_draws(app, schema, prestus_root, subject):
    (subject / "m2m_sub-004" / "final_tissues.nii.gz").unlink()
    m, viewer, beam = _pages(schema, prestus_root, subject, placed=False)
    viewer._write("trans_pos", np.array([0.0, 0.0, 70.0]))
    viewer._write("focus_pos", np.array([0.0, 0.0, 10.0]))
    beam.refresh()
    assert not beam.stack.is_empty()
    assert "Segmentation" in beam.info.text() and "no profile yet" in beam.info.text()


def test_prediction_is_relative_and_computed_on_request(app, schema, prestus_root, subject):
    m, viewer, beam = _pages(schema, prestus_root, subject)
    bridge = _FakeBridge(state="off")
    beam.bridge_getter = lambda: bridge
    beam.refresh()
    assert bridge.calls == []                       # MATLAB not running: a redraw never starts it
    assert not beam.compute.isHidden() and not beam.rb_absolute.isEnabled()
    beam.compute.click()
    assert bridge.calls == ["open", "view"]
    assert "Free-water focus" in beam.info.text() and "60.0 mm" in beam.info.text()
    assert beam.compute.isHidden()
    beam.refresh()
    assert bridge.calls == ["open", "view"]         # cached per device and focal depth
    m.set("transducer.focal_distance_ep", 50)
    bridge.state = "ready"                          # running: computed for the new depth at once
    beam.refresh()
    assert len(bridge.calls) == 4


def test_run_profile_absolute_and_stale(app, schema, prestus_root, subject, tmp_path):
    """A run's axial_profile CSV (plot_axial_profile.m) gives absolute profiles, no MATLAB."""
    m, viewer, beam = _pages(schema, prestus_root, subject)
    viewer._match_focal_distance()
    ep = m.get("transducer.focal_distance_ep")
    img = tmp_path / "sim" / "sub-004" / "img"
    img.mkdir(parents=True)
    rows = ["depth_from_exit_plane_mm,pressure_Pa,intensity_Wcm2,tissue,freewater_oneil_Wcm2"]
    for d in range(-20, 120, 2):
        oneil = "NaN" if d < 0 else f"{30 * math.exp(-((d - ep) / 12) ** 2):.4f}"
        rows.append(f"{d},1,{5 * math.exp(-((d - ep + 3) / 10) ** 2):.4f},brain,{oneil}")
    (img / "sub-004_layered_axial_profile.csv").write_text("\n".join(rows) + "\n")
    (tmp_path / "sim" / "sub-004" / "sub-004_layered_targeting.csv").write_text(
        f"subject_id,focal_distance_ep_planned_mm\n4,{ep}\n")
    bridge = _FakeBridge(state="ready")
    beam.bridge_getter = lambda: bridge
    m.set("path.sim", str(tmp_path / "sim"))
    beam.refresh()
    assert beam.rb_absolute.isEnabled() and beam.compute.isHidden()
    assert "Free-water focus" in beam.info.text() and "another placement" not in beam.info.text()
    _click(beam, beam.prof_ax)
    assert bridge.calls == []                       # nothing to compute in MATLAB
    beam.rb_absolute.click()
    assert beam.scale == "absolute" and "Simulated peak in brain" in beam.info.text()
    m.set("transducer.focal_distance_ep", ep + 10)   # settings changed since the run
    beam.refresh()
    assert "another placement or focal depth" in beam.info.text()


def test_shown_tab_reads_run_outputs_in_the_background(app, schema, prestus_root, subject,
                                                      tmp_path, monkeypatch):
    """With the window shown, refresh() draws at once from what it has and reads the run's
    medium off the GUI thread, then redraws with it."""
    import threading
    import time

    from prestus_nexus.ui import beam_view
    m, viewer, beam = _pages(schema, prestus_root, subject)
    nii = tmp_path / "sim" / "sub-004" / "nii"
    nii.mkdir(parents=True)
    medium = np.where(_radius() <= 50, 2, 1).astype(np.uint8)          # 1 = water
    nib.save(nib.Nifti1Image(medium, AFFINE), str(nii / "sub-004_layered_T1w_medium_masks.nii.gz"))
    m.set("path.sim", str(tmp_path / "sim"))
    gate, threads = threading.Event(), []
    real = beam_view._read_run_outputs

    def slow(*args):
        threads.append(threading.current_thread())
        gate.wait(10)                                   # the share is slow
        return real(*args)
    monkeypatch.setattr(beam_view, "_read_run_outputs", slow)
    beam.show()
    try:
        beam.refresh()                                  # returns without waiting for the read
        assert beam.medium is None and "SimNIBS labels" in beam.info.text()
        gate.set()
        t0 = time.monotonic()
        while beam.medium is None and time.monotonic() - t0 < 10:
            app.processEvents()
            time.sleep(0.02)
        assert beam.medium is not None and threads[0] is not threading.main_thread()
        assert "last run" in beam.info.text()
    finally:
        beam.hide()


def test_profile_from_payload_reads_one_depth():
    z = [0.0, 1.0, 2.0]
    out = {"results": {"geometric": {"z": z, "Imod": [[1.0], [3.0], [2.0]], "focus": 1.0}}}
    prof = profile_from_payload(out)              # one depth sent as a flat column
    assert prof["model"] == "geometric" and list(prof["I"]) == pytest.approx([1 / 3, 1, 2 / 3])
    assert profile_from_payload({"results": {}}) is None and profile_from_payload(None) is None
