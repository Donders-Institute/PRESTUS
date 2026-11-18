"""calibration_plots: the arrays read from library entries, and the widgets offscreen."""

import math
import os
from pathlib import Path

import pytest

from prestus_nexus import yamlio

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from prestus_nexus.ui import calibration_plots as cp  # noqa: E402

CAL = Path(__file__).resolve().parents[2] / "config" / "calibration"


def test_phase_unwrap_turns_wraparound_into_a_small_step():
    entry = {"refinement": {"depths_ep_mm": [30, 40, 50],
                            "elem_phase_correction_deg": [[354, 10], [26, 12], [20, 14]]}}
    a = cp.entry_arrays(entry)
    ring1 = a["phase_deg"][0]
    assert ring1[0] == pytest.approx(-6) and ring1[1] - ring1[0] == pytest.approx(32)
    assert [(j, c) for j, c, _ in a["phase_jumps"]] == [(1, 0)]          # 32° > 30° marked
    assert a["cal_range"] == (30, 50) and a["drive_source"] == ""


def test_sub_ring_columns_merge_and_status_follows_blend():
    entry = {"meta": {"quality": {"geometric_r2_median": 0.7, "in_sample_r2_median": 0.97, "flags": []}},
             "global_model": {"fit_depths_ep_mm": [40, 50], "fit_r2": [0.8, 0.95],
                              "amplitude_by_depth": {"depths_ep_mm": [40, 50], "elem_p_amp": [9e4, 1e5],
                                                     "isppa_wcm2": 30}},
             "refinement": {"depths_ep_mm": [40, 50], "elem_p_amp": [1e5, 1.1e5], "isppa_wcm2": 30,
                            "elem_phase_correction_deg": [[1, 1, 5], [2, 2, 6]], "blend_mm": 5}}
    a = cp.entry_arrays(entry)
    assert len(a["phase_deg"]) == 2 and a["channel_labels"] == ["ring 1", "ring 2"]
    assert a["drive_source"] == "refinement" and a["drive_pa"] == [1e5, 1.1e5]
    assert [cp.depth_status(a, d) for d in (45, 53, 60)] == ["within", "blend", "outside"]
    text = cp.summary_text(a, [45])
    assert text.startswith("2 depths, 40.0–50.0 mm · median R² 0.97") and "45 mm: within range" in text
    assert cp.r2_colour_key(0.95) == "green" and cp.r2_colour_key(0.85) == "amber" \
        and cp.r2_colour_key(0.5) == "red" and cp.r2_colour_key(math.nan) == "ink4"


def test_empty_and_legacy_entries_degrade():
    a = cp.entry_arrays(None)
    assert a["cal_range"] is None and a["phase_deg"] == [] and cp.summary_text(a) == \
        "no calibrated depths in this entry"
    legacy = cp.entry_arrays({"global_model": {"depths_ep_mm": [35, 45], "elem_phase_correction_deg": [1, 2]}})
    assert legacy["cal_range"] == (35, 45) and legacy["phase_depths"] == []


@pytest.mark.parametrize("name", ["CTX_500_026_SC_203_035", "CTX_250_009_SC_203_035",
                                  "IS_PCD15287_01001_IGT_128_ch_comb_10_ch",
                                  "IS_PCD15473_01003_IGT_128_ch_comb_10_ch"])
def test_shipped_entries(name):
    path = CAL / f"{name}.yaml"
    if not path.is_file():
        pytest.skip("shipped calibration entry not available")
    a = cp.entry_arrays(yamlio.load_file(path))
    if a["calibration_type"] == "geometric":
        assert a["cal_range"] is None
        return
    assert a["cal_range"] and len(a["fit_depths"]) == len(a["fit_r2"])
    assert all(abs(c[j] - c[j - 1]) <= 180 for c in a["phase_deg"] for j in range(1, len(c)))
    if name.startswith("IS_PCD15287"):
        assert len(a["phase_deg"]) == 10                   # 30 sub-ring columns, 10 rings


def test_widget_offscreen():
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    w = cp.CalibrationCoverage()
    w.resize(420, 320)
    assert w.stack.currentWidget() is w.empty and "No calibration entry" in w.empty.text()
    entry = {"meta": {"quality": {"in_sample_r2_median": 0.97}},
             "global_model": {"fit_depths_ep_mm": [40, 50], "fit_r2": [0.85, 0.95]},
             "refinement": {"depths_ep_mm": [40, 50], "elem_p_amp": [1e5, 1.1e5], "isppa_wcm2": 30,
                            "elem_phase_correction_deg": [[354, 1], [26, 2]]}}
    w.set_entry(entry, [45, 70], "test")
    if w.canvas is not None:
        assert w.stack.currentWidget() is w.canvas
        for view in cp.VIEWS:
            w.set_view(view)
            w.canvas.draw()
    assert "70 mm" in w.summary_text() and "outside" in w.summary_text()
    s = cp.R2Strip()
    s.resize(400, 70)
    got = []
    s.selected.connect(got.append)
    s.set_data([40, 50], [0.95, 0.7], allow_all=True, caption="the fit")
    s.grab()
    _all, rects = s._geometry()
    from PyQt6.QtCore import QEvent, QPointF, Qt
    from PyQt6.QtGui import QMouseEvent
    ev = QMouseEvent(QEvent.Type.MouseButtonPress, QPointF(rects[1].center()), QPointF(rects[1].center()),
                     Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier)
    s.mousePressEvent(ev)
    assert got == [50.0] and s.current == 50.0
    app.processEvents()
