"""Grid orientation under grid.mode (ports of mean_beam_axis.m, preproc_head.m and
preproc_align_to_mean_focal_axis.m), against hand-computed cases."""

import math

import numpy as np
import pytest

from prestus_nexus import multitransducer as mt
from prestus_nexus.model import ConfigModel


def tx(trans, focus, **kw):
    return {"trans_pos": list(trans), "focus_pos": list(focus), **kw}


def test_single_beam_is_its_own_axis():
    txs = [tx([0, 0, 0], [0, 0, 10])]
    for mode in ("focal_axis", "transducer_axis"):
        z, fell = mt.grid_z_axis(txs, mode)
        np.testing.assert_allclose(z, [0, 0, 1])
        assert not fell
    assert mt.beam_tilts(txs, [0, 0, 1]) == [pytest.approx(0)]


def test_two_beams_same_side_give_bisector():
    # both 40 deg from -z, converging from above (80 deg apart): bisector -z, each tilted 40 deg
    s, c = math.sin(math.radians(40)), math.cos(math.radians(40))
    txs = [tx([-10 * s, 0, 10 * c], [0, 0, 0]), tx([10 * s, 0, 10 * c], [0, 0, 0])]
    z, fell = mt.grid_z_axis(txs, "focal_axis")
    np.testing.assert_allclose(z, [0, 0, -1], atol=1e-12)
    assert not fell
    assert mt.beam_tilts(txs, z) == [pytest.approx(40), pytest.approx(40)]
    # transducer_axis: beam 1 on z, beam 2 at the full 80 deg
    z1, _ = mt.grid_z_axis(txs, "transducer_axis")
    np.testing.assert_allclose(z1, [s, 0, -c])
    assert mt.beam_tilts(txs, z1) == [pytest.approx(0, abs=1e-6), pytest.approx(80)]


def test_two_opposed_beams_tilt_by_half_the_axis_angle():
    # beam 1 along +x; beam 2 points against it, its axis 60 deg from beam 1's axis
    c, s = math.cos(math.radians(60)), math.sin(math.radians(60))
    txs = [tx([-10, 0, 0], [0, 0, 0]), tx([10 * c, 0, -10 * s], [0, 0, 0])]
    z, fell = mt.grid_z_axis(txs, "focal_axis")
    assert not fell
    np.testing.assert_allclose(z, [math.cos(math.radians(30)), 0, -0.5], atol=1e-12)
    assert z @ np.array([1, 0, 0]) > 0                    # points along transducer 1's beam
    assert mt.beam_tilts(txs, z) == [pytest.approx(30), pytest.approx(30)]


def test_three_spread_beams_fall_back_to_transducer_1():
    # axes x, +y, -y: the y beams cancel, mean norm 1/3 < 0.5
    txs = [tx([-10, 0, 0], [0, 0, 0]), tx([0, -10, 0], [0, 0, 0]), tx([0, 10, 0], [0, 0, 0])]
    _, norm, _ = mt.mean_beam_axis(mt.beam_units(txs))
    assert norm == pytest.approx(1 / 3)
    z, fell = mt.grid_z_axis(txs, "focal_axis")
    assert fell
    np.testing.assert_allclose(z, [1, 0, 0])
    # grid.focal_axis_min_norm below the norm keeps the mean axis (here also x)
    z, fell = mt.grid_z_axis(txs, "focal_axis", min_norm=0.3)
    assert not fell
    np.testing.assert_allclose(z, [1, 0, 0])
    # orthogonal triple: norm 1/sqrt(3) = 0.577 > 0.5 keeps the mean, 0.6 falls back
    ortho = [tx([-1, 0, 0], [0, 0, 0]), tx([0, -1, 0], [0, 0, 0]), tx([0, 0, -1], [0, 0, 0])]
    z, fell = mt.grid_z_axis(ortho)
    assert not fell
    np.testing.assert_allclose(z, np.ones(3) / math.sqrt(3))
    assert mt.grid_z_axis(ortho, min_norm=0.6)[1]


def test_ras_plus_keeps_ras_z():
    txs = [tx([-10, 0, 10], [0, 0, 0]), tx([10, 0, 10], [0, 0, 0])]     # 90 deg apart
    z, fell = mt.grid_z_axis(txs, "ras_plus")
    np.testing.assert_allclose(z, [0, 0, 1])
    assert not fell
    assert mt.beam_tilts(txs, z) == [pytest.approx(45), pytest.approx(45)]


def test_removed_and_unknown_modes():
    txs = [tx([0, 0, 0], [0, 0, 10])]
    with pytest.raises(mt.GridModeError, match="principal_axis"):
        mt.grid_z_axis(txs, "principal_axis")
    np.testing.assert_allclose(mt.grid_z_axis(txs, "nonsense")[0], [0, 0, 1])   # -> focal_axis


def test_align_to_focus_false_uses_plus_z_for_transducer_1():
    txs = [tx([0, 0, 0], [10, 0, 0], align_to_focus=False), tx([0, 0, 0], [0, 0, 5])]
    z, _ = mt.grid_z_axis(txs, "transducer_axis")
    np.testing.assert_allclose(z, [0, 0, 1])


def test_missing_positions_raise():
    with pytest.raises(mt.GridModeError, match="transducer 2"):
        mt.grid_z_axis([tx([0, 0, 0], [0, 0, 1]), {"trans_pos": [], "focus_pos": []}])
    # transducer_axis only needs transducer 1
    z, _ = mt.grid_z_axis([tx([0, 0, 0], [0, 0, 1]), {"trans_pos": [], "focus_pos": []}],
                          "transducer_axis")
    np.testing.assert_allclose(z, [0, 0, 1])


def test_orientation_from_model_flags_tilted_beams_on_coarse_grids(schema):
    m = ConfigModel(schema)
    m.set("transducer.trans_pos", [-10, 0, 10])
    m.set("transducer.focus_pos", [0, 0, 0])
    m.set("transducer.freq_hz", 500000)
    m.add_transducer()
    m.set("transducer.trans_pos", [10, 0, 10], 1)
    m.set("transducer.focus_pos", [0, 0, 0], 1)
    m.set("grid.resolution_mm", 1.0)                 # 3 points per wavelength in water
    o = mt.orientation(m)
    assert not o["error"]
    assert o["tilts"] == [pytest.approx(45), pytest.approx(45)]
    assert o["coarse_tilted"] == [0, 1]
    m.set("grid.resolution_mm", 0.4)                 # 7.5 points per wavelength
    assert mt.orientation(m)["coarse_tilted"] == []
    m.set("grid.mode", "principal_axis")
    assert "removed" in mt.orientation(m)["error"]


# ── top-view schematic ────────────────────────────────────────────────────────

def test_schematic_top_view_left_right_front():
    """Left and right stances land on the left and right, the farthest on the outline, the
    foci inside; z is dropped; a missing position gives None."""
    txs = [{"trans_pos": [40, 100, 150], "focus_pos": [90, 100, 100]},     # left
           {"trans_pos": [160, 100, 150], "focus_pos": [110, 100, 100]},   # right
           {"trans_pos": [100, 170, 120], "focus_pos": [100, 120, 100]},   # front
           {"trans_pos": None, "focus_pos": None}]
    out = mt.schematic_top_view(txs)
    (l_t, l_f), (r_t, r_f), (a_t, a_f), missing = out
    assert missing is None
    assert l_t[0] < 0 < r_t[0] and a_t[1] > 0
    ell = lambda v: math.hypot(v[0] / mt.HEAD_RX, v[1] / mt.HEAD_RY)   # noqa: E731
    assert max(ell(l_t), ell(r_t), ell(a_t)) == pytest.approx(1.0)
    assert all(ell(f) < 1 for f in (l_f, r_f, a_f))


def test_schematic_top_view_degenerate():
    assert mt.schematic_top_view([{"trans_pos": None}]) == [None]
    same = mt.schematic_top_view([{"trans_pos": [1, 2, 3], "focus_pos": [1, 2, 3]}])
    assert np.allclose(same[0][0], 0) and np.allclose(same[0][1], 0)     # no division by zero
    no_focus = mt.schematic_top_view([{"trans_pos": [0, 0, 0]}, {"trans_pos": [10, 0, 0]}])
    assert no_focus[0][1] is None and no_focus[0][0][0] < 0 < no_focus[1][0][0]
