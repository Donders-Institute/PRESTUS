import math

import numpy as np
import pytest

from prestus_nexus import coupling as cp


def _beam(apex, focus, isppa=30.0):
    return cp.Beam(np.array(apex, float), np.array(focus, float), aperture_mm=60.0,
                   freq_hz=500e3, isppa_wcm2=isppa)


def test_at_target_range_over_relative_phase():
    t = cp.at_target([30, 30])
    assert t["independent"] == pytest.approx(60)              # intensities add
    assert t["coherent_max"] == pytest.approx(120)            # (√I1 + √I2)²: 4× one beam
    assert t["coherent_min"] == pytest.approx(0)              # equal amplitudes cancel
    t = cp.at_target([30, 10, 5])
    a = np.sqrt([30, 10, 5])
    assert t["coherent_min"] == pytest.approx((a[0] - a[1] - a[2]) ** 2)


def test_spot_widths_and_fringes():
    b = _beam([0, 0, -60], [0, 0, 0])
    lat, ax = b.widths_mm()
    lam = 1500 / 500e3 * 1e3                                   # 3 mm in water
    assert lat == pytest.approx(1.02 * lam * 1.0) and ax == pytest.approx(7 * lam * 1.0)
    assert cp.fringe_spacing_mm(90, lam) == pytest.approx(lam / (2 * math.sin(math.radians(45))))
    assert cp.fringe_spacing_mm(0, lam) is None                # parallel beams: no fringes


def test_maps_coherent_peak_and_independent_sum():
    beams = [_beam([-60, 0, 0], [0, 0, 0]), _beam([0, -60, 0], [0, 0, 0])]   # 90° apart, same focus
    assert cp.beam_angle_deg(beams) == pytest.approx(90)
    m = cp.field_maps(beams, extent_mm=10, n=101)
    c = 50                                                     # the shared focus
    assert m["independent"][c, c] == pytest.approx(60)
    assert m["coherent"][c, c] == pytest.approx(120)           # in phase at the foci
    assert np.any(m["coherent"] < 0.5 * m["independent"])        # destructive fringes nearby
