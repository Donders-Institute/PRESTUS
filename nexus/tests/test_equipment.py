"""Equipment library, calibration coverage and the model's library layer."""

import pytest

from prestus_nexus import checks
from prestus_nexus.equipment import EquipmentLibrary
from prestus_nexus.model import ConfigModel


@pytest.fixture(scope="module")
def library(prestus_root):
    return EquipmentLibrary(prestus_root)


@pytest.fixture
def model(schema, library):
    m = ConfigModel(schema, library)
    m.set("transducer.serial", "CTX_500_026")
    m.set("transducer.combo.ds_serial", "SC_203_035")
    return m


def test_devices_load(library):
    dev = library.get("CTX_500_026")
    assert dev is not None and dev.type == "annular" and dev.freq_hz == 500e3
    assert dev.min_foc < dev.max_foc
    assert "SC_203_035" in dev.combos


def test_calibration_lookup_uses_combo_name(library, tmp_path):
    cal = library.calibration("CTX_500_026", "SC_203_035")
    assert cal.found and not cal.fallback and cal.path.name == "CTX_500_026_SC_203_035.yaml"
    assert cal.intensities                              # e.g. i30 -> 30 W/cm²
    assert not library.calibration("CTX_500_026", "NOPE").found
    # a user library with an entry takes precedence (transducer.library_path)
    (tmp_path / "CTX_500_026.yaml").write_text("meta: {calibration_type: measured}\n")
    user = library.calibration("CTX_500_026", "NOPE", str(tmp_path))
    assert user.found and user.fallback and user.calibration_type == "measured"


def test_library_values_are_the_baseline(model, library):
    geom = library.get("CTX_500_026").geometry
    assert model.get("transducer.type") == "annular"
    assert model.get("transducer.freq_hz") == geom["freq_hz"]
    assert model.get("transducer.annular.elem_od_mm") == geom["annular"]["elem_od_mm"]
    assert model.from_library("transducer.annular.curv_radius_mm")
    # setting the library value again adds no override
    model.set("transducer.freq_hz", geom["freq_hz"])
    assert "freq_hz" not in model.tx[0]


def test_inline_values_override_library(model):
    model.set("transducer.freq_hz", 250e3)
    assert model.get("transducer.freq_hz") == 250e3
    assert not model.from_library("transducer.freq_hz")
    assert model.effective()["transducer"]["freq_hz"] == 250e3
    # an empty study value counts as unset (resolve_transducer_from_serial)
    model.tx[0]["freq_hz"] = []
    assert model.get("transducer.freq_hz") == 500e3


def test_library_checks(model):
    model.set("transducer.focal_distance_ep", 95)
    found = {(i.level, i.path) for i in checks.library_issues(model)}
    assert ("warning", "transducer.focal_distance_ep") in found     # above max_foc
    model.set("transducer.combo.ds_serial", "NOPE")
    found = {(i.level, i.path) for i in checks.library_issues(model)}
    assert ("error", "transducer.combo.ds_serial") in found         # no calibration entry
    model.set("transducer.serial", "NOT_A_DEVICE")
    found = {(i.level, i.path) for i in checks.library_issues(model)}
    assert ("error", "transducer.serial") in found


def test_lookup_needs_focal_distance(model):
    found = {(i.level, i.path) for i in checks.library_issues(model)}
    assert ("error", "transducer.focal_distance_ep") in found
    model.set("transducer.focal_distance_ep", 60)
    assert not [i for i in checks.library_issues(model) if i.level == "error"]


def test_depths_mm_reads_the_current_schema(library):
    """Regression guard: depths_mm must read refinement.depths_ep_mm (the schema every
    real per-depth entry in config/calibration/ uses), not global_model.depths_ep_mm
    (which doesn't exist there) -- the latter silently returned [] for the whole
    per-depth library."""
    cal = library.calibration("CTX_250_001", "SC_105_010")
    assert cal.found and cal.calibration_type == "empirical"
    assert cal.depths_mm, "depths_mm is empty: the field read has drifted from the real schema again"
    assert cal.fit_depths_mm and len(cal.fit_depths_mm) == len(cal.fit_r2)
    assert cal.blend_mm == 5.0


def test_at_depth_in_range_and_extrapolated(library):
    cal = library.calibration("CTX_250_001", "SC_105_010")
    lo, hi = min(cal.depths_mm), max(cal.depths_mm)
    mid = cal.at_depth((lo + hi) / 2)
    assert mid.in_range and mid.phase_model == "refined" and mid.bracket is not None
    far = cal.at_depth(hi + 100)
    assert not far.in_range and far.phase_model == "static"
    assert far.nearest_mm == hi and far.extrapolation_mm == pytest.approx(100)
    just_past = cal.at_depth(hi + cal.blend_mm / 2)
    assert just_past.phase_model == "blend"


def test_at_depth_without_calibration_is_none():
    from prestus_nexus.equipment import Calibration
    cal = Calibration(path=None, key="nope")
    assert cal.at_depth(50) is None
