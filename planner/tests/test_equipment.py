"""Equipment library, calibration coverage and the model's library layer."""

import pytest

from prestus_planner import checks
from prestus_planner.equipment import EquipmentLibrary
from prestus_planner.model import ConfigModel


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
