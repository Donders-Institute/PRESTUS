"""Equipment library, calibration coverage and the model's library layer."""

import pytest

from prestus_nexus import checks, equipment
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


def test_models_listed_before_their_units(library):
    order = list(library.devices)
    shipped = [s for s, d in library.devices.items() if not d.user]
    units = 0
    for serial in shipped:
        dev = library.devices[serial]
        group = equipment.model_of(serial, dev.model, dev.file, library.devices)
        if group != serial:
            units += 1
            assert order.index(group) < order.index(serial), serial
            between = order[order.index(group):order.index(serial)]
            assert all(equipment.model_of(s, library.devices[s].model, library.devices[s].file,
                                          library.devices) == group for s in between), serial
    assert units > 0                                                  # the check ran on units
    i = order.index("CTX_500")                                        # a unit without 'model:'
    assert order[i:i + 3] == ["CTX_500", "CTX_500_006", "CTX_500_026"]   # joins its folder's model


def test_models_first_orders_any_group():
    D = equipment.Device
    devs = {s: D(serial=s, name="", manufact="", n_elem=1, min_foc=None, max_foc=None,
                 geometry={}, combos={}, model=m, file=f"/x/{folder}/tran-{s}_transducer.yaml")
            for s, m, folder in [("A_1", "A", "model-A"), ("A", "A", "model-A"),
                                 ("A_2", "", "model-A"), ("B", "", "model-B"),
                                 ("C_9", "C", "flat"), ("Z", "", "flat")]}
    assert list(equipment.models_first(devs)) == ["A", "A_1", "A_2", "B", "C_9", "Z"]


def test_calibration_lookup_uses_combo_name(library, tmp_path):
    cal = library.calibration("CTX_500_026", "SC_203_035")
    assert cal.found and not cal.fallback and cal.path.name == "tran-CTX_500_026_ds-SC_203_035_calib.yaml"
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


def test_parsed_entries_are_cached(prestus_root, monkeypatch):
    """A second start reuses the parsed equipment YAMLs (network drives: ~8 s)."""
    first = EquipmentLibrary(prestus_root)
    from prestus_nexus import yamlio

    def no_read(path):
        raise AssertionError(f"read {path} although it is unchanged")
    monkeypatch.setattr(yamlio, "load_file", no_read)
    second = EquipmentLibrary(prestus_root)
    assert set(second.devices) == set(first.devices) and first.devices


def test_changed_entry_is_read_again(tmp_path):
    eq = tmp_path / "config" / "equipment"
    eq.mkdir(parents=True)
    entry = eq / "TX_1.yaml"
    entry.write_text("type: 'transducer'\nserial: 'TX_1'\nname: 'first'\n")
    assert EquipmentLibrary(tmp_path).get("TX_1").name == "first"
    entry.write_text("type: 'transducer'\nserial: 'TX_1'\nname: 'second one'\n")
    assert EquipmentLibrary(tmp_path).get("TX_1").name == "second one"


def test_unknown_unit_falls_back_to_its_model(library):
    assert library.model_for("CTX_500_031") == "CTX_500"
    assert library.model_for("CTX_250_2ch_007") == "CTX_250_2ch"      # longest model prefix
    assert library.model_for("CTX_500_026") is None                   # in the library itself
    assert library.model_for("CTX_500") is None
    assert library.model_for("NOT_A_DEVICE_001") is None
    assert library.lookup("CTX_500_031").serial == "CTX_500"
    assert library.lookup("CTX_500_026").serial == "CTX_500_026"
    for serial, dev in library.devices.items():                       # every entry names its model
        if serial != "CTX_500_006":
            assert dev.model in library.devices, serial


def test_model_fallback_is_a_warning_and_resolves_geometry(schema, library):
    m = ConfigModel(schema, library)
    m.set("transducer.serial", "CTX_500_031")
    m.set("transducer.combo.ds_serial", "SC_203_035")
    m.set("transducer.focal_distance_ep", 50)
    assert m.effective()["transducer"]["annular"]["curv_radius_mm"] == 62.94
    issues = checks.library_issues(m)
    found = {(i.level, i.path) for i in issues}
    assert ("warning", "transducer.serial") in found
    assert not any(i.level == "error" for i in issues), [i.message for i in issues]


def test_unit_inherits_model_steering_range(library):
    model, unit = library.get("IS_PCD15287"), library.get("IS_PCD15287_01002")
    assert (unit.min_foc, unit.max_foc) == (model.min_foc, model.max_foc)
    assert unit.min_foc < 7.1 and unit.max_foc > 93.92     # hardware range spans the measured depths


def test_deprecated_device_is_labelled(library):
    from dataclasses import replace
    dev = replace(library.get("CTX_500_026"), status="deprecated", superseded_by="CTX_500")
    assert dev.deprecated and dev.label().endswith("(deprecated)")


def test_excluded_device_is_labelled(library):
    from dataclasses import replace
    dev = replace(library.get("IS_PCD15473_01001"), status="excluded")
    assert dev.excluded and dev.label().endswith("(excluded)")
    assert not library.get("IS_PCD15473_01001").excluded


def test_own_library_folder_adds_equipment(schema, prestus_root, tmp_path):
    """Equipment entries in transducer.library_path add a unit that inherits its model."""
    own = tmp_path / "my_library"
    own.mkdir()
    (own / "tran-IS_PCD15473_09999_transducer.yaml").write_text(
        "type: transducer\nserial: IS_PCD15473_09999\nmodel: IS_PCD15473\nname: my unit\n"
        "combos:\n  IGT_32_ch_comb_10_ch: {}\n")
    lib = EquipmentLibrary(prestus_root)
    m = ConfigModel(schema, lib)
    m.set("transducer.serial", "IS_PCD15473_09999")
    assert m.device(0).serial == "IS_PCD15473"           # not shipped: falls back to its model
    m.set("transducer.library_path", str(own))
    dev = m.device(0)
    assert dev is not None and dev.user and dev.label().endswith("(own library)")
    assert dev.geometry == lib.get("IS_PCD15473").geometry     # from the shipped model
    assert dev.combos["IGT_32_ch_comb_10_ch"].available_ch == 32  # from ds/
    m.set("transducer.library_path", "")
    assert m.device(0).serial == "IS_PCD15473"


def _unit(folder, name):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "tran-IS_PCD15473_09998_transducer.yaml").write_text(
        f"type: transducer\nserial: IS_PCD15473_09998\nmodel: IS_PCD15473\nname: {name}\n")


def test_own_libraries_search_order(schema, prestus_root, tmp_path, monkeypatch):
    """Study library_path before path.transducer_library before PRESTUS_TRANSDUCER_LIBRARY."""
    study, project, machine = tmp_path / "study", tmp_path / "project", tmp_path / "machine"
    for folder in (study, project, machine):
        _unit(folder, folder.name)
    (machine / "profiles").mkdir()
    (machine / "profiles" / "p.csv").write_text("x")
    monkeypatch.setenv("PRESTUS_TRANSDUCER_LIBRARY", str(machine))
    lib = EquipmentLibrary(prestus_root)
    m = ConfigModel(schema, lib)
    m.set("transducer.serial", "IS_PCD15473_09998")
    assert m.device(0).name == "machine"
    m.set("path.transducer_library", str(project))
    assert m.device(0).name == "project"
    m.set("transducer.library_path", str(study))
    assert m.device(0).name == "study"
    assert m.user_library_folders() == [study, project, machine]
    assert lib.profile_file("p.csv", "IS_PCD15473_09998") == machine / "profiles" / "p.csv"
    assert lib.profile_file("nope.csv") is None


def test_short_names(prestus_root):
    """Model then unit number for library devices and driving systems; unknown serials and
    serials without a unit number unchanged."""
    from prestus_nexus.equipment import EquipmentLibrary, short_serial
    lib = EquipmentLibrary(prestus_root)
    assert lib.short_name("CTX_500_026") == "CTX-500 #026"
    assert lib.short_name("IS_PCD15287_01001") == "IS-PCD15287 #01001"
    assert lib.short_name("SC_203_035") == "SC-203 #035"
    assert lib.short_name("IGT_32_ch_comb_10_ch") == "IGT_32_ch_comb_10_ch"
    assert lib.short_name("NOT_IN_LIBRARY_7") == "NOT_IN_LIBRARY_7"
    assert lib.short_name("") == ""
    assert short_serial("CTX_250_014", "CTX_250_2ch") == "CTX-250 #014"
