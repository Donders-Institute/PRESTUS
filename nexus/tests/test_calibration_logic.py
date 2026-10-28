"""Calibration page logic without Qt: depth lists, study depths, meta.fit settings."""

import pytest

from prestus_nexus import calibration as cal
from prestus_nexus.equipment import EquipmentLibrary
from prestus_nexus.model import ConfigModel
from prestus_nexus.study import Manifest


@pytest.mark.parametrize("text, expected", [
    ("20:5:35", [20, 25, 30, 35]), ("30:32", [30, 31, 32]), ("42.5, 30 50", [30, 42.5, 50]),
    ("30;30.001", [30]), ("", None), ("a, b", None), ("40:5:20", None), ("1:2:3:4", None),
])
def test_parse_depths(text, expected):
    assert cal.parse_depths(text) == expected


def test_study_depths_from_config_and_manifest(schema, prestus_root):
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    m.set("transducer.serial", "CTX_250_001")
    m.set("transducer.combo.ds_serial", "SC_105_010")
    m.set("transducer.focal_distance_ep", 42.0)
    ti = m.add_transducer()
    m.set("transducer.serial", "CTX_500_006", ti)
    m.set("transducer.focal_distance_ep", 30.0, ti)
    man = Manifest()
    man.add_column("transducer.focal_distance_ep")
    man.add_column("transducer(2).focal_distance_ep")
    man.add_row({"subject_id": "1", "transducer.focal_distance_ep": "55",
                 "transducer(2).focal_distance_ep": "20"})
    man.add_row({"subject_id": "2", "transducer.focal_distance_ep": ""})
    assert cal.study_depths(m, "CTX_250_001", "SC_105_010", man) == [42.0, 55.0]
    assert cal.study_depths(m, "CTX_250_001", "", man) == []        # other driving system
    # transducer 2 inherits combo.ds_serial from transducer 1, as in PRESTUS
    assert cal.study_depths(m, "CTX_500_006", "SC_105_010", man) == [20.0, 30.0]


def test_target_depths_join_the_settings(schema, prestus_root):
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    m.set("transducer.serial", "CTX_250_001")
    m.set("transducer.focal_distance_ep", 42.0)
    ti = m.add_transducer()
    m.set("transducer.serial", "CTX_500_006", ti)
    planned = {0: 47.834, 1: 30.0}
    targets = cal.target_depths(m, "CTX_250_001", "", planned.get)
    assert targets == {0: 47.8}                                     # other device left out
    assert cal.target_depths(m, "CTX_250_001", "", lambda _ti: None) == {}
    settings = cal.study_depths(m, "CTX_250_001", "")
    assert cal.desired_depths(settings, targets) == [42.0, 47.8]
    assert cal.desired_depths([42.0], {0: 42.0}) == [42.0]
    assert cal.describe_depths(settings, targets) == "focal depth setting 42 mm; transducer 1 target 47.8 mm"
    assert cal.describe_depths([42.0], {0: 42.0}).endswith("(= setting)")


def test_settings_from_record_restores_the_fit():
    rec = {"tool": "calibration_gui", "mode": "per_depth", "calibrate_at": "desired",
           "depths_ep_mm": [35, 45], "isppa_wcm2": 20,
           "fit_options": {"FreeParams": ["phase", "c_steer"], "ForwardModel": "oneil",
                           "NearFieldSkipMm": 12, "HuberDelta": 0.5, "NStarts": 8, "Domain": "amplitude",
                           "Lambda": 0.001, "ZShiftBounds": [-3, 3], "PhaseConstraint": "linear",
                           "FocalWeight": 2},
           "profile_options": {"Align": "peak", "MaxExtrapolationMm": 4}}
    s, depths = cal.settings_from_record(rec, cal.FitSettings())
    assert (s.free_params, s.forward_model, s.domain) == (["phase"], "oneil", "amplitude")   # c_steer (removed static model) dropped
    assert (s.near_field_skip_mm, s.huber_delta, s.n_starts, s.isppa_wcm2) == (12, 0.5, 8, 20)
    assert (s.align, s.max_extrapolation_mm, s.calibrate_at) == ("peak", 4, "desired")
    assert s.advanced == {"PhaseConstraint": "linear", "FocalWeight": 2}       # defaults dropped
    assert depths == [35.0, 45.0]
    assert cal.settings_from_record(None, s) == (s, None)


def test_single_free_param_from_yaml_scalar():
    s, depths = cal.settings_from_record({"mode": "per_depth", "fit_options": {"FreeParams": "phase"}},
                                         cal.FitSettings())
    assert s.free_params == ["phase"] and depths is None


def test_ring_and_extra_options_round_trip():
    s = cal.FitSettings(advanced={"NearFieldPenalty": 0, "AmpBounds": [0.2, 5]},
                        extra={"Tolerance": 0.5, "GridDims": "auto", "ApplyCorrection": False},
                        ring_options={"NSub": 6, "OffAxis": "/data/planes.zip", "Kind": "transverse"})
    b = s.bridge_settings()
    assert b["tolerance"] == 0.5 and b["advanced"] == {"NearFieldPenalty": 0, "AmpBounds": [0.2, 5]}
    ro = b["ring_profile_options"]
    assert ro["NSub"] == 6 and ro["Smooth"] == 0.3 and ro["OffAxis"] == "/data/planes.zip"
    assert ro["OffAxisOptions"]["Kind"] == "transverse" and ro["OffAxisOptions"]["RMaxMm"] == 30
    assert "OffAxis" not in cal.ring_profile_options({})                 # axial only
    assert s.kwave_options() == {"GridDims": "auto", "ApplyCorrection": False}
    # meta.fit as calibration_fit writes it: ring options nested, k-Wave options kept from the page
    rec = {"mode": "per_depth", "fit_options": {"NearFieldPenalty": 0, "Seed": 1},
           "profile_options": {"Align": "label", "Tolerance": 0.5},
           "ring_profile": {"n_sub": 6, "options": ro}}
    r, _ = cal.settings_from_record(rec, s)
    assert r.advanced == {"NearFieldPenalty": 0}
    assert r.extra == {"GridDims": "auto", "ApplyCorrection": False, "Tolerance": 0.5}
    assert r.ring_options == {"NSub": 6, "OffAxis": "/data/planes.zip", "Kind": "transverse"}


def test_presets_match_the_advanced_keys():
    for name, p in cal.PRESETS.items():
        main = {"FreeParams", "ForwardModel", "Domain", "HuberDelta", "NearFieldSkipMm"}
        assert main <= set(p), name
        assert set(p) - main <= set(cal.ADVANCED_DEFAULTS), name
