"""A custom annular transducer is driven by its focal depth setting and target ISPPA;
element phases and amplitudes are optional (fill_annular_drive.m)."""

from prestus_nexus import checks
from prestus_nexus.model import ConfigModel


def _custom(schema, **extra):
    m = ConfigModel(schema)
    for k, v in {"transducer.type": "annular", "transducer.freq_hz": 500e3,
                 "transducer.annular.elem_n": 2, "transducer.annular.elem_id_mm": [0.0, 30.0],
                 "transducer.annular.elem_od_mm": [29.0, 42.0],
                 "transducer.annular.curv_radius_mm": 62.9, **extra}.items():
        m.set(k, v)
    return m


def _drive(m):
    return {(i.level, i.path) for i in checks.drive_issues(m)}


def test_depth_and_isppa_are_enough(schema):
    m = _custom(schema, **{"transducer.focal_distance_ep": 50, "transducer.target_isppa_wcm2": 20})
    assert _drive(m) == set()


def test_without_depth_or_phases_there_is_no_steering(schema):
    m = _custom(schema, **{"transducer.target_isppa_wcm2": 20})
    assert ("error", "transducer.focal_distance_ep") in _drive(m)
    m.set("transducer.annular.elem_phase_deg", [0.0, 40.0])          # calibrated phases instead
    assert ("error", "transducer.focal_distance_ep") not in _drive(m)


def test_without_isppa_or_amplitude_there_is_no_drive_level(schema):
    m = _custom(schema, **{"transducer.focal_distance_ep": 50, "transducer.target_isppa_wcm2": []})
    assert ("error", "transducer.target_isppa_wcm2") in _drive(m)
    m.set("transducer.annular.elem_p_amp", [1e5, 1e5])
    assert _drive(m) == set()


def test_library_devices_are_not_checked_here(schema):
    m = _custom(schema, **{"transducer.serial": "CTX_500_026"})
    assert _drive(m) == set()


def test_phases_and_amplitudes_are_advanced_drive_fields(schema):
    for key in ("transducer.annular.elem_phase_deg", "transducer.annular.elem_p_amp"):
        node = schema.get(key)
        assert node.advanced and node.group == "tx_drive"
    assert schema.get("transducer.focal_distance_ep").group == "tx_drive"
    assert not schema.get("transducer.focal_distance_ep").advanced
