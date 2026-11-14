"""checks.multitransducer_issues, heating_issues and elastic_issues (rules mirrored from
grid_transducer_location.m, source_create.m, combine_async_intensity.m, grid_axisymmetry.m,
thermal_parameters.m, thermal_grid_setup.m, thermal_stage.m, medium_elastic_properties.m)."""

from prestus_nexus import checks
from prestus_nexus.model import ConfigModel


def _two(schema, **values):
    m = ConfigModel(schema)
    m.add_transducer()
    for k, v in values.items():
        m.set(k.replace("__", "."), v)
    return m


def _paths(issues, level=None):
    return [(i.path, i.level) for i in issues if level is None or i.level == level]


def test_non_layered_media_keep_transducer_1(schema):
    m = _two(schema, simulation__medium="water")
    assert ("simulation.medium", "warning") in _paths(checks.multitransducer_issues(m))
    m.set("simulation.medium", "layered")
    assert not checks.multitransducer_issues(m)
    one = ConfigModel(schema)
    one.set("simulation.medium", "water")
    assert not checks.multitransducer_issues(one)


def test_axisymmetric_with_several_transducers_or_async_or_ras_plus(schema):
    m = _two(schema, grid__axisymmetric=1)
    assert ("grid.axisymmetric", "error") in _paths(checks.multitransducer_issues(m))
    m.set("simulation.transducer_coupling", "async")
    msgs = [i.message for i in checks.multitransducer_issues(m)]
    assert any("combine" in s for s in msgs)
    one = ConfigModel(schema)
    one.set("grid.axisymmetric", 1)
    assert not checks.multitransducer_issues(one)
    one.set("grid.mode", "ras_plus")
    assert ("grid.mode", "error") in _paths(checks.multitransducer_issues(one))


def test_heating_needs_timing(schema):
    m = ConfigModel(schema)
    assert not checks.heating_issues(m)                      # heating off by default
    m.set("modules.run_heating_sims", 1)
    errs = _paths(checks.heating_issues(m), "error")
    assert ("timing.pri", "error") in errs and ("timing.ptri", "error") in errs
    for k, v in {"timing.pri": 0.005, "timing.ptri": 1.0, "timing.pd": 0.0005,
                 "timing.ptd": 0.3, "timing.ptrd": 80}.items():
        m.set(k, v)
    assert not checks.heating_issues(m)
    m.set("timing.pd", 0)
    assert ("timing.pd", "warning") in _paths(checks.heating_issues(m))
    m.set("timing.pd", 0.0005)
    m.set("grid.thermal_fov_mm", [40, 40])
    assert ("grid.thermal_fov_mm", "error") in _paths(checks.heating_issues(m))


def test_thermal_fov_against_acoustic_fov(schema):
    m = ConfigModel(schema)
    m.set("grid.thermal_use_acoustic_fov", 1)
    assert ("grid.thermal_use_acoustic_fov", "info") in _paths(checks.heating_issues(m))
    m.set("grid.acoustic_fov_diameter_mm", 30)
    m.set("grid.acoustic_fov_length_mm", 90)
    assert not checks.heating_issues(m)
    m.set("grid.thermal_fov_mm", [50, 50, 60])
    assert ("grid.thermal_fov_mm", "info") in _paths(checks.heating_issues(m))
    m.set("grid.thermal_fov_mm", [20, 20, 60])
    assert not checks.heating_issues(m)


def test_elastic_shear_ppw(schema):
    m = ConfigModel(schema)
    m.set("acoustic_backend", "kwave_elastic")
    m.set("transducer.freq_hz", 500000)
    m.set("grid.resolution_mm", 1.0)
    issues = checks.elastic_issues(m)
    assert [i.path for i in issues] == ["grid.resolution_mm"]
    m.set("grid.resolution_mm", 0.2)
    assert not checks.elastic_issues(m)
    m.set("acoustic_backend", "kwave")
    m.set("grid.resolution_mm", 1.0)
    assert not checks.elastic_issues(m)
