"""Uncertainty variants: which medium files are used and what each variant's properties are."""

import pytest

from prestus_nexus import checks
from prestus_nexus import uncertainty as unc
from prestus_nexus.model import ConfigModel
from prestus_nexus.yamlio import load_file


def _shipped(prestus_root, variant):
    return prestus_root / "config" / "uncertainty" / f"config_medium_{variant}.yaml"


def test_empty_keys_use_shipped_files(schema, prestus_root):
    m = ConfigModel(schema)
    for variant in unc.VARIANTS:
        assert unc.variant_file(m, prestus_root, variant) == _shipped(prestus_root, variant)
    sets = unc.property_sets(m, prestus_root)
    assert not sets.errors and set(sets.sets) == {"default", "conservative", "liberal"}


def test_file_name_and_absolute_path(schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("simulation.uncertainty_conservative", "my_bounds.yaml")
    m.set("simulation.uncertainty_liberal", str(tmp_path / "lib.yaml"))
    assert unc.variant_file(m, prestus_root, "conservative") == \
        prestus_root / "config" / "uncertainty" / "my_bounds.yaml"
    assert unc.variant_file(m, prestus_root, "liberal") == tmp_path / "lib.yaml"


def test_variant_values_merge_onto_default(schema, prestus_root):
    m = ConfigModel(schema)
    m.set("medium_properties.skin.sound_speed", 1650)          # an edited default
    sets = unc.property_sets(m, prestus_root)
    shipped = load_file(_shipped(prestus_root, "conservative"))["medium_properties"]
    assert sets.value("conservative", "brain", "sound_speed") == shipped["brain"]["sound_speed"]
    assert not sets.inherited("conservative", "brain", "sound_speed")
    # skin is not in the shipped files: every variant uses the (edited) default
    for variant in ("default", *unc.VARIANTS):
        assert sets.value(variant, "skin", "sound_speed") == 1650
    assert sets.inherited("liberal", "skin", "sound_speed")


def test_attenuation_at_transducer_frequency(schema, prestus_root):
    m = ConfigModel(schema)
    m.set("transducer.freq_hz", 250e3)
    sets = unc.property_sets(m, prestus_root)
    assert sets.freq_from_transducer and sets.freq_mhz == pytest.approx(0.25)
    brain = sets.sets["default"]["brain"]
    assert sets.value("default", "brain", "attenuation") == \
        pytest.approx(brain["alpha_coeff"] * 0.25 ** brain["alpha_power"])


def test_missing_or_invalid_file_is_an_issue(schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("simulation.uncertainty_liberal", str(tmp_path / "missing.yaml"))
    assert checks.uncertainty_issues(m) == []                     # uncertainty off: not used
    m.set("simulation.uncertainty", True)
    bad = tmp_path / "bad.yaml"
    bad.write_text("grid:\n  resolution_mm: 1\n")
    m.set("simulation.uncertainty_conservative", str(bad))
    found = {(i.path, i.level): i.message for i in checks.uncertainty_issues(m)}
    assert "file not found" in found[("simulation.uncertainty_liberal", "error")]
    assert "no medium_properties" in found[("simulation.uncertainty_conservative", "error")]
