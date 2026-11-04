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
    assert not sets.errors and set(sets.sets) == {"default", "temp_high", "temp_low", "trans_high", "trans_low"}


def test_file_name_and_absolute_path(schema, prestus_root, tmp_path):
    m = ConfigModel(schema)
    m.set("simulation.uncertainty_temp_low", "my_bounds.yaml")
    m.set("simulation.uncertainty_trans_high", str(tmp_path / "th.yaml"))
    assert unc.variant_file(m, prestus_root, "temp_low") == \
        prestus_root / "config" / "uncertainty" / "my_bounds.yaml"
    assert unc.variant_file(m, prestus_root, "trans_high") == tmp_path / "th.yaml"


def test_variant_values_merge_onto_default(schema, prestus_root):
    m = ConfigModel(schema)
    m.set("medium_properties.skin.sound_speed", 1650)          # an edited default
    sets = unc.property_sets(m, prestus_root)
    shipped = load_file(_shipped(prestus_root, "temp_low"))["medium_properties"]
    assert sets.value("temp_low", "brain", "sound_speed") == shipped["brain"]["sound_speed"]
    assert not sets.inherited("temp_low", "brain", "sound_speed")
    # skin is not in the shipped files: every variant uses the (edited) default
    for variant in ("default", *unc.VARIANTS):
        assert sets.value(variant, "skin", "sound_speed") == 1650
    assert sets.inherited("temp_high", "skin", "sound_speed")


def test_trans_files_are_acoustic_only(prestus_root):
    # the trans variants run without heating: their shipped files set no thermal property
    for variant in unc.PAIRS["trans"]:
        listed = load_file(_shipped(prestus_root, variant))["medium_properties"]
        for tissue, props in listed.items():
            assert not set(props) & set(unc.THERMAL), (variant, tissue)


def test_pairs_bound_in_opposite_directions(schema, prestus_root):
    # brain sound speed raises heating and intensity in opposite directions (uncertainty_variants.m)
    sets = unc.property_sets(ConfigModel(schema), prestus_root)
    assert sets.value("temp_high", "brain", "sound_speed") == sets.value("trans_low", "brain", "sound_speed")
    assert sets.value("temp_low", "brain", "sound_speed") == sets.value("trans_high", "brain", "sound_speed")


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
    m.set("simulation.uncertainty_trans_low", str(tmp_path / "missing.yaml"))
    assert checks.uncertainty_issues(m) == []                     # uncertainty off: not used
    m.set("simulation.uncertainty", True)
    bad = tmp_path / "bad.yaml"
    bad.write_text("grid:\n  resolution_mm: 1\n")
    m.set("simulation.uncertainty_temp_high", str(bad))
    found = {(i.path, i.level): i.message for i in checks.uncertainty_issues(m)}
    assert "file not found" in found[("simulation.uncertainty_trans_low", "error")]
    assert found[("simulation.uncertainty_trans_low", "error")].startswith("Trans-L variant")
    assert "no medium_properties" in found[("simulation.uncertainty_temp_high", "error")]


def test_pct_band_and_fraction_outside(schema, prestus_root):
    pytest.importorskip("PyQt6")
    from prestus_nexus.ui.pct_page import fraction_outside, skull_band
    sets = unc.property_sets(ConfigModel(schema), prestus_root)
    bounds = {v: sets.listed[v] for v in unc.VARIANTS}
    lo, hi = skull_band(bounds, "sound_speed")
    assert (lo, hi) == (min(sets.value(v, "skull", "sound_speed") for v in unc.VARIANTS),
                        max(sets.value(v, "skull", "sound_speed") for v in unc.VARIANTS))
    # bin centres 5, 15, 25, 35: the first and last lie outside [10, 30]
    assert fraction_outside([0, 10, 20, 30, 40], [1, 2, 3, 4], 10, 30) == pytest.approx(0.5)
    assert fraction_outside([0, 1], [0], 0, 1) is None


def test_absorption_field_replaces_inherited_one(schema, prestus_root):
    # medium_absorption_layer.m: the Temp files set skull absorption_coeff, which clears the default fraction
    m = ConfigModel(schema)
    m.set("transducer.freq_hz", 500e3)
    sets = unc.property_sets(m, prestus_root)
    listed = sets.listed["temp_high"]["skull"]
    assert sets.value("temp_high", "skull", "absorption_fraction") is None
    assert sets.value("temp_high", "skull", "absorption") == \
        pytest.approx(listed["absorption_coeff"] * 0.5 ** listed["alpha_power"])
    brain = sets.sets["temp_low"]["brain"]
    assert sets.value("temp_low", "brain", "absorption") == pytest.approx(
        brain["absorption_fraction"] * brain["alpha_coeff"] * 0.5 ** brain["alpha_power"])
    assert sets.inherited("trans_high", "skull", "absorption")


def _stats_json(path, predicted=False, source="native", n=10):
    import json
    data = {"n_voxels": n, "properties": {}}
    if predicted:
        data.update(predicted=True, source=source)
    path.write_text(json.dumps(data))
    return path


def test_pct_stats_file_prefers_the_default_run(tmp_path):
    pytest.importorskip("PyQt6")
    import os
    from prestus_nexus.ui.pct_page import stats_candidates
    cache = tmp_path / "cache"
    assert stats_candidates(cache, 7, "layered", "") == (None, None)       # no cache folder
    cache.mkdir()
    pred = _stats_json(cache / "sub-007_layered_pct_prediction_desc-a.json", predicted=True)
    assert stats_candidates(cache, 7, "layered", "_desc-a") == (None, pred)
    assert stats_candidates(cache, 7, "layered", "") == (None, None)       # prediction matched exactly
    other = _stats_json(cache / "sub-007_layered_medium_stats_desc-b.json")
    assert stats_candidates(cache, 7, "layered", "_desc-a") == (other, pred)   # newest real stats
    own = _stats_json(cache / "sub-007_layered_medium_stats_desc-a.json")
    os.utime(own, (1, 1))                                                   # older, but the default's own
    assert stats_candidates(cache, 7, "layered", "_desc-a") == (own, pred)


def test_pct_stats_title_says_which_it_is(tmp_path):
    pytest.importorskip("PyQt6")
    from prestus_nexus.ui.pct_page import stats_title
    real = tmp_path / "sub-001_layered_medium_stats.json"
    assert "mapped skull of the default run" in stats_title({"n_voxels": 5}, real)
    assert "n=5 skull voxels" in stats_title({"n_voxels": 5}, real)
    native = stats_title({"predicted": True, "source": "native", "n_voxels": 3}, real)
    assert "predicted from the pseudo-CT (native voxels: the pCT image before cropping" in native
    grid = stats_title({"predicted": True, "source": "grid"}, real)
    assert "grid voxels: the exact preprocessed grid" in grid and "—" not in grid + native


def test_pct_constant_mapping_draws_one_bar_and_counts_exactly():
    pytest.importorskip("PyQt6")
    from prestus_nexus.ui.pct_page import histogram_bars, outside_share
    # jsonencode writes a one-bin histogram as scalars; histcounts of a constant spans a unit bin
    prop = {"min": 4.7, "max": 4.7, "histogram": {"counts": 120, "edges": [4.5, 5.5]}}
    assert histogram_bars(prop) == ([4.7], [120], [0.0])
    assert outside_share(prop, 4.0, 4.8) == 0.0          # the value, not the bin centre (5.0), decides
    assert outside_share(prop, 1.0, 4.0) == 1.0
    spread = {"min": 0, "max": 40, "histogram": {"counts": [1, 2, 3, 4], "edges": [0, 10, 20, 30, 40]}}
    assert histogram_bars(spread) == ([5, 15, 25, 35], [1, 2, 3, 4], [10, 10, 10, 10])
    assert outside_share(spread, 10, 30) == pytest.approx(0.5)
    assert histogram_bars({"histogram": {"counts": [], "edges": []}}) is None


def test_bridge_allows_nexus_pct():
    from prestus_nexus.bridge import ALLOWED
    assert "nexus_pct" in ALLOWED
