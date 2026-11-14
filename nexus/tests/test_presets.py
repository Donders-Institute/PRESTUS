"""run_presets.yaml / presets.py: named hpc.* bundles for the Run tab."""

import pytest

from prestus_nexus import presets
from prestus_nexus.model import ConfigModel


def test_shipped_presets_load():
    ps = presets.load_presets()
    names = [p.name for p in ps]
    assert names == ["Local", "Donders HPC — default", "Donders HPC — GPU, high-fidelity",
                     "Snellius — A100", "Snellius — H100"]
    assert all(p.description for p in ps)
    # every override is under hpc.* or platform (the apply() guard assumes this)
    for p in ps:
        assert p.overrides and all(k in ("hpc", "platform") or k.startswith("hpc.") for k in p.overrides)


def test_snellius_presets_use_the_exact_partition_selectors():
    # functions/hpc/hpc_submit_job.m's extract_snellius_parameters switches on this
    # literal string; a typo here would silently fall through to its 'otherwise' error.
    by_name = {p.name: p for p in presets.load_presets()}
    assert by_name["Snellius — A100"].overrides["hpc.partition"] == "gpu_a100"
    assert by_name["Snellius — H100"].overrides["hpc.partition"] == "gpu_h100"
    assert by_name["Snellius — A100"].overrides["hpc.name"] == "snellius"


def test_apply_sets_hpc_fields(schema):
    m = ConfigModel(schema)
    preset = presets.Preset("test", "test preset", {"hpc.name": "snellius", "hpc.partition": "gpu_a100"})
    presets.apply(m, preset)
    assert m.get("hpc.name") == "snellius" and m.get("hpc.partition") == "gpu_a100"


def test_apply_rejects_keys_outside_hpc(schema):
    m = ConfigModel(schema)
    preset = presets.Preset("bad", "", {"grid.resolution_mm": 1.0})
    with pytest.raises(ValueError, match="hpc"):
        presets.apply(m, preset)


def test_donders_default_resets_what_the_gpu_preset_sets(schema):
    # presets only write their own keys, so 'default' must clear the GPU preset's
    # partition/gres/limits or switching back would silently keep gpu40g.
    by_name = {p.name: p for p in presets.load_presets()}
    m = ConfigModel(schema)
    presets.apply(m, by_name["Donders HPC — GPU, high-fidelity"])
    presets.apply(m, by_name["Donders HPC — default"])
    assert m.get("hpc.partition") == "" and m.get("hpc.gpu") == ""
    assert m.get("hpc.timelimit") == "04:00:00" and m.get("hpc.memorylimit") == 20


def test_current_profile_follows_the_settings(schema):
    by_name = {p.name: p for p in presets.load_presets()}
    m = ConfigModel(schema)
    # the shipped hpc defaults are the Donders default profile
    assert presets.current(m, list(by_name.values())).name == "Donders HPC — default"
    presets.apply(m, by_name["Local"])
    assert m.get("platform") == "matlab"
    assert presets.current(m, list(by_name.values())).name == "Local"
    m.set("platform", "slurm")
    assert not presets.matches(m, by_name["Local"])
    m.set("hpc.memorylimit", 33)
    m.set("hpc.gpu", "gpu:1")
    assert presets.current(m, list(by_name.values())) is None        # matches no profile
