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
    # every override is under hpc.* (the apply() guard assumes this)
    for p in ps:
        assert all(k == "hpc" or k.startswith("hpc.") for k in p.overrides)


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
