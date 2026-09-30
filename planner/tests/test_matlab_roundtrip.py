"""Planner-written YAML must load in MATLAB exactly as the planner predicts.

Writes a study config with the planner, loads it with PRESTUS'
load_parameters in MATLAB, and compares every overridden value (and the
values transducer 2 inherits from transducer 1) with ConfigModel.effective().

Needs MATLAB: set PRESTUS_MATLAB to the matlab binary, or have it on PATH or
in /Applications. Skipped otherwise.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from prestus_planner.equipment import EquipmentLibrary
from prestus_planner.model import TX, ConfigModel, iter_leaves, values_equal


def _find_matlab() -> str | None:
    cand = [os.environ.get("PRESTUS_MATLAB"), shutil.which("matlab")]
    cand += sorted(glob.glob("/Applications/MATLAB_R20*.app/bin/matlab"), reverse=True)
    return next((c for c in cand if c and Path(c).exists()), None)


MATLAB = _find_matlab()
pytestmark = [pytest.mark.matlab,
              pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")]

# load_parameters rewrites or derives these; they are not a YAML round-trip question
# (elem_p_amp: validate_annular_transducer expands a scalar to one value per element)
DERIVED = {"io.output_affix", "path.seg", "transducer.focal_distance_offset",
           "transducer.focal_distance_bowl", "transducer.align_to_focus",
           "transducer.annular.elem_p_amp"}


def build_model(schema) -> ConfigModel:
    m = ConfigModel(schema)
    values = {
        "subject_id": 7,
        "platform": "matlab",
        "simulation.medium": "water",
        "simulation.code_type": "matlab_cpu",
        "path.anat": "/data/it's here",                 # quote inside string
        "path.sim": "/data/sims",
        "path.t1_pattern": "T1_%1$03d_run-*.nii.gz",    # MATLAB sprintf + glob
        "startup.modules_to_load": ["matlab/R2023b", "fsl/6.0.7"],
        "hpc.timelimit": "12:00:00",                    # sexagesimal if unquoted
        "hpc.gpu": "yes",                               # boolean if unquoted
        "hpc.job_prefix": "",                           # empty string over a non-empty default
        "hpc.wait_for_job": True,
        "hpc.memorylimit": 64,
        "io.save_MNI": 0,
        "grid.resolution_mm": 0.25,
        "grid.default_dims": [96, 96, 256],
        "grid.pml_size": 20,
        "timing.pd": 0.02,
        "timing.pri": 1e-1,
        "medium_properties.water.alpha_coeff": 2.2e-3,
        "placement.manual.coordinate_space": "ras",
        "transducer.name": "CTX-test",
        "transducer.type": "annular",
        "transducer.freq_hz": 5e5,
        "transducer.target_isppa_wcm2": [10, 20.5],
        "transducer.trans_pos": [-62.1, 12.4, 48.0],
        "transducer.focus_pos": [-21.5, 10.2, 2.3],
        "transducer.focal_distance_ep": 60,
        "transducer.annular.elem_n": 4,
        "transducer.annular.elem_id_mm": [0, 32.9, 46, 55.9],
        "transducer.annular.elem_od_mm": [32.3, 45.4, 55.4, 64],
        "transducer.annular.elem_p_amp": 100000,
        "transducer.annular.elem_phase_deg": [0, 12.5, 25, 37.5],
        "transducer.annular.curv_radius_mm": 63.2,
        "transducer.annular.dist_geom_ep_mm": 52.38,
    }
    for key, value in values.items():
        m.set(key, value)
    m.add_transducer()
    m.set("transducer.trans_pos", [60.0, 10.0, 45.5], 1)
    m.set("transducer.focus_pos", [20.0, 8.0, 2.0], 1)
    return m


def _matlab_key(path: str, ti: int | None) -> str:
    if ti is None:
        return path
    return f"{TX}({ti + 1})" + path[len(TX):]


def _normalise(expected, actual):
    """JSON from MATLAB: 1-element string arrays encode as scalars."""
    if isinstance(expected, list) and len(expected) == 1 and not isinstance(actual, list):
        return [actual]
    return actual


def test_load_parameters_reads_what_planner_writes(schema, prestus_root, tmp_path):
    model = build_model(schema)
    study = model.save(tmp_path / "config_roundtrip.yaml", header="round-trip test")

    effective = model.effective()
    txs = effective[TX]
    checks: dict[str, object] = {}
    for path, _ in iter_leaves(model.overrides):
        if path not in DERIVED:
            node = effective
            for part in path.split("."):
                node = node[part]
            checks[_matlab_key(path, None)] = node
    # every key set on transducer 1 must arrive on both transducers
    keys_t1 = [f"{TX}.{p}" for p, _ in iter_leaves(model.tx[0])]
    keys_t2 = [f"{TX}.{p}" for p, _ in iter_leaves(model.tx[1])]
    for ti, keys in ((0, keys_t1), (1, keys_t1 + keys_t2)):
        for path in keys:
            if path in DERIVED:
                continue
            node = txs[ti]
            for part in path.split(".")[1:]:
                node = node[part]
            checks[_matlab_key(path, ti)] = node

    keys_file = tmp_path / "keys.json"
    keys_file.write_text(json.dumps(sorted(checks)))
    out_file = tmp_path / "out.json"
    script = textwrap.dedent(f"""
        addpath(genpath('{prestus_root / "functions"}'));
        p = load_parameters('{study.name}', '{study.parent}');
        keys = jsondecode(fileread('{keys_file}'));
        res = containers.Map();
        for i = 1:numel(keys)
            res(keys{{i}}) = eval(['p.' keys{{i}}]);
        end
        fid = fopen('{out_file}', 'w'); fprintf(fid, '%s', jsonencode(res)); fclose(fid);
    """)
    (tmp_path / "roundtrip_check.m").write_text(script)
    proc = subprocess.run([MATLAB, "-batch", "roundtrip_check"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    got = json.loads(out_file.read_text())

    mismatches = {}
    for key, expected in checks.items():
        actual = _normalise(expected, got.get(key))
        if not values_equal(expected, actual):
            mismatches[key] = (expected, got.get(key))
    assert not mismatches, mismatches


def test_serial_transducer_resolves_like_the_planner(schema, prestus_root, tmp_path):
    """A library device: MATLAB's resolved geometry equals the planner's view."""
    model = ConfigModel(schema, EquipmentLibrary(prestus_root))
    for key, value in {"subject_id": 1, "path.sim": str(tmp_path), "simulation.medium": "water",
                       "transducer.serial": "CTX_500_026",
                       "transducer.combo.ds_serial": "SC_203_035",
                       "transducer.focal_distance_ep": 60}.items():
        model.set(key, value)
    study = model.save(tmp_path / "config_serial.yaml")
    tx = model.effective()[TX]
    keys = {"transducer.type": tx["type"], "transducer.freq_hz": tx["freq_hz"],
            **{f"transducer.annular.{k}": v for k, v in tx["annular"].items()
               if k in ("elem_n", "elem_id_mm", "elem_od_mm", "curv_radius_mm", "dist_geom_ep_mm")}}
    out = tmp_path / "out.json"
    script = textwrap.dedent(f"""
        addpath(genpath('{prestus_root / "functions"}'));
        p = load_parameters('{study.name}', '{study.parent}');
        keys = jsondecode(fileread('{tmp_path / "keys.json"}'));
        res = containers.Map();
        for i = 1:numel(keys)
            res(keys{{i}}) = eval(['p.' keys{{i}}]);
        end
        res('n_phases') = numel(p.transducer.annular.elem_phase_deg);
        fid = fopen('{out}', 'w'); fprintf(fid, '%s', jsonencode(res)); fclose(fid);
    """)
    (tmp_path / "keys.json").write_text(json.dumps(sorted(keys)))
    (tmp_path / "serial_check.m").write_text(script)
    proc = subprocess.run([MATLAB, "-batch", "serial_check"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    got = json.loads(out.read_text())
    mismatches = {k: (v, got.get(k)) for k, v in keys.items() if not values_equal(v, got.get(k))}
    assert not mismatches, mismatches
    assert got["n_phases"] == tx["annular"]["elem_n"]      # phases resolved from the library


def test_ras_placement_hits_the_intended_voxel(schema, prestus_root, tmp_path):
    """coordinate_space: ras -> PRESTUS placement stage -> the voxel the planner shows."""
    nib = pytest.importorskip("nibabel")
    import numpy as np

    from prestus_planner import placement as pl
    A = np.array([[-1.1, 0.05, 0, 90.3], [0.02, 1.2, 0.1, -126.7],
                  [0, -0.08, 0.9, -72.2], [0, 0, 0, 1]])
    data = np.zeros((40, 50, 60), np.float32)
    data[10, 20, 30] = 1                         # MATLAB element [11 21 31]
    data[20, 25, 40] = 2                         # MATLAB element [21 26 41]
    anat = tmp_path / "anat"
    anat.mkdir()
    nib.save(nib.Nifti1Image(data, A), str(anat / "sub-001_T1w.nii.gz"))
    vol = pl.Volume(anat / "sub-001_T1w.nii.gz")
    trans_ras = vol.to_ras([10, 20, 30])[0]
    focus_ras = vol.to_ras([20, 25, 40])[0]

    model = ConfigModel(schema)
    for key, value in {"subject_id": 1, "path.anat": str(anat), "path.sim": str(tmp_path),
                       "placement.manual.coordinate_space": "ras",
                       "transducer.trans_pos": [round(float(c), 3) for c in trans_ras],
                       "transducer.focus_pos": [round(float(c), 3) for c in focus_ras]}.items():
        model.set(key, value)
    study = model.save(tmp_path / "config_ras.yaml")
    out = tmp_path / "out.json"
    script = textwrap.dedent(f"""
        addpath(genpath('{prestus_root / "functions"}'));
        p = load_parameters('{study.name}', '{study.parent}');
        p = preproc_transducer_placement(p);
        img = niftiread(fullfile(p.path.anat, 'sub-001_T1w.nii.gz'));
        t = p.transducer.trans_pos; f = p.transducer.focus_pos;
        res = struct('trans', t, 'focus', f, 'v_trans', img(t(1), t(2), t(3)), ...
                     'v_focus', img(f(1), f(2), f(3)));
        fid = fopen('{out}', 'w'); fprintf(fid, '%s', jsonencode(res)); fclose(fid);
    """)
    (tmp_path / "ras_check.m").write_text(script)
    proc = subprocess.run([MATLAB, "-batch", "ras_check"], cwd=tmp_path,
                          capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    got = json.loads(out.read_text())
    assert got["trans"] == [11, 21, 31] == pl.ras_to_prestus_voxel(vol, trans_ras)
    assert got["focus"] == [21, 26, 41] == pl.ras_to_prestus_voxel(vol, focus_ras)
    assert (got["v_trans"], got["v_focus"]) == (1, 2)
