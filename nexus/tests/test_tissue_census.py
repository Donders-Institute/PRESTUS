"""Which tissues a segmentation gives voxels to (tissue_census.py mirrors check_layers.m)."""

import os

import numpy as np
import pytest

from prestus_nexus import tissue_census as tc

LAYERS = {"water": [0, 3, 6, 9, 10], "brain": [1, 2], "skin": [5], "skull": [4],
          "skull_cortical": [7], "skull_trabecular": [8]}
NAMES = ["water", "skin", "skull", "brain", "skull_trabecular", "skull_cortical"]
# a charm segmentation: no combined skull label 4
VOLUMES = {0: 8000.0, 1: 560.0, 2: 790.0, 3: 256.0, 5: 2500.0, 6: 20.0, 7: 560.0, 8: 113.0,
           9: 45.0, 10: 2.0}


def test_charm_segmentation_leaves_single_skull_unused():
    use = tc.tissue_use(LAYERS, NAMES, VOLUMES, pct=False)
    assert use["skull"].state == "absent" and not use["skull"].modelled
    assert use["skull_cortical"].modelled and use["skull_cortical"].volume_ml == 560.0
    assert use["brain"].volume_ml == 1350.0
    # water: labels no other layer claims, background excluded from the volume
    assert use["water"].state == "remainder"
    assert use["water"].labels == (3, 6, 9, 10)
    assert use["water"].volume_ml == pytest.approx(256 + 20 + 45 + 2)


def test_pct_merges_split_skull_into_single_layer():
    use = tc.tissue_use(LAYERS, NAMES, VOLUMES, pct=True)
    assert use["skull"].modelled and use["skull"].labels == (4, 7, 8)
    assert use["skull"].volume_ml == pytest.approx(673.0)
    assert use["skull_cortical"].state == use["skull_trabecular"].state == "merged"
    assert use["water"].labels == (3, 6, 9, 10)


def test_tissue_without_layer():
    layers = {k: v for k, v in LAYERS.items() if k != "skin"}
    use = tc.tissue_use(layers, NAMES, VOLUMES, pct=False)
    assert use["skin"].state == "no_layer"
    assert 5 in use["water"].labels                  # unclaimed skin label becomes water


def test_label_volumes_from_nifti(tmp_path):
    nib = pytest.importorskip("nibabel")
    seg = np.zeros((10, 10, 10), dtype=np.uint8)
    seg[:2] = 5
    seg[2:3] = 7
    f = tmp_path / "final_tissues.nii.gz"
    nib.save(nib.Nifti1Image(seg, np.diag([2.0, 2.0, 2.0, 1.0])), str(f))       # 8 µL voxels
    assert tc.cached_volumes(f) is None
    vols = tc.label_volumes(f)
    assert vols == {0: pytest.approx(700 * 0.008), 5: pytest.approx(200 * 0.008),
                    7: pytest.approx(100 * 0.008)}
    assert tc.cached_volumes(f) is vols
    os.utime(f, (1, 1))                              # changed file: counted again
    assert tc.cached_volumes(f) is None


def test_volume_text():
    assert tc.volume_text(2521.6) == "2.52 L"
    assert tc.volume_text(560.9) == "561 mL"
    assert tc.volume_text(0.4) == "0.4 mL"
    assert tc.volume_text(None) == ""
