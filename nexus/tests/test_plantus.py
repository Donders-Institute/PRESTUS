"""PlanTUS page logic: reading plantus_metrics.mat, colours, picking and the vertex readout."""

import numpy as np
import pytest
from scipy.io import savemat

from prestus_nexus import mesh3d as m3
from prestus_nexus import plantus as pt
from prestus_nexus.model import ConfigModel

TARGET = np.array([0.0, 10.0, 20.0])


def head():
    """A 'skin' sphere of radius 80 mm around the origin."""
    return m3.sphere([0.0, 0.0, 0.0], 80.0, n=20)


def metrics(verts, faces):
    """The variables position_transducer_plantus (plantus_common.save_metrics) writes."""
    n = len(verts)
    dist = np.linalg.norm(verts - TARGET, axis=1)
    avoid = np.ones(n)
    avoid[verts[:, 2] < -40] = 0                       # the bottom of the head is excluded
    score = (dist - dist.min()) / np.ptp(dist)
    score[avoid == 0] = np.inf
    best = int(np.argmin(score))
    return {
        "skin_coordinates": verts, "skin_normals": verts / 80.0, "skin_faces": faces.astype(np.int32),
        "target_center": TARGET,
        "skin_target_distances": dist, "skin_target_angles": np.full(n, 5.0),
        "skin_target_intersections": np.full(n, 3.0), "skin_skull_angles": np.full(n, 8.0),
        "skull_thickness": np.where(verts[:, 0] > 0, 6.0, np.nan), "avoidance": avoid,
        "score": score, "best_vertex": best, "selected_vertex": best, "selection_mode": "auto",
        "min_distance": 30.0, "max_distance": 90.0, "optimal_distance": 60.0,
        "weights": np.array([0.2] * 5),
        "weight_names": np.array(["skin_target_distances", "skin_target_angles",
                                  "skin_target_intersections", "skin_skull_angles",
                                  "skull_thickness"], dtype=object),
    }


def write_metrics(model: ConfigModel):
    path = pt.metrics_path(model)
    path.parent.mkdir(parents=True, exist_ok=True)
    savemat(path, metrics(*head()))
    (path.parent / pt.TX_CONFIG_FILENAME).write_text(
        "transducer_diameter: 64\nplane_offset: 11.5\nadditional_offset: 5\n", encoding="utf-8")
    return path


@pytest.fixture
def maps(tmp_path):
    path = tmp_path / pt.METRICS_FILENAME
    savemat(path, metrics(*head()))
    return pt.load(path)


def test_load_reads_maps_and_selection(maps):
    verts, faces = head()
    assert maps.n == len(verts) and maps.faces.shape == faces.shape
    assert np.allclose(maps.target, TARGET)
    assert maps.best == maps.selected and maps.selection_mode == "auto"
    assert set(maps.maps) == {k for k, *_ in pt.MAPS}
    assert maps.weights["skull_thickness"] == pytest.approx(0.2)
    assert (maps.optimal_distance, maps.min_distance, maps.max_distance) == (60.0, 30.0, 90.0)


def test_readout_rank_and_exclusion(maps):
    rows = dict(maps.readout(maps.best))
    assert rows["Composite score"].endswith(f"(rank 1 of {int(np.isfinite(maps.score).sum())})")
    excluded = int(np.flatnonzero(maps.avoidance == 0)[0])
    assert not maps.eligible(excluded) and maps.rank(excluded) is None
    assert "avoidance area" in dict(maps.readout(excluded))["Composite score"]
    left = int(np.flatnonzero(maps.coords[:, 0] < 0)[0])
    assert dict(maps.readout(left))["Skull thickness"] == "n/a"          # x < 0 side has no value
    assert "not on the skin" in maps.readout(maps.n + 3)[0][1]


def test_goodness_puts_the_best_end_bright(maps):
    g = pt.goodness(maps, "score")
    assert np.isnan(g[maps.avoidance == 0]).all()
    assert g[maps.best] == pytest.approx(1.0)
    d = pt.goodness(maps, "skin_target_distances")                      # closeness to 60 mm
    dev = np.abs(maps.maps["skin_target_distances"] - 60.0)
    ok = np.isfinite(d)
    assert d[ok][np.argmin(dev[ok])] == pytest.approx(d[ok].max())
    rgba = pt.colours(g)
    assert rgba.shape == (maps.n, 4)
    assert np.allclose(rgba[maps.avoidance == 0], pt.EXCLUDED_RGBA)


def test_pick_vertex_finds_the_vertex_under_the_ray(maps):
    v = maps.best
    p = maps.coords[v]
    origin = p * 3.0                                    # outside the head, on the same radial line
    assert pt.pick_vertex(origin, -p, maps.coords, maps.faces) == v
    assert pt.pick_vertex(origin, p, maps.coords, maps.faces) is None   # pointing away: no hit
    # the first surface along the ray, not the far side
    hit = pt.ray_hit(np.array([0.0, 0.0, 300.0]), np.array([0.0, 0.0, -1.0]), maps.coords, maps.faces)
    assert hit is not None and hit[1][2] == pytest.approx(80.0, abs=1.0)


def test_preview_puts_the_transducer_behind_the_skin(maps):
    tx = {"transducer_diameter": 64, "plane_offset": 11.5, "additional_offset": 5}
    p = pt.preview(maps, maps.best, tx)
    skin = maps.coords[maps.best]
    assert np.allclose(p["axis"], (TARGET - skin) / np.linalg.norm(TARGET - skin))
    assert np.linalg.norm(p["centre"] - skin) == pytest.approx(16.5)
    assert np.dot(p["centre"] - skin, p["axis"]) < 0                    # outside the head
    assert p["diameter"] == 64


def test_metrics_path_is_the_placement_run_folder(schema, tmp_path):
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path / "sim"))
    m.set("subject_id", 7)
    m.set("placement.plantus.target_name", "M1")
    assert pt.metrics_path(m) == tmp_path / "sim" / "placement" / "sub-007" / "plantus" / "M1" / pt.METRICS_FILENAME
    write_metrics(m)
    assert pt.read_tx_config(pt.metrics_path(m).parent)["plane_offset"] == 11.5
    assert pt.read_tx_config(tmp_path) == {}
