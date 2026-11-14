"""
plantus_common.py  –  shared PlanTUS pipeline for both PRESTUS launchers
========================================================================
prestus_plantus_launcher.py (headless) and prestus_plantus_gui_launcher.py
(interactive) differ only in how the skin vertex is selected. Everything
else — config parsing, per-vertex metrics, composite scoring, the metric
maps and scenes written to disk, and plantus_metrics.mat — lives here so
both modes produce the same outputs.
"""

import os
import math

import numpy as np
import scipy.io

METRICS_FILENAME = "plantus_metrics.mat"
# Written by plantus_evaluation.py from the metrics file
EVALUATION_FILENAME = "plantus_evaluation.json"
SWEEP_FILENAME = "plantus_weight_sensitivity.csv"
MAPS_FILENAME = "plantus_criteria_maps.png"
SENSITIVITY_FILENAME = "plantus_weight_sensitivity.png"

# CHARM tissue labels (final_tissues.nii.gz) used by the skull-thickness metric
TISSUE_LABELS_FILENAME = "final_tissues.nii.gz"
BONE_LABELS = (4, 7, 8)          # bone, compact bone, spongy bone
INTRACRANIAL_LABELS = (1, 2, 3)  # WM, GM, CSF
SKULL_THICKNESS_DEPTH_MM = 50.0  # ray length from the skin
SKULL_THICKNESS_STEP_MM = 0.25   # sampling step along the ray


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def ensure_list(v):
    """Return v as a list; wraps scalars that YAML parsed as a single number."""
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return list(v)
    return [v]


def load_yaml(path):
    try:
        import yaml
        with open(path, "r") as fh:
            return yaml.safe_load(fh)
    except ImportError:
        pass
    # Minimal hand-rolled parser for flat key: value YAML
    cfg = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                continue
            k, _, v = line.partition(":")
            k = k.strip()
            v = v.strip()
            try:
                if v.startswith("["):
                    import ast
                    cfg[k] = ast.literal_eval(v)
                elif "." in v:
                    cfg[k] = float(v)
                elif v.lower() in ("true", "false"):
                    cfg[k] = v.lower() == "true"
                else:
                    cfg[k] = int(v)
            except (ValueError, SyntaxError):
                cfg[k] = v
    return cfg


def _as_bool(v):
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes")
    return bool(v)


def read_config(cfg):
    """Return the placement settings from a PlanTUSTxConfig.yaml dict."""
    s = {
        "max_distance":        float(cfg.get("max_distance",        70.0)),
        "min_distance":        float(cfg.get("min_distance",        30.0)),
        "transducer_diameter": float(cfg.get("transducer_diameter", 60.0)),
        "max_angle":           float(cfg.get("max_angle",           10.0)),
        "plane_offset":        float(cfg.get("plane_offset",         0.0)),
        "additional_offset":   float(cfg.get("additional_offset",    0.0)),
        "focal_distance_list": ensure_list(cfg.get("focal_distance_list", [])),
        "flhm_list":           ensure_list(cfg.get("flhm_list",           [])),
        "weights": [float(cfg.get("weight_skin_target_distances",     0.2)),
                    float(cfg.get("weight_skin_target_angles",        0.2)),
                    float(cfg.get("weight_skin_target_intersections", 0.2)),
                    float(cfg.get("weight_skin_skull_angles",         0.2)),
                    float(cfg.get("weight_skull_thickness",           0.2))],
        "save_metrics":        _as_bool(cfg.get("save_metrics", True)),
        "evaluate":            _as_bool(cfg.get("evaluate", True)),
        "connectome_wb_path":  str(cfg.get("connectome_wb_path", "") or "").strip(),
        # PlanTUS >= 2.0 only (Taubin passes on the skin surface; 0 = raw mesh)
        "skin_smoothing_iterations": int(cfg.get("skin_smoothing_iterations", 2) or 0),
    }
    s["optimal_distance"] = float(cfg.get(
        "optimal_distance", (s["min_distance"] + s["max_distance"]) / 2))
    return s


def resolve_plantus_root(cli_root, launcher_file):
    """Return the PlanTUS installation directory (contains code/PlanTUS.py).

    Priority: --plantus_root > sibling of PRESTUS > <project>/tools/PlanTUS.
    Returns (root_or_None, candidates_tried).
    """
    here = os.path.dirname(os.path.abspath(launcher_file))
    candidates = []
    if cli_root:
        candidates.append(cli_root)
    candidates.append(os.path.join(here, "..", "..", "..", "PlanTUS"))
    candidates.append(os.path.join(here, "..", "..", "..", "..", "tools", "PlanTUS"))
    for cand in candidates:
        cand = os.path.normpath(cand)
        if os.path.isfile(os.path.join(cand, "code", "PlanTUS.py")):
            return cand, candidates
    return None, candidates


def plantus_major_version(plantus_root):
    """2 for PlanTUS >= 2.0 (no Connectome Workbench), else 1.

    Read from the source rather than by importing it, so the dependency check
    can run before the import (which itself fails on a missing module).
    """
    path = os.path.join(plantus_root, "code", "PlanTUS.py")
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            src = fh.read()
    except OSError:
        return 1
    return 2 if "def convert_simnibs_mesh_to_surfaces(" in src else 1


def check_dependencies(version=1, need_wb_view=False):
    """Return a list of missing PlanTUS runtime dependencies (empty if none).

    PlanTUS 1.x shells out to wb_command via os.system and ignores failures,
    so a missing Workbench only surfaces later as unrelated fslmaths/GIFTI
    errors; nilearn is imported lazily deep inside the metric computation.
    Checking up front turns those into one actionable message. PlanTUS 2.x
    needs no Workbench; need_wb_view is set by the Workbench picker.
    """
    import shutil
    import importlib.util
    missing = []
    tools = (["wb_command"] if version < 2 else []) + (["wb_view"] if need_wb_view else [])
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(f"{tool} (Connectome Workbench): set "
                           "placement.plantus.connectome_wb_path to its bin/ "
                           "directory, or install connectome-workbench into the "
                           "Python env")
    modules = ("nilearn", "vtk", "h5py") + (("trimesh", "pandas") if version >= 2 else ())
    for mod in modules:
        if importlib.util.find_spec(mod) is None:
            missing.append(f"Python module '{mod}': pip install {mod} into "
                           "the env at placement.plantus.env_path")
    return missing


def stage_mask(mask_path, output_path):
    """Return the mask path inside output_path, copying it there if needed.

    PlanTUS expects the mask in its output directory (PlanTUS_wrapper.py
    convention); PRESTUS already writes it there.
    """
    import shutil
    dst = os.path.join(output_path, os.path.basename(mask_path))
    if not os.path.exists(dst) or not os.path.samefile(mask_path, dst):
        shutil.copy(mask_path, dst)
    return dst


def scene_templates(plantus_root):
    """Return (planning, placement) scene template paths; '' if missing."""
    d = os.path.join(plantus_root, "resources", "scene_templates")
    out = []
    for name in ("TUSTransducerPlacementPlanning_TEMPLATE.scene",
                 "TUSTransducerPlacement_TEMPLATE.scene"):
        p = os.path.join(d, name)
        out.append(p if os.path.isfile(p) else "")
    return tuple(out)


# ---------------------------------------------------------------------------
# Composite scoring
# ---------------------------------------------------------------------------

def _safe_normalise(arr, invert=False):
    """Normalise arr to [0,1]; NaN entries become 1 (worst score)."""
    out = np.array(arr, dtype=float)
    valid = np.isfinite(out)
    if valid.sum() == 0:
        return np.ones_like(out)
    mn, mx = out[valid].min(), out[valid].max()
    if mx > mn:
        out = (out - mn) / (mx - mn)
    else:
        out[:] = 0.0
    if invert:
        out = 1.0 - out
    out[~valid] = 1.0
    return out


def compute_composite_score(distances, angles, intersections,
                             skull_angles, avoidance,
                             w_dist, w_ang, w_inter, w_skull, w_thick,
                             max_distance, min_distance=0.0,
                             optimal_distance=None, skull_thickness=None):
    """Return a per-vertex composite cost (lower = better).

    Vertices with avoidance==0, distance<min_distance, or distance>max_distance
    get score=inf.  Thicker skull is worse; vertices where no bone was found
    (NaN) get the worst thickness sub-score.  When skull_thickness is None or
    all-NaN (metrics files written before the criterion existed, or no label
    volume) w_thick is ignored and the other four weights are renormalised.

    The distance sub-score penalises deviation from optimal_distance (defaults
    to max_distance when not provided) rather than raw distance, so that
    vertices where the acoustic focus lands near the target are preferred.
    """
    if optimal_distance is None:
        optimal_distance = max_distance

    # Distance score: deviation from optimal depth (lower deviation = better)
    dist_dev = np.abs(np.array(distances, dtype=float) - optimal_distance)
    s_dist  = _safe_normalise(dist_dev)           # smaller deviation → lower score
    s_ang   = _safe_normalise(angles)
    s_inter = _safe_normalise(intersections, invert=True)
    s_skull = _safe_normalise(skull_angles)
    if skull_thickness is None or not np.isfinite(skull_thickness).any():
        w_thick = 0.0
        s_thick = 0.0
    else:
        s_thick = _safe_normalise(skull_thickness)
    total_w = w_dist + w_ang + w_inter + w_skull + w_thick
    if total_w <= 0:
        total_w = 1.0
    score = (w_dist  * s_dist +
             w_ang   * s_ang  +
             w_inter * s_inter +
             w_skull * s_skull +
             w_thick * s_thick) / total_w

    # Mask out unavailable / out-of-range vertices
    score[np.array(avoidance) == 0] = np.inf
    score[np.array(distances) < min_distance] = np.inf
    score[np.array(distances) > max_distance] = np.inf

    return score


def score_metrics(m, s):
    """Composite score for the metrics dict m under settings s.

    Returns (score, best_vertex); best_vertex is None if no vertex is eligible.
    """
    score = compute_composite_score(
        m["skin_target_distances"], m["skin_target_angles"],
        m["skin_target_intersections"], m["skin_skull_angles"],
        m["avoidance"], *s["weights"],
        s["max_distance"], min_distance=s["min_distance"],
        optimal_distance=s["optimal_distance"],
        skull_thickness=m.get("skull_thickness"))
    if not np.isfinite(score).any():
        return score, None
    return score, int(np.argmin(score))


def vertex_summary(m, s, score, vertex):
    """Lines describing one skin vertex: its criteria, eligibility and rank.

    Shown after each click in wb_view and in the review dialog.
    """
    n = len(score)
    if not 0 <= vertex < n:
        return [f"Vertex {vertex} is not on the skin surface ({n} vertices)."]

    def fmt(v, unit):
        return "n/a" if not np.isfinite(v) else f"{v:.1f} {unit}"

    d = float(m["skin_target_distances"][vertex])
    lines = [
        f"Skin-target distance  {fmt(d, 'mm')} "
        f"(optimal {s['optimal_distance']:.1f}, window "
        f"{s['min_distance']:.1f}-{s['max_distance']:.1f} mm)",
        f"Beam angle at skin    {fmt(float(m['skin_target_angles'][vertex]), 'deg')}",
        f"Beam-target overlap   {fmt(float(m['skin_target_intersections'][vertex]), 'mm')}",
        f"Skin-skull angle      {fmt(float(m['skin_skull_angles'][vertex]), 'deg')}",
        f"Skull thickness       {fmt(float(m['skull_thickness'][vertex]), 'mm')}",
    ]
    n_eligible = int(np.isfinite(score).sum())
    if np.isfinite(score[vertex]):
        rank = int((score < score[vertex]).sum()) + 1
        lines.append(f"Composite score       {score[vertex]:.4f} "
                     f"(rank {rank} of {n_eligible} eligible)")
    else:
        reasons = []
        if m["avoidance"][vertex] == 0:
            reasons.append("in the avoidance mask (eyes, ears, air cavities)")
        if d < s["min_distance"] or d > s["max_distance"]:
            reasons.append("outside the distance window")
        lines.append("Not eligible under the composite score"
                     + (": " + "; ".join(reasons) if reasons else "") + ".")
    return lines


# ---------------------------------------------------------------------------
# Per-vertex metrics
# ---------------------------------------------------------------------------

def _write_metric(PlanTUS, name, values, output_path):
    """Write a skin metric map (masked by the avoidance mask) as .func.gii."""
    PlanTUS.create_metric_from_pseudo_nifti(
        name, values, output_path + "/skin.surf.gii")
    PlanTUS.mask_metric(output_path + "/" + name + "_skin.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/" + name + "_skin.func.gii", "CORTEX_LEFT")


def _intersection_lengths(intersections):
    vals = []
    for segs in intersections:
        if len(segs) <= 1:
            vals.append(0)
        elif len(segs) in (2, 3):
            vals.append(np.linalg.norm(np.array(segs[1]) - np.array(segs[0])))
        elif len(segs) == 4:
            vals.append(np.linalg.norm(np.array(segs[1]) - np.array(segs[0])) +
                        np.linalg.norm(np.array(segs[3]) - np.array(segs[2])))
        else:
            vals.append(np.nan)
    return np.array(vals)


def skull_thickness_along_normals(points, normals, labels_path,
                                  max_depth=SKULL_THICKNESS_DEPTH_MM,
                                  step=SKULL_THICKNESS_STEP_MM, chunk=2048):
    """Bone path length [mm] along the inward skin normal of each vertex.

    Samples the CHARM label volume (nearest voxel) from the skin point inwards
    and counts bone samples before the first intracranial sample (WM, GM,
    CSF).  This follows the same ray as the skin-skull angle.  NaN where the
    ray does not reach intracranial tissue within max_depth (face, neck, rays
    grazing the skull base, which would otherwise score as very thin skull)
    or meets no bone on the way.  The skull surface mesh is not used: it
    lacks the inner skull boundary, so entry/exit pairs cannot be read from
    it.
    """
    import nibabel as nib

    img = nib.load(labels_path)
    labels = np.asarray(img.dataobj).squeeze().astype(np.int32)
    vox_from_ras = np.linalg.inv(img.affine)
    depths = np.arange(0.0, max_depth + step / 2, step)
    shape = np.array(labels.shape)
    bone_lut = np.isin(np.arange(labels.max() + 1), BONE_LABELS)
    icr_lut = np.isin(np.arange(labels.max() + 1), INTRACRANIAL_LABELS)

    points = np.asarray(points, dtype=float)
    normals = np.asarray(normals, dtype=float)
    thickness = np.full(len(points), np.nan)
    for start in range(0, len(points), chunk):
        p = points[start:start + chunk]
        n = normals[start:start + chunk]
        # normals point away from the head, so step along -n
        ras = p[:, None, :] - depths[None, :, None] * n[:, None, :]
        ijk = np.rint(ras @ vox_from_ras[:3, :3].T + vox_from_ras[:3, 3]).astype(int)
        inside = np.all((ijk >= 0) & (ijk < shape), axis=2)
        lab = np.zeros(inside.shape, dtype=np.int32)
        lab[inside] = labels[ijk[inside, 0], ijk[inside, 1], ijk[inside, 2]]
        bone = bone_lut[lab]
        icr = icr_lut[lab]
        # bone samples before the first intracranial one
        reached = icr.any(axis=1)
        first_icr = icr.argmax(axis=1)
        before = np.arange(len(depths))[None, :] < first_icr[:, None]
        n_bone = (bone & before).sum(axis=1)
        thickness[start:start + chunk] = np.where(
            reached & (n_bone > 0), n_bone * step, np.nan)
    return thickness


def is_v2(PlanTUS):
    """PlanTUS >= 2.0: Workbench-free, renamed mesh conversion, own viewer."""
    return hasattr(PlanTUS, "convert_simnibs_mesh_to_surfaces")


def mesh_to_surface(PlanTUS, mesh_path, tags, name, output_path, smoothing=0):
    """<name>.surf.gii and <name>.stl from tagged SimNIBS mesh elements.

    PlanTUS 2.x renamed convert_simnibs_mesh_to_surface(s) and added skin
    smoothing; 1.x has no smoothing, so it is ignored there.
    """
    if is_v2(PlanTUS):
        PlanTUS.convert_simnibs_mesh_to_surfaces(
            mesh_path, tags, name, output_path, smoothing_iterations=int(smoothing))
    else:
        PlanTUS.convert_simnibs_mesh_to_surface(mesh_path, tags, name, output_path)


def compute_metrics(PlanTUS, mesh_path, mask_path, output_path, s, log):
    """Run PlanTUS steps 2–4 and write the metric maps into output_path.

    mask_path must already be inside output_path. Returns a dict of per-vertex
    arrays (skin coordinates/normals/faces, target centre, the five criteria
    and the avoidance mask).  Skull thickness is read from the CHARM label
    volume next to the mesh; without it the criterion is all-NaN, its weight
    is ignored in scoring, and a warning is logged.
    """
    import nibabel as nib

    log("Converting SimNIBS mesh to surfaces...")
    mesh_to_surface(PlanTUS, mesh_path, [1005], "skin", output_path,
                    s.get("skin_smoothing_iterations", 0))
    PlanTUS.add_structure_information(output_path + "/skin.surf.gii",  "CORTEX_LEFT")
    mesh_to_surface(PlanTUS, mesh_path, [1007, 1008], "skull", output_path)
    PlanTUS.add_structure_information(output_path + "/skull.surf.gii", "CORTEX_RIGHT")

    log("Computing avoidance mask...")
    PlanTUS.create_avoidance_mask(
        mesh_path, output_path + "/skin.surf.gii", s["transducer_diameter"] / 2)

    log("Computing skin-to-target metrics...")
    target_center = PlanTUS.roi_center_of_gravity(mask_path)

    distances = PlanTUS.distance_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)
    _write_metric(PlanTUS, "distances", distances, output_path)
    PlanTUS.threshold_metric(output_path + "/distances_skin.func.gii", s["max_distance"])
    PlanTUS.mask_metric(output_path + "/distances_skin_thresholded.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/distances_skin_thresholded.func.gii", "CORTEX_LEFT")

    skin_coordinates, skin_normals = PlanTUS.compute_surface_metrics(
        output_path + "/skin.surf.gii")
    vectors = PlanTUS.vectors_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)
    angles = np.abs(np.array([
        math.degrees(PlanTUS.angle_between_vectors(vectors[i], skin_normals[i]))
        for i in np.arange(len(vectors))]))
    _write_metric(PlanTUS, "angles", angles, output_path)

    log("Computing skin-target intersection metric...")
    mask_name = os.path.basename(mask_path).replace(".nii.gz", "").replace(".nii", "")
    PlanTUS.stl_from_nii(mask_path, 0.25)
    intersections = _intersection_lengths(PlanTUS.compute_vector_mesh_intersections(
        skin_coordinates, skin_normals,
        output_path + "/" + mask_name + "_3Dmodel.stl", 200))
    _write_metric(PlanTUS, "target_intersection", intersections, output_path)

    log("Computing skin-skull angle metric...")
    skull_coordinates, skull_normals = PlanTUS.compute_surface_metrics(
        output_path + "/skull.surf.gii")
    skin_skull_intersections = PlanTUS.compute_vector_mesh_intersections(
        skin_coordinates, skin_normals, output_path + "/skull.stl", 40)
    skull_angle_list = []
    for i in np.arange(len(skin_coordinates)):
        try:
            intersection_coordinate = skin_skull_intersections[i][0]
            ED_skull = np.linalg.norm(skull_coordinates - intersection_coordinate, axis=1)
            closest = int(np.argmin(ED_skull))
            skull_angle_list.append(math.degrees(
                PlanTUS.angle_between_vectors(skin_normals[i], skull_normals[closest])))
        except IndexError:
            # no skull within the 40 mm ray: NaN scores as worst (was 0 = best)
            skull_angle_list.append(np.nan)
    skull_angles = np.array(skull_angle_list)
    _write_metric(PlanTUS, "skin_skull_angles", skull_angles, output_path)

    log("Computing skull thickness metric...")
    labels_path = os.path.join(os.path.dirname(mesh_path), TISSUE_LABELS_FILENAME)
    if os.path.isfile(labels_path):
        thickness = skull_thickness_along_normals(
            skin_coordinates, skin_normals, labels_path)
        log(f"  transcranial path found for {int(np.isfinite(thickness).sum())} / "
            f"{len(thickness)} skin vertices")
    else:
        log(f"Warning: {labels_path} not found; skull thickness is not "
            "available and its weight is ignored.")
        thickness = np.full(len(skin_coordinates), np.nan)
    _write_metric(PlanTUS, "skull_thickness", thickness, output_path)

    avoidance_path = output_path + "/avoidance_skin.func.gii"
    if os.path.isfile(avoidance_path):
        avoidance = nib.load(avoidance_path).darrays[0].data
    else:
        log("Warning: avoidance mask not found, using all-ones (no exclusion).")
        avoidance = np.ones(len(skin_coordinates))

    return {
        "skin_coordinates":          np.asarray(skin_coordinates, dtype=float),
        "skin_normals":              np.asarray(skin_normals, dtype=float),
        "skin_faces":                np.asarray(
            nib.load(output_path + "/skin.surf.gii").darrays[1].data, dtype=np.int32),
        "target_center":             np.asarray(target_center, dtype=float),
        "skin_target_distances":     np.asarray(distances, dtype=float),
        "skin_target_angles":        np.asarray(angles, dtype=float),
        "skin_target_intersections": np.asarray(intersections, dtype=float),
        "skin_skull_angles":         np.asarray(skull_angles, dtype=float),
        "skull_thickness":           np.asarray(thickness, dtype=float),
        "avoidance":                 np.asarray(avoidance, dtype=float),
    }


def write_planning_scene(PlanTUS, template, output_path, t1_path, mask_filename,
                         with_thickness=True):
    """Write <output_path>/scene.scene showing the metric maps (for wb_view).

    with_thickness also loads the skull-thickness map (add_scene_metric).
    """
    if not template:
        return None
    names = [
        'SKIN_SURFACE_FILENAME',   'SKIN_SURFACE_FILEPATH',
        'SKULL_SURFACE_FILENAME',  'SKULL_SURFACE_FILEPATH',
        'DISTANCES_FILENAME',      'DISTANCES_FILEPATH',
        'INTERSECTION_FILENAME',   'INTERSECTION_FILEPATH',
        'ANGLES_FILENAME',         'ANGLES_FILEPATH',
        'ANGLES_SKIN_SKULL_FILENAME', 'ANGLES_SKIN_SKULL_FILEPATH',
        'DISTANCES_MAX_FILENAME',  'DISTANCES_MAX_FILEPATH',
        'T1_FILENAME',             'T1_FILEPATH',
        'MASK_FILENAME',           'MASK_FILEPATH',
    ]
    values = [
        'skin.surf.gii',             './skin.surf.gii',
        'skull.surf.gii',            './skull.surf.gii',
        'distances_skin.func.gii',   './distances_skin.func.gii',
        'target_intersection_skin.func.gii', './target_intersection_skin.func.gii',
        'angles_skin.func.gii',      './angles_skin.func.gii',
        'skin_skull_angles_skin.func.gii', './skin_skull_angles_skin.func.gii',
        'distances_skin_thresholded.func.gii', './distances_skin_thresholded.func.gii',
        os.path.basename(t1_path),   os.path.relpath(t1_path, output_path),
        mask_filename,               './' + mask_filename,
    ]
    scene_path = output_path + "/scene.scene"
    PlanTUS.create_scene(template, scene_path, names, values)
    if with_thickness and os.path.isfile(
            os.path.join(output_path, "skull_thickness_skin.func.gii")):
        add_scene_metric(scene_path, "./skull_thickness_skin.func.gii")
    return scene_path


def add_scene_metric(scene_path, metric_path, structure="CORTEX_LEFT"):
    """Add a metric file to a scene's spec file so wb_view loads it.

    The PlanTUS planning template lists four criterion maps; the skull
    thickness map is added this way. wb_view loads every selected file of the
    scene's spec file, so the map appears as an overlay choice in each tab.
    Returns False (scene unchanged) if the spec-file list is not found.
    """
    import re

    with open(scene_path, "r", encoding="utf-8") as f:
        txt = f.read()
    m = re.search(r'(<ObjectArray Type="class" Name="dataFilesArray" Length=")(\d+)(">)'
                  r'(.*?)(\n([ \t]*)</ObjectArray>)', txt, re.S)
    if m is None:
        return False
    n = int(m.group(2))
    ind = m.group(6)
    element = (
        f'\n{ind}    <Element Index="{n}">'
        f'\n{ind}        <Object Type="class" Class="SpecFileDataFile" Name="specFileDataFile" Version="1">'
        f'\n{ind}            <Object Type="enumeratedType" Name="dataFileType">METRIC</Object>'
        f'\n{ind}            <Object Type="enumeratedType" Name="structure">{structure}</Object>'
        f'\n{ind}            <Object Type="pathName" Name="fileName">{metric_path}</Object>'
        f'\n{ind}            <Object Type="boolean" Name="selected">true</Object>'
        f'\n{ind}        </Object>'
        f'\n{ind}    </Element>')
    txt = (txt[:m.start()] + m.group(1) + str(n + 1) + m.group(3) + m.group(4)
           + element + m.group(5) + txt[m.end():])
    with open(scene_path, "w", encoding="utf-8") as f:
        f.write(txt)
    return True


def fix_placement_scene(output_path, vertex, t1_path):
    """Point a vertex's placement scene (vtx<N>/scene.scene) at the T1.

    PlanTUS writes the T1 as '../../../T1.nii.gz' (the m2m folder of its own
    layout), which is not where PRESTUS keeps it. Returns the scene path, or
    None if PlanTUS wrote no scene.
    """
    vtx_dir = os.path.join(output_path, "vtx" + str(vertex))
    scene_path = os.path.join(vtx_dir, "scene.scene")
    if not os.path.isfile(scene_path):
        return None
    with open(scene_path, "r", encoding="utf-8") as f:
        txt = f.read()
    if "../../../T1.nii.gz" not in txt:
        return scene_path       # PlanTUS 2.x writes absolute paths to the T1 it was given
    # The file name is also the scene's name for the T1's display state.
    marker = "\0T1_FILEPATH\0"
    txt = (txt.replace("../../../T1.nii.gz", marker)
              .replace("T1.nii.gz", os.path.basename(t1_path))
              .replace(marker, os.path.relpath(t1_path, vtx_dir)))
    with open(scene_path, "w", encoding="utf-8") as f:
        f.write(txt)
    return scene_path


# ---------------------------------------------------------------------------
# Placement and plantus_metrics.mat
# ---------------------------------------------------------------------------

def place_vertex(PlanTUS, vertex, output_path, mask_path, t1_path, s,
                 placement_scene_template):
    """Run PlanTUS.prepare_acoustic_simulation for one vertex.

    PlanTUS appends to its Localite XML snippet, so a re-run for the same
    vertex would accumulate duplicate markers; the file is removed first.
    PlanTUS 2.x would open its own results dialog after placing; it is
    skipped (the placement is shown in PRESTUS). Returns the vertex's
    placement scene (transducer and focus on skin and T1), or None if PlanTUS
    wrote none.
    """
    mask_name = os.path.basename(mask_path).replace(".nii.gz", "").replace(".nii", "")
    vtx = "vtx" + str(vertex)
    for xml_name in ("position_matrix_" + mask_name + "_" + vtx + "_Localite_XML.txt",   # 1.x
                     mask_name + "_" + vtx + "_TransducerPosition_Localite_dummyXML.txt"):  # 2.x
        xml_path = os.path.join(output_path, vtx, xml_name)
        if os.path.isfile(xml_path):
            os.remove(xml_path)
    args = (vertex, output_path, mask_path, t1_path,
            s["max_distance"], s["min_distance"],
            s["transducer_diameter"], s["max_angle"],
            s["plane_offset"], s["additional_offset"],
            "",  # no transducer surface model (generic model)
            s["focal_distance_list"], s["flhm_list"],
            placement_scene_template)
    if is_v2(PlanTUS):
        PlanTUS.prepare_acoustic_simulation(*args, skip_viewer=True)
    else:
        PlanTUS.prepare_acoustic_simulation(*args)
    return fix_placement_scene(output_path, vertex, t1_path)


def check_maps_for_vertex(PlanTUS, output_path, mask_path, vertex, log, tol_mm=1.0):
    """Errors (list of str) that stop placing a chosen vertex on existing maps.

    The vertex index refers to skin.surf.gii of an earlier full run in
    output_path. It is only valid while that surface and the target are the
    ones the maps were computed for: the target centre must match the one in
    plantus_metrics.mat. The target's 3D model is rebuilt from the current
    mask, as prepare_acoustic_simulation reads it.
    """
    skin = os.path.join(output_path, "skin.surf.gii")
    metrics = os.path.join(output_path, METRICS_FILENAME)
    if not (os.path.isfile(skin) and os.path.isfile(metrics)):
        return [f"No PlanTUS maps in {output_path}: run PlanTUS once without "
                "placement.plantus.vertex (and with save_metrics) first."]
    d = scipy.io.loadmat(metrics, variable_names=["skin_coordinates", "target_center"])
    n = len(d["skin_coordinates"])
    if not 0 <= vertex < n:
        return [f"placement.plantus.vertex = {vertex} is not on the skin surface "
                f"({n} vertices, 0-based)."]
    now = np.asarray(PlanTUS.roi_center_of_gravity(mask_path), dtype=float).ravel()
    then = np.asarray(d["target_center"], dtype=float).ravel()
    shift = float(np.linalg.norm(now - then))
    if shift > tol_mm:
        return [f"The target moved {shift:.1f} mm since the PlanTUS maps were computed; "
                "clear placement.plantus.vertex and run PlanTUS again."]
    PlanTUS.stl_from_nii(mask_path, 0.25)
    log(f"Placing vertex {vertex} on the maps in {output_path}.")
    return []


def remove_stale_metrics(output_path, log):
    """Delete plantus_metrics.mat and its evaluation left by an earlier run.

    They would otherwise describe a different placement than the one this run
    writes (e.g. an old headless result next to a new GUI placement).
    """
    for name in (METRICS_FILENAME, EVALUATION_FILENAME, SWEEP_FILENAME,
                 MAPS_FILENAME, SENSITIVITY_FILENAME):
        p = os.path.join(output_path, name)
        if os.path.isfile(p):
            os.remove(p)
            log(f"Removed {name} from an earlier run.")


def save_metrics(output_path, m, s, score, best_vertex, selected_vertex,
                 selection_mode, log):
    """Write plantus_metrics.mat (per-vertex criteria, score and selection).

    best_vertex is the composite-score optimum; selected_vertex is the vertex
    actually placed (== best_vertex in headless mode). Both are 0-based.
    Failures are reported but do not abort, since the placement is already
    written.
    """
    path = os.path.join(output_path, METRICS_FILENAME)
    try:
        scipy.io.savemat(path, {
            **m,  # skin_faces is 0-based
            "score":            np.asarray(score, dtype=float),
            "best_vertex":      -1 if best_vertex is None else int(best_vertex),
            "selected_vertex":  int(selected_vertex),
            "selection_mode":   selection_mode,
            "min_distance":     s["min_distance"],
            "max_distance":     s["max_distance"],
            "optimal_distance": s["optimal_distance"],
            "weights":          np.array(s["weights"]),
            "weight_names":     np.array(["skin_target_distances",
                                          "skin_target_angles",
                                          "skin_target_intersections",
                                          "skin_skull_angles",
                                          "skull_thickness"], dtype=object),
        })
        log(f"Per-vertex criteria saved: {path}")
    except Exception as exc:  # noqa: BLE001 — metrics are optional
        log(f"Warning: could not save {METRICS_FILENAME}: {exc}")


def run_evaluation(output_path, s, log):
    """Evaluate the placement from plantus_metrics.mat (plantus_evaluation.py).

    Writes plantus_evaluation.json, the weight-sweep CSV and two figures next
    to the metrics file; generate_placement_report.m reads them.  Never
    raises: the placement itself is already written.
    """
    if not s["evaluate"]:
        return
    path = os.path.join(output_path, METRICS_FILENAME)
    if not os.path.isfile(path):
        log("Placement evaluation skipped: no plantus_metrics.mat "
            "(needs save_metrics = true).")
        return
    try:
        import plantus_evaluation
        plantus_evaluation.evaluate(path, output_path, log=log)
    except Exception as exc:  # noqa: BLE001 — evaluation is optional
        log(f"Warning: placement evaluation failed: {exc}")


def update_selected_vertex(output_path, vertex, selection_mode, log):
    """Record a re-placed vertex in an existing plantus_metrics.mat."""
    path = os.path.join(output_path, METRICS_FILENAME)
    if not os.path.isfile(path):
        return
    try:
        d = scipy.io.loadmat(path)
        d = {k: v for k, v in d.items() if not k.startswith("__")}
        d["selected_vertex"] = int(vertex)
        d["selection_mode"] = selection_mode
        scipy.io.savemat(path, d)
        log(f"Selected vertex updated in {path}")
    except Exception as exc:  # noqa: BLE001
        log(f"Warning: could not update {METRICS_FILENAME}: {exc}")
