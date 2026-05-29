#!/usr/bin/env python3
"""
prestus_plantus_launcher.py  –  headless PlanTUS placement for PRESTUS
=======================================================================
Replaces the interactive PlanTUS_wrapper.py with a non-GUI, CLI-driven
pipeline that automatically selects the optimal skin vertex using a
weighted composite score.

Usage
-----
python prestus_plantus_launcher.py \
    <t1_filepath> <simnibs_mesh_filepath> \
    <target_mask_filepath> <tx_config_yaml> \
    [--do_only_trajectory <vertex_index>]

The script writes one *_Localite.mat under:
  <mesh_dir>/PlanTUS/<target_name>/vtx<N>/

Exit codes
----------
0  success
1  error (traceback printed to stderr)
"""

import sys
import os
import argparse
import math
import numpy as np
import scipy.io

# Compatibility shims for NumPy 2.0 removals used by PlanTUS
if not hasattr(np, 'asfarray'):
    np.asfarray = lambda a, dtype=float: np.asarray(a, dtype=dtype)
if not hasattr(np, 'float'):
    np.float = float
if not hasattr(np, 'string_'):
    np.string_ = np.bytes_

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Headless PlanTUS transducer placement")
    p.add_argument("t1_filepath")
    p.add_argument("simnibs_mesh_filepath")
    p.add_argument("target_mask_filepath")
    p.add_argument("tx_config_yaml")
    p.add_argument("--do_only_trajectory", type=int, default=-1,
                   help="Skip metric computation and run prepare_acoustic_simulation "
                        "for a specific vertex index directly (for re-runs).")
    p.add_argument("--plantus_root", default="",
                   help="Path to PlanTUS installation directory (contains code/PlanTUS.py).")
    p.add_argument("--output_dir", default="",
                   help="Directory where PlanTUS output (including *_Localite.mat) is written. "
                        "Defaults to <mesh_dir>/PlanTUS/<target_name>/ when not set.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# YAML config loader (yaml or pyyaml; fall back to manual parse)
# ---------------------------------------------------------------------------

def _ensure_list(v):
    """Return v as a list; wraps scalars that YAML parsed as a single number."""
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
            # Try numeric / list
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
    out[~valid] = 1.0
    if invert:
        out = 1.0 - out
    return out


def compute_composite_score(distances, angles, intersections,
                             skull_angles, avoidance,
                             w_dist, w_ang, w_inter, w_skull, w_thick,
                             max_distance, min_distance=0.0,
                             optimal_distance=None):
    """Return a per-vertex composite cost (lower = better).

    Vertices with avoidance==0, distance<min_distance, or distance>max_distance
    get score=inf.  skull_thickness weight is reserved; if no metric is provided
    the weight is redistributed equally to the other four.

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
    # skull_thickness: no metric available — absorb weight into others
    total_w = w_dist + w_ang + w_inter + w_skull + w_thick
    if total_w <= 0:
        total_w = 1.0
    score = (w_dist  * s_dist +
             w_ang   * s_ang  +
             w_inter * s_inter +
             w_skull * s_skull) / total_w  # w_thick not used

    # Mask out unavailable / out-of-range vertices
    score[np.array(avoidance) == 0] = np.inf
    score[np.array(distances) < min_distance] = np.inf
    score[np.array(distances) > max_distance] = np.inf

    return score


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    t1_path    = args.t1_filepath
    mesh_path  = args.simnibs_mesh_filepath
    mask_path  = args.target_mask_filepath
    yaml_path  = args.tx_config_yaml

    for p, label in [(t1_path, "T1"), (mesh_path, "mesh"),
                     (mask_path, "target mask"), (yaml_path, "YAML config")]:
        if not os.path.isfile(p):
            sys.exit(f"[prestus_plantus_launcher] {label} not found: {p}")

    cfg = load_yaml(yaml_path)

    # Prepend connectome_wb_path so PlanTUS.py can call wb_command by name.
    # If not set, fall back to the directory of the running Python binary so that
    # a self-contained conda env (e.g. PRESTUS_env_4.6.0) requires no extra config.
    wb_path = cfg.get("connectome_wb_path", "") or os.path.dirname(sys.executable)
    os.environ["PATH"] = wb_path + os.pathsep + os.environ.get("PATH", "")

    max_distance        = float(cfg.get("max_distance",        70.0))
    min_distance        = float(cfg.get("min_distance",        30.0))
    transducer_diameter = float(cfg.get("transducer_diameter", 60.0))
    max_angle           = float(cfg.get("max_angle",           10.0))
    plane_offset        = float(cfg.get("plane_offset",         0.0))
    additional_offset   = float(cfg.get("additional_offset",    0.0))
    focal_distance_list = _ensure_list(cfg.get("focal_distance_list",  []))
    flhm_list           = _ensure_list(cfg.get("flhm_list",            []))
    w_dist  = float(cfg.get("weight_skin_target_distances",     0.2))
    w_ang   = float(cfg.get("weight_skin_target_angles",        0.2))
    w_inter = float(cfg.get("weight_skin_target_intersections", 0.2))
    w_skull = float(cfg.get("weight_skin_skull_angles",         0.2))
    w_thick = float(cfg.get("weight_skull_thickness",           0.2))

    target_mask_filename = os.path.split(mask_path)[1]
    target_name = target_mask_filename.replace(".nii.gz", "").replace(".nii", "")

    if args.output_dir:
        output_path = os.path.abspath(args.output_dir)
    else:
        output_path = os.path.join(os.path.split(mesh_path)[0],
                                   "PlanTUS", target_name)
    os.makedirs(output_path, exist_ok=True)

    # Locate PlanTUS code directory.
    # Priority: --plantus_root CLI arg > sibling of PRESTUS > relative fallback.
    candidates_plantus = []
    if args.plantus_root:
        candidates_plantus.append(os.path.join(args.plantus_root, "code"))
    prestus_root = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", ".."))
    candidates_plantus.append(os.path.join(prestus_root, "PlanTUS", "code"))
    candidates_plantus.append(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "..", "tools", "PlanTUS", "code"))

    plantus_code_path = None
    for cand in candidates_plantus:
        cand = os.path.normpath(cand)
        if os.path.isfile(os.path.join(cand, "PlanTUS.py")):
            plantus_code_path = cand
            break
    if plantus_code_path is None:
        sys.exit("[prestus_plantus_launcher] PlanTUS.py not found. Tried:\n  " +
                 "\n  ".join(candidates_plantus))

    sys.path.insert(0, plantus_code_path)
    os.chdir(plantus_code_path)
    import PlanTUS

    # Patch PlanTUS.compute_FLHM_for_focal_distance to handle fewer than 4
    # calibration points (cubic curve_fit requires M >= N = 4).
    # Patch PlanTUS.create_scene to skip gracefully when no template is provided
    # (headless mode passes "" and PlanTUS would crash on shutil.copy of empty path)
    _orig_create_scene = PlanTUS.create_scene
    def _patched_create_scene(scene_template_filepath, *args, **kwargs):
        if not scene_template_filepath:
            return
        return _orig_create_scene(scene_template_filepath, *args, **kwargs)
    PlanTUS.create_scene = _patched_create_scene

    # PlanTUS.prepare_acoustic_simulation unconditionally calls
    # os.system("wb_view ... scene.scene") after create_scene. In headless mode
    # scene.scene is never written (create_scene is skipped above), so wb_view
    # would open an error dialog. Suppress it by no-oping os.system within the
    # PlanTUS module for wb_view calls.
    import os as _os
    _orig_os_system = _os.system
    def _patched_os_system(cmd):
        if "wb_view" in cmd:
            return 0
        return _orig_os_system(cmd)
    # Patch os.system on the shared os module so PlanTUS's internal calls are intercepted
    _os.system = _patched_os_system

    _orig_compute_flhm = PlanTUS.compute_FLHM_for_focal_distance
    def _patched_compute_flhm(focal_distance, focal_distance_list, flhm_list):
        if len(focal_distance_list) >= 4:
            return _orig_compute_flhm(focal_distance, focal_distance_list, flhm_list)
        # Fewer than 4 points: use linear interpolation / extrapolation instead
        return float(np.interp(focal_distance, focal_distance_list, flhm_list))
    PlanTUS.compute_FLHM_for_focal_distance = _patched_compute_flhm

    # ------------------------------------------------------------------ #
    # Fast path: re-run for a specific vertex without recomputing metrics  #
    # ------------------------------------------------------------------ #
    if args.do_only_trajectory >= 0:
        vtx = args.do_only_trajectory
        print(f"[prestus_plantus_launcher] --do_only_trajectory {vtx}: "
              f"skipping metric computation.")
        PlanTUS.prepare_acoustic_simulation(
            vtx, output_path, mask_path, t1_path,
            max_distance, min_distance,
            transducer_diameter, max_angle, plane_offset, additional_offset,
            "",  # no transducer surface model
            focal_distance_list, flhm_list,
            "")  # no placement scene template
        print("[prestus_plantus_launcher] Done.")
        return

    # ------------------------------------------------------------------ #
    # Step 1 – Convert SimNIBS mesh to skin and skull surfaces             #
    # ------------------------------------------------------------------ #
    print("[prestus_plantus_launcher] Converting SimNIBS mesh to surfaces...")
    import shutil
    dst_mask = os.path.join(output_path, target_mask_filename)
    if os.path.abspath(mask_path) != os.path.abspath(dst_mask):
        shutil.copy(mask_path, dst_mask)
    mask_path = dst_mask

    PlanTUS.convert_simnibs_mesh_to_surface(mesh_path, [1005], "skin",  output_path)
    PlanTUS.add_structure_information(output_path + "/skin.surf.gii",  "CORTEX_LEFT")
    PlanTUS.convert_simnibs_mesh_to_surface(mesh_path, [1007, 1008], "skull", output_path)
    PlanTUS.add_structure_information(output_path + "/skull.surf.gii", "CORTEX_RIGHT")

    # ------------------------------------------------------------------ #
    # Step 2 – Avoidance mask                                              #
    # ------------------------------------------------------------------ #
    print("[prestus_plantus_launcher] Computing avoidance mask...")
    PlanTUS.create_avoidance_mask(
        mesh_path, output_path + "/skin.surf.gii", transducer_diameter / 2)

    # ------------------------------------------------------------------ #
    # Step 3 – Per-vertex metrics                                          #
    # ------------------------------------------------------------------ #
    print("[prestus_plantus_launcher] Computing skin-to-target metrics...")

    target_center = PlanTUS.roi_center_of_gravity(mask_path)

    skin_target_distances = PlanTUS.distance_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)

    skin_coordinates, skin_normals = PlanTUS.compute_surface_metrics(
        output_path + "/skin.surf.gii")
    skin_target_vectors = PlanTUS.vectors_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)

    skin_target_angles = []
    for i in np.arange(len(skin_target_vectors)):
        skin_target_angles.append(math.degrees(
            PlanTUS.angle_between_vectors(skin_target_vectors[i], skin_normals[i])))
    skin_target_angles = np.abs(np.array(skin_target_angles))

    print("[prestus_plantus_launcher] Computing skin-target intersection metric...")
    PlanTUS.stl_from_nii(mask_path, 0.25)
    skin_target_intersections = PlanTUS.compute_vector_mesh_intersections(
        skin_coordinates, skin_normals,
        output_path + "/" + target_name + "_3Dmodel.stl", 200)

    skin_target_intersection_values = []
    for i in np.arange(len(skin_target_intersections)):
        segs = skin_target_intersections[i]
        if len(segs) <= 1:
            skin_target_intersection_values.append(0)
        elif len(segs) in (2, 3):
            skin_target_intersection_values.append(
                np.linalg.norm(np.array(segs[1]) - np.array(segs[0])))
        elif len(segs) == 4:
            skin_target_intersection_values.append(
                np.linalg.norm(np.array(segs[1]) - np.array(segs[0])) +
                np.linalg.norm(np.array(segs[3]) - np.array(segs[2])))
        else:
            skin_target_intersection_values.append(np.nan)
    skin_target_intersection_values = np.array(skin_target_intersection_values)

    print("[prestus_plantus_launcher] Computing skin-skull angle metric...")
    skull_coordinates, skull_normals = PlanTUS.compute_surface_metrics(
        output_path + "/skull.surf.gii")
    skin_skull_intersections = PlanTUS.compute_vector_mesh_intersections(
        skin_coordinates, skin_normals, output_path + "/skull.stl", 40)

    skin_skull_angle_list = []
    for i in np.arange(len(skin_coordinates)):
        try:
            intersection_coordinate = skin_skull_intersections[i][0]
            ED_skull = np.linalg.norm(skull_coordinates - intersection_coordinate, axis=1)
            closest = int(np.argmin(ED_skull))
            angle = math.degrees(
                PlanTUS.angle_between_vectors(skin_normals[i], skull_normals[closest]))
            skin_skull_angle_list.append(angle)
        except Exception:
            skin_skull_angle_list.append(0)
    skin_skull_angles = np.array(skin_skull_angle_list)

    # ------------------------------------------------------------------ #
    # Step 4 – Load avoidance mask from saved func.gii                    #
    # ------------------------------------------------------------------ #
    import nibabel as nib
    avoidance_path = output_path + "/avoidance_skin.func.gii"
    if os.path.isfile(avoidance_path):
        avoidance_img = nib.load(avoidance_path)
        avoidance_mask_arr = avoidance_img.darrays[0].data
    else:
        print("[prestus_plantus_launcher] Warning: avoidance mask not found, "
              "using all-ones (no exclusion).")
        avoidance_mask_arr = np.ones(len(skin_coordinates))

    # ------------------------------------------------------------------ #
    # Step 5 – Composite scoring and vertex selection                      #
    # ------------------------------------------------------------------ #
    print("[prestus_plantus_launcher] Computing composite score...")
    optimal_distance = float(cfg.get("optimal_distance", (min_distance + max_distance) / 2))
    score = compute_composite_score(
        skin_target_distances, skin_target_angles,
        skin_target_intersection_values, skin_skull_angles,
        avoidance_mask_arr,
        w_dist, w_ang, w_inter, w_skull, w_thick,
        max_distance, min_distance=min_distance,
        optimal_distance=optimal_distance)

    n_valid = np.sum(np.isfinite(score))
    print(f"  Valid candidate vertices: {n_valid} / {len(score)}")
    if n_valid == 0:
        sys.exit("[prestus_plantus_launcher] No valid skin vertex found "
                 "(all excluded by avoidance mask or max_distance). "
                 "Consider increasing max_distance or mask_radius_mm.")

    best_vertex = int(np.argmin(score))
    print(f"  Best vertex: {best_vertex}  (score = {score[best_vertex]:.4f})")
    print(f"    distance to target : {skin_target_distances[best_vertex]:.1f} mm")
    print(f"    skin-target angle  : {skin_target_angles[best_vertex]:.1f} deg")
    print(f"    skin-skull angle   : {skin_skull_angles[best_vertex]:.1f} deg")

    # ------------------------------------------------------------------ #
    # Step 6 – Generate Localite position matrix                          #
    # ------------------------------------------------------------------ #
    print(f"[prestus_plantus_launcher] Generating placement for vertex {best_vertex}...")
    PlanTUS.prepare_acoustic_simulation(
        best_vertex, output_path, mask_path, t1_path,
        max_distance, min_distance,
        transducer_diameter, max_angle, plane_offset, additional_offset,
        "",  # no transducer surface model
        focal_distance_list, flhm_list,
        "")  # no placement scene template

    print("[prestus_plantus_launcher] Done.")


if __name__ == "__main__":
    main()
