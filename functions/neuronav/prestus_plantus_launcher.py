#!/usr/bin/env python3
"""
prestus_plantus_launcher.py  –  headless PlanTUS placement for PRESTUS
=======================================================================
Replaces the interactive PlanTUS_wrapper.py with a non-GUI, CLI-driven
pipeline that automatically selects the optimal skin vertex using a
weighted composite score. Metrics, scoring and outputs are shared with the
GUI launcher (plantus_common.py), so both modes write the same files.

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
import numpy as np

import plantus_compat
import plantus_common as pc
# Re-exported for scripts written against the pre-plantus_common launcher,
# which imported compute_composite_score from this module.
from plantus_common import load_yaml, compute_composite_score  # noqa: F401

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


def log(msg):
    print(f"[prestus_plantus_launcher] {msg}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    t1_path    = os.path.abspath(args.t1_filepath)
    mesh_path  = os.path.abspath(args.simnibs_mesh_filepath)
    mask_path  = os.path.abspath(args.target_mask_filepath)
    yaml_path  = os.path.abspath(args.tx_config_yaml)

    for p, label in [(t1_path, "T1"), (mesh_path, "mesh"),
                     (mask_path, "target mask"), (yaml_path, "YAML config")]:
        if not os.path.isfile(p):
            sys.exit(f"[prestus_plantus_launcher] {label} not found: {p}")

    s = pc.read_config(load_yaml(yaml_path))

    # Prepend connectome_wb_path so PlanTUS.py can call wb_command by name.
    # If not set, fall back to the directory of the running Python binary so that
    # a self-contained conda env (e.g. PRESTUS_env_4.6.0) requires no extra config.
    wb_path = s["connectome_wb_path"] or os.path.dirname(sys.executable)
    os.environ["PATH"] = wb_path + os.pathsep + os.environ.get("PATH", "")

    target_name = os.path.basename(mask_path).replace(".nii.gz", "").replace(".nii", "")

    if args.output_dir:
        output_path = os.path.abspath(args.output_dir)
    else:
        output_path = os.path.join(os.path.split(mesh_path)[0],
                                   "PlanTUS", target_name)
    os.makedirs(output_path, exist_ok=True)
    # Before any exit, so a run that fails (missing PlanTUS or dependencies)
    # cannot leave an earlier run's metrics behind as if they were its own.
    # The --do_only_trajectory path reuses the maps, not these files.
    if args.do_only_trajectory < 0:
        pc.remove_stale_metrics(output_path, log)

    plantus_root, tried = pc.resolve_plantus_root(args.plantus_root, __file__)
    if plantus_root is None:
        sys.exit("[prestus_plantus_launcher] PlanTUS.py not found. Tried:\n  " +
                 "\n  ".join(os.path.join(t, "code") for t in tried))
    version = pc.plantus_major_version(plantus_root)
    log(f"PlanTUS {version}.x at {plantus_root}")
    if pc.deprecation_notice(version):
        log(pc.deprecation_notice(version))
    missing = pc.check_dependencies(version)
    if missing:
        sys.exit("[prestus_plantus_launcher] Missing PlanTUS dependencies "
                 f"(python: {sys.executable}):\n  " + "\n  ".join(missing))
    plantus_code_path = os.path.join(plantus_root, "code")
    planning_template, placement_template = pc.scene_templates(plantus_root)

    sys.path.insert(0, plantus_code_path)
    os.chdir(plantus_code_path)
    import PlanTUS

    # Patch PlanTUS.create_scene to skip gracefully when no template is found
    # (PlanTUS would crash on shutil.copy of an empty path).
    _orig_create_scene = PlanTUS.create_scene
    def _patched_create_scene(scene_template_filepath, *args, **kwargs):
        if not scene_template_filepath:
            return
        return _orig_create_scene(scene_template_filepath, *args, **kwargs)
    PlanTUS.create_scene = _patched_create_scene

    # PlanTUS 1.x prepare_acoustic_simulation unconditionally calls
    # os.system("wb_view ... scene.scene") after writing the placement scene.
    # Headless mode writes the same scenes as the GUI mode (so they can be
    # opened later) but must not open a viewer: no-op wb_view calls.
    import os as _os
    _orig_os_system = _os.system
    def _patched_os_system(cmd):
        if "wb_view" in cmd:
            return 0
        return _orig_os_system(cmd)
    # Patch os.system on the shared os module so PlanTUS's internal calls are intercepted
    _os.system = _patched_os_system

    plantus_compat.patch_compute_flhm(PlanTUS)

    # ------------------------------------------------------------------ #
    # Fast path: re-run for a specific vertex without recomputing metrics  #
    # ------------------------------------------------------------------ #
    if args.do_only_trajectory >= 0:
        vtx = args.do_only_trajectory
        log(f"--do_only_trajectory {vtx}: skipping metric computation.")
        mask_path = pc.stage_mask(mask_path, output_path)
        errors = pc.check_maps_for_vertex(PlanTUS, output_path, mask_path, vtx, log)
        if errors:
            sys.exit("[prestus_plantus_launcher] " + "\n  ".join(errors))
        pc.place_vertex(PlanTUS, vtx, output_path, mask_path, t1_path, s, placement_template)
        if s["save_metrics"]:
            pc.update_selected_vertex(output_path, vtx, "manual", log)
        pc.run_evaluation(output_path, s, log)
        log("Done.")
        return

    mask_path = pc.stage_mask(mask_path, output_path)

    # ------------------------------------------------------------------ #
    # Steps 1–4 – surfaces, avoidance mask, per-vertex metrics            #
    # ------------------------------------------------------------------ #
    m = pc.compute_metrics(PlanTUS, mesh_path, mask_path, output_path, s, log)
    pc.write_planning_scene(PlanTUS, planning_template, output_path, t1_path,
                            os.path.basename(mask_path))

    # ------------------------------------------------------------------ #
    # Step 5 – Composite scoring and vertex selection                      #
    # ------------------------------------------------------------------ #
    log("Computing composite score...")
    score, best_vertex = pc.score_metrics(m, s)

    n_valid = np.sum(np.isfinite(score))
    print(f"  Valid candidate vertices: {n_valid} / {len(score)}")
    if best_vertex is None:
        sys.exit("[prestus_plantus_launcher] No valid skin vertex found "
                 "(all excluded by avoidance mask or max_distance). "
                 "Consider increasing max_distance or mask_radius_mm.")

    print(f"  Best vertex: {best_vertex}  (score = {score[best_vertex]:.4f})")
    print(f"    distance to target : {m['skin_target_distances'][best_vertex]:.1f} mm")
    print(f"    skin-target angle  : {m['skin_target_angles'][best_vertex]:.1f} deg")
    print(f"    skin-skull angle   : {m['skin_skull_angles'][best_vertex]:.1f} deg")
    print(f"    skull thickness    : {m['skull_thickness'][best_vertex]:.1f} mm")

    # ------------------------------------------------------------------ #
    # Step 6 – Generate Localite position matrix                          #
    # ------------------------------------------------------------------ #
    log(f"Generating placement for vertex {best_vertex}...")
    pc.place_vertex(PlanTUS, best_vertex, output_path, mask_path, t1_path, s,
                    placement_template)

    # Per-vertex criteria, so placements can be re-scored offline (e.g.
    # weight sweeps) without recomputing the ray casts.
    if s["save_metrics"]:
        pc.save_metrics(output_path, m, s, score, best_vertex, best_vertex,
                        "auto", log)
    pc.run_evaluation(output_path, s, log)

    log("Done.")


if __name__ == "__main__":
    main()
