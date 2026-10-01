#!/usr/bin/env python3
"""
prestus_plantus_gui_launcher.py  –  interactive (GUI) PlanTUS session for PRESTUS
==============================================================================
Reads the same YAML config produced by PRESTUS, sets up all variables
expected by PlanTUS_wrapper.py, then runs the interactive pipeline that opens
Connectome Workbench.  The user clicks a skin vertex in wb_view, confirms
placement, and the script writes a *_Localite.mat under the PlanTUS output
directory.  PRESTUS reads that file once the subprocess exits.

Metrics, scoring and outputs are shared with the headless launcher
(plantus_common.py): the composite score is computed here too, and
plantus_metrics.mat records both its optimum and the vertex placed.

Confirmation (--selection):
  prompt      yes/no console prompt after each click (needs a terminal stdin)
  last_click  last selected vertex is placed when wb_view closes (used by
              PRESTUS, which does not connect stdin)
  auto        (default) prompt if stdin is a terminal, else last_click

Usage
-----
python prestus_plantus_gui_launcher.py \
    <t1_filepath> <simnibs_mesh_filepath> \
    <target_mask_filepath> <tx_config_yaml> \
    --plantus_root <path_to_PlanTUS_installation> \
    [--selection auto|prompt|last_click]

Requirements
------------
- pynput          (pip install pynput)
- Connectome Workbench wb_view reachable via PATH or connectome_wb_path YAML field
- All PlanTUS Python dependencies (nibabel, scipy, nilearn, ...)
"""

import sys
import os
import argparse
import re
import subprocess
import threading

import numpy as np

import plantus_compat
import plantus_common as pc


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Interactive (GUI) PlanTUS session launched from PRESTUS")
    p.add_argument("t1_filepath")
    p.add_argument("simnibs_mesh_filepath")
    p.add_argument("target_mask_filepath")
    p.add_argument("tx_config_yaml")
    p.add_argument("--plantus_root", default="",
                   help="Path to PlanTUS installation directory "
                        "(must contain code/PlanTUS.py and resources/).")
    p.add_argument("--output_dir", default="",
                   help="Output directory (default: <mesh_dir>/PlanTUS/<mask_stem>, "
                        "the PlanTUS_wrapper.py location).")
    p.add_argument("--selection", choices=["auto", "prompt", "last_click"],
                   default="auto",
                   help="How a clicked vertex is confirmed. 'prompt': yes/no at "
                        "the console (needs an interactive terminal). "
                        "'last_click': the last selected vertex is placed when "
                        "wb_view closes (no stdin needed). 'auto' (default): "
                        "'prompt' if stdin is a terminal, else 'last_click'.")
    return p.parse_args()


def log(msg):
    print(f"[prestus_plantus_gui_launcher] {msg}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Line-buffer stdout so progress and instructions reach the caller (e.g.
    # MATLAB's run_command with Echo) as they are printed, not at exit.
    sys.stdout.reconfigure(line_buffering=True)
    args = parse_args()

    # Absolute paths: the script chdirs into the PlanTUS code folder below.
    t1_path    = os.path.abspath(args.t1_filepath)
    mesh_path  = os.path.abspath(args.simnibs_mesh_filepath)
    mask_path  = os.path.abspath(args.target_mask_filepath)
    yaml_path  = os.path.abspath(args.tx_config_yaml)

    for p, label in [(t1_path, "T1"), (mesh_path, "mesh"),
                     (mask_path, "target mask"), (yaml_path, "YAML config")]:
        if not os.path.isfile(p):
            sys.exit(f"[prestus_plantus_gui_launcher] {label} not found: {p}")

    s = pc.read_config(pc.load_yaml(yaml_path))

    target_mask_filename = os.path.basename(mask_path)
    target_roi_name = target_mask_filename.replace(".nii.gz", "").replace(".nii", "")

    # Resolved before the chdir below. Same output location as the headless launcher, so PRESTUS finds the
    # *Localite.mat in the work_dir it passes as --output_dir.
    if args.output_dir:
        output_path = os.path.abspath(args.output_dir)
    else:
        output_path = os.path.join(os.path.split(mesh_path)[0],
                                   "PlanTUS", target_roi_name)
    os.makedirs(output_path, exist_ok=True)

    # ---- Resolve PlanTUS root ------------------------------------------------
    plantus_main_folder, tried = pc.resolve_plantus_root(args.plantus_root, __file__)
    if plantus_main_folder is None:
        sys.exit("[prestus_plantus_gui_launcher] PlanTUS.py not found. Tried:\n  " +
                 "\n  ".join(tried))
    plantus_code_path = os.path.join(plantus_main_folder, "code")
    planning_template, placement_template = pc.scene_templates(plantus_main_folder)

    sys.path.insert(0, plantus_code_path)
    os.chdir(plantus_code_path)
    import PlanTUS

    # Same fix as the headless launcher: PlanTUS's cubic FLHM fit needs >= 4
    # calibration points, but the geometric fallback yields only one.
    plantus_compat.patch_compute_flhm(PlanTUS)

    pc.remove_stale_metrics(output_path, log)
    target_roi_filepath = pc.stage_mask(mask_path, output_path)

    # Connectome Workbench path
    if s["connectome_wb_path"]:
        # Prepend to PATH so wb_view / wb_command are found
        os.environ["PATH"] = s["connectome_wb_path"] + os.pathsep + os.environ.get("PATH", "")

    missing = pc.check_dependencies()
    if missing:
        sys.exit("[prestus_plantus_gui_launcher] Missing PlanTUS dependencies "
                 f"(python: {sys.executable}):\n  " + "\n  ".join(missing))

    # ---- Metric computation and scoring (shared with headless) -------------
    m = pc.compute_metrics(PlanTUS, mesh_path, target_roi_filepath, output_path, s, log)
    score, best_vertex = pc.score_metrics(m, s)
    if best_vertex is None:
        log("Note: no vertex is eligible under the composite score "
            "(avoidance mask / distance range); any clicked vertex can still be placed.")
    else:
        log(f"Composite-score optimum (headless choice): vertex {best_vertex} "
            f"(score = {score[best_vertex]:.4f})")

    # ---- Build and open planning scene in wb_view ---------------------------
    if not planning_template:
        sys.exit("[prestus_plantus_gui_launcher] Planning scene template not found in "
                 + os.path.join(plantus_main_folder, "resources", "scene_templates"))
    scene_path = pc.write_planning_scene(PlanTUS, planning_template, output_path,
                                         t1_path, target_mask_filename)

    # Selection mode.  "prompt" asks yes/no at the console after each click,
    # which needs an interactive stdin.  When launched from MATLAB (via
    # run_command) stdin is a pipe nobody writes to, so input() would block
    # forever -- there, "last_click" is used instead: every click updates the
    # selection and the last selected vertex is placed once wb_view closes.
    selection = args.selection
    if selection == "auto":
        selection = "prompt" if sys.stdin is not None and sys.stdin.isatty() else "last_click"

    print("[prestus_plantus_gui_launcher] Opening Connectome Workbench...")
    print("  Click on a skin vertex in wb_view to select the transducer position.")
    if selection == "prompt":
        print("  Type 'yes' at the prompt to confirm and generate placement, then")
        print("  close wb_view to finish.")
    else:
        print("  Click again to change the selection; the last selected vertex is")
        print("  used. Close wb_view to confirm (close without clicking to abort).")

    command = ["wb_view", "-logging", "FINER", str(scene_path)]
    pattern = re.compile(r"Switched vertex to triangle nearest vertex\s+(\.\d+)")
    triangle_number = None
    process_line = False

    from pynput import mouse

    def on_click(x, y, button, pressed):
        nonlocal process_line
        if pressed:
            process_line = True

    listener = mouse.Listener(on_click=on_click)
    listener.start()

    # wb_view's stdout is discarded: only stderr (the logging stream) is read,
    # and an unread stdout pipe could fill up and stall wb_view.
    process = subprocess.Popen(
        command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        cwd=output_path, text=True)

    placed = []

    def place(vertex):
        print(f"Generating placement for vertex {vertex}...")
        pc.place_vertex(PlanTUS, vertex, output_path, target_roi_filepath,
                        t1_path, s, placement_template)
        placed.append(vertex)
        if 0 <= vertex < len(score) and np.isfinite(score[vertex]):
            rank = int((score < score[vertex]).sum()) + 1
            print(f"Placement written for vertex {vertex} "
                  f"(composite score {score[vertex]:.4f}, rank {rank} "
                  f"of {int(np.isfinite(score).sum())} eligible).")
        else:
            print(f"Placement written for vertex {vertex} "
                  f"(not eligible under the composite score).")

    def read_output():
        nonlocal triangle_number, process_line
        while True:
            output = process.stderr.readline()
            if output == "" and process.poll() is not None:
                break
            if process_line:
                match = pattern.search(output)
                if match:
                    triangle_number = int(match.group(1).replace(".", ""))
                    print(f"\nSwitched to vertex {triangle_number}")
                    if selection == "prompt":
                        response = input(
                            f"Generate transducer placement for vertex "
                            f"{triangle_number}? (yes/no): ").strip().lower()
                        if response == "yes":
                            place(triangle_number)
                        else:
                            print("No action taken — click another vertex.")
                    else:
                        print("  Selected. Close wb_view to confirm, or click another vertex.")
                    process_line = False

    output_thread = threading.Thread(target=read_output)
    output_thread.start()

    process.wait()
    output_thread.join()
    listener.stop()

    if triangle_number is None:
        sys.exit("[prestus_plantus_gui_launcher] wb_view closed without vertex selection. "
                 "No Localite.mat written.")

    if selection == "last_click":
        place(triangle_number)

    # Exit non-zero when nothing was confirmed, so PRESTUS does not fall back
    # on a *Localite.mat left over from an earlier run.
    if not placed:
        sys.exit("[prestus_plantus_gui_launcher] No placement confirmed. "
                 "No Localite.mat written.")

    # PRESTUS uses the most recently written placement, i.e. the last one.
    if s["save_metrics"]:
        pc.save_metrics(output_path, m, s, score, best_vertex, placed[-1],
                        "gui", log)
    pc.run_evaluation(output_path, s, log)

    print("[prestus_plantus_gui_launcher] Done.")


if __name__ == "__main__":
    main()
