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
Each click prints the vertex's criteria, score and rank.

Confirmation (--selection):
  review      when wb_view closes, the last selected vertex is placed, its
              transducer is shown in wb_view and a dialog asks to use it,
              choose another vertex (wb_view reopens, no recomputation) or
              cancel (needs PyQt5, shipped with SimNIBS; used by PRESTUS)
  prompt      yes/no console prompt after each click (needs a terminal stdin)
  last_click  last selected vertex is placed when wb_view closes
  auto        (default) prompt if stdin is a terminal, else review if PyQt5
              is available, else last_click

Usage
-----
python prestus_plantus_gui_launcher.py \
    <t1_filepath> <simnibs_mesh_filepath> \
    <target_mask_filepath> <tx_config_yaml> \
    --plantus_root <path_to_PlanTUS_installation> \
    [--selection auto|review|prompt|last_click]

Requirements
------------
- pynput          (pip install pynput)
- Connectome Workbench wb_view reachable via PATH or connectome_wb_path YAML field
- All PlanTUS Python dependencies (nibabel, scipy, nilearn, ...)
"""

import sys
import os
import argparse
import html
import re
import shlex
import subprocess
import threading
import time

import numpy as np

import plantus_compat
import plantus_common as pc

# A wb_view that exits this soon without a click is taken as having failed to
# load the scene (e.g. a Workbench version that rejects an added map).
SCENE_FAILURE_S = 15


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
    p.add_argument("--selection", choices=["auto", "review", "prompt", "last_click"],
                   default="auto",
                   help="How a clicked vertex is confirmed. 'review': the last "
                        "selected vertex is placed when wb_view closes and a "
                        "dialog asks to use it or choose again (needs PyQt5). "
                        "'prompt': yes/no at the console (needs an interactive "
                        "terminal). 'last_click': the last selected vertex is "
                        "placed when wb_view closes. 'auto' (default): 'prompt' "
                        "if stdin is a terminal, else 'review' if PyQt5 is "
                        "available, else 'last_click'.")
    return p.parse_args()


def log(msg):
    print(f"[prestus_plantus_gui_launcher] {msg}")


def have_qt():
    try:
        import PyQt5.QtWidgets  # noqa: F401
        return True
    except Exception:  # noqa: BLE001 — any import failure means no dialog
        return False


def review_dialog(vertex, lines):
    """Ask whether to use the placed vertex: 'use', 'again' or 'cancel'."""
    from PyQt5 import QtCore, QtWidgets

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(["PlanTUS"])
    box = QtWidgets.QMessageBox()
    box.setWindowTitle("PlanTUS: review placement")
    box.setIcon(QtWidgets.QMessageBox.Question)
    box.setTextFormat(QtCore.Qt.RichText)
    box.setText(f"<b>Transducer placed at vertex {vertex}.</b><br>"
                "Connectome Workbench shows the transducer and focus.")
    box.setInformativeText("<pre>" + html.escape("\n".join(lines)) + "</pre>")
    use = box.addButton("Use this position", QtWidgets.QMessageBox.AcceptRole)
    again = box.addButton("Choose another vertex", QtWidgets.QMessageBox.ActionRole)
    cancel = box.addButton("Cancel", QtWidgets.QMessageBox.RejectRole)
    box.setDefaultButton(use)
    box.setEscapeButton(cancel)
    box.setWindowFlags(box.windowFlags() | QtCore.Qt.WindowStaysOnTopHint)
    box.raise_()
    box.activateWindow()
    box.exec_()
    app.processEvents()
    clicked = box.clickedButton()
    if clicked is use:
        return "use"
    if clicked is again:
        return "again"
    return "cancel"


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

    # prepare_acoustic_simulation ends with os.system("wb_view ... &") on the
    # placement scene, before PRESTUS has pointed that scene at the T1. Its
    # wb_view calls are recorded instead and the viewer is opened below.
    import os as _os
    _orig_os_system = _os.system
    wb_view_calls = []

    def _patched_os_system(cmd):
        if "wb_view" in cmd:
            wb_view_calls.append(cmd)
            return 0
        return _orig_os_system(cmd)
    _os.system = _patched_os_system

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
    with_thickness = True

    # Selection mode.  "prompt" asks yes/no at the console after each click,
    # which needs an interactive stdin.  When launched from MATLAB (via
    # run_command) stdin is a pipe nobody writes to, so input() would block
    # forever -- there, "review" (a Qt dialog) or "last_click" is used.
    selection = args.selection
    if selection == "auto":
        if sys.stdin is not None and sys.stdin.isatty():
            selection = "prompt"
        else:
            selection = "review" if have_qt() else "last_click"
    elif selection == "review" and not have_qt():
        log("PyQt5 not available in this Python; the last selected vertex is "
            "placed without review.")
        selection = "last_click"

    placed = []
    viewer = None

    def open_viewer(placement_scene):
        """Show a placement scene (transducer, focus) in a separate wb_view."""
        nonlocal viewer
        close_viewer()
        if placement_scene:
            viewer = subprocess.Popen(
                ["wb_view", "-logging", "OFF", placement_scene],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                cwd=os.path.dirname(placement_scene))
        wb_view_calls.clear()

    def close_viewer():
        nonlocal viewer
        if viewer is not None and viewer.poll() is None:
            viewer.terminate()
        viewer = None

    def summary(vertex):
        lines = pc.vertex_summary(m, s, score, vertex)
        if best_vertex is not None and vertex != best_vertex:
            lines.append(f"Composite optimum     vertex {best_vertex}")
        return lines

    def place(vertex):
        print(f"Generating placement for vertex {vertex}...")
        placement_scene = pc.place_vertex(PlanTUS, vertex, output_path,
                                          target_roi_filepath, t1_path, s,
                                          placement_template)
        placed.append(vertex)
        print(f"Placement written for vertex {vertex}.")
        return placement_scene

    def planning_session():
        """Open the planning scene; return the last clicked vertex (or None)."""
        command = ["wb_view", "-logging", "FINER", str(scene_path)]
        pattern = re.compile(r"Switched vertex to triangle nearest vertex\s+(\.\d+)")
        clicked = None
        process_line = False

        from pynput import mouse

        def on_click(x, y, button, pressed):
            nonlocal process_line
            if pressed:
                process_line = True

        listener = mouse.Listener(on_click=on_click)
        listener.start()

        # wb_view's stdout is discarded: only stderr (the logging stream) is
        # read, and an unread stdout pipe could fill up and stall wb_view.
        started = time.monotonic()
        process = subprocess.Popen(
            command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
            cwd=output_path, text=True)

        def read_output():
            nonlocal clicked, process_line
            while True:
                output = process.stderr.readline()
                if output == "" and process.poll() is not None:
                    break
                if process_line:
                    match = pattern.search(output)
                    if match:
                        clicked = int(match.group(1).replace(".", ""))
                        print(f"\nSwitched to vertex {clicked}")
                        for line in summary(clicked):
                            print("  " + line)
                        if selection == "prompt":
                            response = input(
                                f"Generate transducer placement for vertex "
                                f"{clicked}? (yes/no): ").strip().lower()
                            if response == "yes":
                                open_viewer(place(clicked))
                            else:
                                print("No action taken — click another vertex.")
                        else:
                            print("  Selected. Close wb_view to continue, or click another vertex.")
                        process_line = False

        output_thread = threading.Thread(target=read_output)
        output_thread.start()
        process.wait()
        output_thread.join()
        listener.stop()
        failed = (clicked is None and process.returncode != 0
                  and time.monotonic() - started < SCENE_FAILURE_S)
        return clicked, failed

    print("[prestus_plantus_gui_launcher] Opening Connectome Workbench...")
    print("  Click on a skin vertex in wb_view to select the transducer position.")
    print("  Each click prints the vertex's criteria and rank here.")
    if selection == "prompt":
        print("  Type 'yes' at the prompt to confirm and generate placement, then")
        print("  close wb_view to finish.")
    elif selection == "review":
        print("  Close wb_view to place the last selected vertex and review it;")
        print("  close without clicking to abort.")
    else:
        print("  Click again to change the selection; the last selected vertex is")
        print("  used. Close wb_view to confirm (close without clicking to abort).")

    last = None
    accepted = None
    while True:
        clicked, failed = planning_session()
        if failed and with_thickness:
            # Retry once with the template's four maps only.
            log("wb_view exited at once; reopening without the skull-thickness map.")
            scene_path = pc.write_planning_scene(PlanTUS, planning_template, output_path,
                                                 t1_path, target_mask_filename,
                                                 with_thickness=False)
            with_thickness = False
            continue
        if clicked is not None:
            last = clicked
        if last is None:
            sys.exit("[prestus_plantus_gui_launcher] wb_view closed without vertex selection. "
                     "No Localite.mat written.")
        if selection == "prompt":
            break
        if selection == "last_click":
            open_viewer(place(last))
            break

        # review: place, show the transducer, ask
        placement_scene = place(last)
        open_viewer(placement_scene)
        if placement_scene is None:
            log("PlanTUS wrote no placement scene; the transducer is not shown.")
        choice = review_dialog(last, summary(last))
        if choice == "use":
            accepted = last
            break
        close_viewer()
        if choice == "cancel":
            sys.exit("[prestus_plantus_gui_launcher] Placement cancelled. "
                     "PRESTUS does not use the vertices placed during review.")
        print("[prestus_plantus_gui_launcher] Reopening wb_view: click another vertex.")

    # Exit non-zero when nothing was confirmed, so PRESTUS does not fall back
    # on a *Localite.mat left over from an earlier run.
    if not placed:
        sys.exit("[prestus_plantus_gui_launcher] No placement confirmed. "
                 "No Localite.mat written.")

    # PRESTUS uses the most recently written placement: the last one placed,
    # which in review mode is the accepted vertex.
    final = accepted if accepted is not None else placed[-1]
    for line in summary(final):
        print("  " + line)
    if s["save_metrics"]:
        pc.save_metrics(output_path, m, s, score, best_vertex, final,
                        "gui", log)
    pc.run_evaluation(output_path, s, log)

    print("[prestus_plantus_gui_launcher] Done.")


if __name__ == "__main__":
    main()
