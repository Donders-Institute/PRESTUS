#!/usr/bin/env python3
"""
prestus_plantus_gui_launcher.py  –  interactive (GUI) PlanTUS session for PRESTUS
==============================================================================
Reads the same YAML config produced by PRESTUS, sets up all variables
expected by PlanTUS_wrapper.py, then runs the interactive pipeline that opens
Connectome Workbench.  The user clicks a skin vertex in wb_view, confirms
placement, and the script writes a *_Localite.mat under the PlanTUS output
directory.  PRESTUS reads that file once the subprocess exits.

Usage
-----
python prestus_plantus_gui_launcher.py \
    <t1_filepath> <simnibs_mesh_filepath> \
    <target_mask_filepath> <tx_config_yaml> \
    --plantus_root <path_to_PlanTUS_installation>

Requirements
------------
- pynput          (pip install pynput)
- Connectome Workbench wb_view reachable via PATH or connectome_wb_path YAML field
- All PlanTUS Python dependencies (nibabel, scipy, nilearn, ...)
"""

import sys
import os
import argparse
import shutil
import math
import numpy as np
import re
import subprocess
import threading


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
    return p.parse_args()


# ---------------------------------------------------------------------------
# YAML loader (same as headless launcher)
# ---------------------------------------------------------------------------

def load_yaml(path):
    try:
        import yaml
        with open(path, "r") as fh:
            return yaml.safe_load(fh)
    except ImportError:
        pass
    cfg = {}
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" not in line:
                continue
            k, _, v = line.partition(":")
            k, v = k.strip(), v.strip()
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
# Main
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
            sys.exit(f"[prestus_plantus_gui_launcher] {label} not found: {p}")

    cfg = load_yaml(yaml_path)

    # ---- Resolve PlanTUS root ------------------------------------------------
    candidates_root = []
    if args.plantus_root:
        candidates_root.append(args.plantus_root)
    prestus_root = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "..", ".."))
    candidates_root.append(os.path.join(prestus_root, "PlanTUS"))
    candidates_root.append(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "..", "tools", "PlanTUS"))

    plantus_main_folder = None
    for cand in candidates_root:
        cand = os.path.normpath(cand)
        if os.path.isfile(os.path.join(cand, "code", "PlanTUS.py")):
            plantus_main_folder = cand
            break
    if plantus_main_folder is None:
        sys.exit("[prestus_plantus_gui_launcher] PlanTUS.py not found. Tried:\n  " +
                 "\n  ".join(candidates_root))

    plantus_code_path = os.path.join(plantus_main_folder, "code")

    sys.path.insert(0, plantus_code_path)
    os.chdir(plantus_code_path)
    import PlanTUS

    # ---- Transducer / placement variables -----------------------------------
    max_distance        = float(cfg.get("max_distance",        70.0))
    min_distance        = float(cfg.get("min_distance",        30.0))
    transducer_diameter = float(cfg.get("transducer_diameter", 60.0))
    max_angle           = float(cfg.get("max_angle",           10.0))
    plane_offset        = float(cfg.get("plane_offset",         0.0))
    additional_offset   = float(cfg.get("additional_offset",    0.0))
    focal_distance_list = list(cfg.get("focal_distance_list",  []))
    flhm_list           = list(cfg.get("flhm_list",            []))

    target_mask_filename = os.path.split(mask_path)[1]
    target_roi_name = target_mask_filename.replace(".nii.gz", "").replace(".nii", "")

    output_path = os.path.join(os.path.split(mesh_path)[0],
                               "PlanTUS", target_roi_name)
    os.makedirs(output_path, exist_ok=True)

    # Copy mask to output directory (PlanTUS_wrapper.py convention)
    shutil.copy(mask_path, output_path + "/")
    target_roi_filepath = output_path + "/" + target_mask_filename

    # Scene templates and transducer model
    planning_scene_template_filepath = os.path.join(
        plantus_main_folder,
        "resources", "scene_templates",
        "TUSTransducerPlacementPlanning_TEMPLATE.scene")
    placement_scene_template_filepath = os.path.join(
        plantus_main_folder,
        "resources", "scene_templates",
        "TUSTransducerPlacement_TEMPLATE.scene")
    transducer_surface_model_filepath = ""  # use generic model

    # Connectome Workbench path
    connectome_wb_path = str(cfg.get("connectome_wb_path", "")).strip()
    if connectome_wb_path:
        # Prepend to PATH so wb_view / wb_command are found
        os.environ["PATH"] = connectome_wb_path + os.pathsep + os.environ.get("PATH", "")

    # ---- Metric computation (same as wrapper.py) ----------------------------
    print("[prestus_plantus_gui_launcher] Converting SimNIBS mesh to surfaces...")
    PlanTUS.convert_simnibs_mesh_to_surface(mesh_path, [1005], "skin",  output_path)
    PlanTUS.add_structure_information(output_path + "/skin.surf.gii",  "CORTEX_LEFT")
    PlanTUS.convert_simnibs_mesh_to_surface(mesh_path, [1007, 1008], "skull", output_path)
    PlanTUS.add_structure_information(output_path + "/skull.surf.gii", "CORTEX_RIGHT")

    PlanTUS.create_avoidance_mask(
        mesh_path, output_path + "/skin.surf.gii", transducer_diameter / 2)

    target_center = PlanTUS.roi_center_of_gravity(target_roi_filepath)

    skin_target_distances = PlanTUS.distance_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)
    PlanTUS.create_metric_from_pseudo_nifti(
        "distances", skin_target_distances, output_path + "/skin.surf.gii")
    PlanTUS.mask_metric(output_path + "/distances_skin.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/distances_skin.func.gii", "CORTEX_LEFT")
    PlanTUS.threshold_metric(output_path + "/distances_skin.func.gii", max_distance)
    PlanTUS.mask_metric(output_path + "/distances_skin_thresholded.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/distances_skin_thresholded.func.gii", "CORTEX_LEFT")

    skin_coordinates, skin_normals = PlanTUS.compute_surface_metrics(
        output_path + "/skin.surf.gii")
    skin_target_vectors = PlanTUS.vectors_between_surface_and_point(
        output_path + "/skin.surf.gii", target_center)
    skin_target_angles = []
    for i in np.arange(len(skin_target_vectors)):
        skin_target_angles.append(math.degrees(
            PlanTUS.angle_between_vectors(skin_target_vectors[i], skin_normals[i])))
    skin_target_angles = np.abs(np.array(skin_target_angles))
    PlanTUS.create_metric_from_pseudo_nifti(
        "angles", skin_target_angles, output_path + "/skin.surf.gii")
    PlanTUS.mask_metric(output_path + "/angles_skin.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/angles_skin.func.gii", "CORTEX_LEFT")

    PlanTUS.stl_from_nii(target_roi_filepath, 0.25)
    skin_target_intersections = PlanTUS.compute_vector_mesh_intersections(
        skin_coordinates, skin_normals,
        output_path + "/" + target_roi_name + "_3Dmodel.stl", 200)

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
    PlanTUS.create_metric_from_pseudo_nifti(
        "target_intersection", skin_target_intersection_values,
        output_path + "/skin.surf.gii")
    PlanTUS.mask_metric(output_path + "/target_intersection_skin.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/target_intersection_skin.func.gii", "CORTEX_LEFT")

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
            skin_skull_angle_list.append(math.degrees(
                PlanTUS.angle_between_vectors(skin_normals[i], skull_normals[closest])))
        except Exception:
            skin_skull_angle_list.append(0)
    skin_skull_angles = np.array(skin_skull_angle_list)
    PlanTUS.create_metric_from_pseudo_nifti(
        "skin_skull_angles", skin_skull_angles, output_path + "/skin.surf.gii")
    PlanTUS.mask_metric(output_path + "/skin_skull_angles_skin.func.gii",
                        output_path + "/avoidance_skin.func.gii")
    PlanTUS.add_structure_information(
        output_path + "/skin_skull_angles_skin.func.gii", "CORTEX_LEFT")

    # ---- Build and open planning scene in wb_view ---------------------------
    scene_variable_names = [
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
    scene_variable_values = [
        'skin.surf.gii',             './skin.surf.gii',
        'skull.surf.gii',            './skull.surf.gii',
        'distances_skin.func.gii',   './distances_skin.func.gii',
        'target_intersection_skin.func.gii', './target_intersection_skin.func.gii',
        'angles_skin.func.gii',      './angles_skin.func.gii',
        'skin_skull_angles_skin.func.gii', './skin_skull_angles_skin.func.gii',
        'distances_skin_thresholded.func.gii', './distances_skin_thresholded.func.gii',
        'T1.nii.gz',                 '../../T1.nii.gz',
        target_mask_filename,        './' + target_mask_filename,
    ]

    scene_path = output_path + "/scene.scene"
    PlanTUS.create_scene(planning_scene_template_filepath, scene_path,
                         scene_variable_names, scene_variable_values)

    print("[prestus_plantus_gui_launcher] Opening Connectome Workbench...")
    print("  Click on a skin vertex in wb_view to select the transducer position.")
    print("  Type 'yes' at the prompt to confirm and generate placement, then")
    print("  close wb_view or press Enter to finish.")

    command = "wb_view -logging FINER " + scene_path
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

    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        shell=True, cwd=output_path, text=True)

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
                    response = input(
                        f"Generate transducer placement for vertex "
                        f"{triangle_number}? (yes/no): ").strip().lower()
                    if response == "yes":
                        print(f"Generating placement for vertex {triangle_number}...")
                        PlanTUS.prepare_acoustic_simulation(
                            triangle_number,
                            output_path,
                            target_roi_filepath,
                            t1_path,
                            max_distance,
                            min_distance,
                            transducer_diameter,
                            max_angle,
                            plane_offset,
                            additional_offset,
                            transducer_surface_model_filepath,
                            focal_distance_list,
                            flhm_list,
                            placement_scene_template_filepath)
                        print(f"Placement written for vertex {triangle_number}.")
                    else:
                        print("No action taken — click another vertex.")
                    process_line = False

    output_thread = threading.Thread(target=read_output)
    output_thread.start()

    process.wait()
    output_thread.join()
    listener.stop()

    if triangle_number is None:
        sys.exit("[prestus_plantus_gui_launcher] wb_view closed without vertex selection. "
                 "No Localite.mat written.")

    print("[prestus_plantus_gui_launcher] Done.")


if __name__ == "__main__":
    main()
