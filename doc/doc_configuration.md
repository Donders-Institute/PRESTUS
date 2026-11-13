# Configuring a Study

All PRESTUS settings live in [YAML](https://yaml.org/) config files, and the config decides what a run simulates: the head, the transducer and its placement, and the protocol. You never write a complete config: PRESTUS starts from its defaults, and your files contain only what differs. Every parameter is listed in [Parameters](doc_parameters.md).

On this page:

- [How a Config Is Assembled](#how-a-config-is-assembled): the four layers `load_parameters` merges, and how to override settings per run.
- [Project Config Folder](#project-config-folder): setting up a project folder that pins its defaults.
- [Example Study Config](#example-study-config): an example config for head simulations.
- [Configs From Earlier PRESTUS Versions](#configs-from-earlier-prestus-versions): how renamed and retired keys are handled.
- [Editing Configs With a GUI](#editing-configs-with-a-gui): the Nexus, which writes the same files.
- [Prompts for Missing Tool Paths](#prompts-for-missing-tool-paths): when PRESTUS asks for missing tool paths.

!!! note "At a glance"
    - **When you need it:** To set up a study, or to change settings for one subject or run.
    - **Before you start:** PRESTUS installed ([Installation](doc_installation.md)).

## How a Config Is Assembled

`load_parameters` builds one `parameters` struct from up to four layers, each merged onto the previous one:

| Layer | Source | Notes |
|---|---|---|
| 1. Toolbox defaults | `config/config_default.yaml` in the PRESTUS folder | Always loaded first. Never edit it for a study. |
| 2. Project defaults | a copy of `config_default.yaml` in your project config folder (optional) | Pins the defaults your project started with (see below). |
| 3. Study config(s) | the YAML file(s) or struct you pass in | Only the values your study changes. |
| 4. Overrides | a trailing `{'dotted.key', value; ...}` cell | Per-run changes, e.g. the subject or one row of a study manifest. |

```matlab
parameters = load_parameters();                                       % defaults only
parameters = load_parameters('config_projectX.yaml', '/path/to/project/config');
parameters = load_parameters('config_a.yaml', dir_a, 'config_b.yaml', dir_b);   % merged left to right
parameters = load_parameters('config_projectX.yaml', cfg_dir, ...
                             {'subject_id', 3; 'io.output_affix', 'M1'});         % plus overrides
```

The merge is deep: a study file that sets only `grid.resolution_mm` keeps every other `grid` setting at its default.

??? note "Merge rules and validation"
    Lists of transducers (`transducer(2)`, …) keep their own fields. One exception: `medium_properties.<tissue>.absorption_fraction` and `absorption_coeff` exclude each other, so a layer that sets one empties the other (see [Heat deposition](doc_simulations-thermal.md#heat-deposition)).

    Before merging, keys from earlier PRESTUS versions are converted to their current names ([see below](#configs-from-earlier-prestus-versions)). After merging, `load_parameters` validates the result: it resolves `transducer.serial` against the [transducer library](doc_transducer_library.md), checks the timing protocol and the backend combination, and turns `io.output_affix` into a `_desc-<label>` file-name entity. Keys PRESTUS does not know produce the warning `prestus:unknownConfigKey`, with a suggestion when it looks like a typo.

You can also change the struct in [MATLAB](https://www.mathworks.com/products/matlab.html) after loading (`parameters.subject_id = 2;`), as in [Your first simulation](doc_getting-started.md). The exception is the transducer library keys (`transducer.serial`, `combo`, `library_path`, and the `focal_distance_ep` they resolve phases for): `load_parameters` looks up geometry, phases and drive from them, so set them in a file or as an override rather than in the struct afterwards.

## Project Config Folder

Keep your configs in a project folder outside the PRESTUS toolbox, and start it with `prestus_config_init`:

```matlab
prestus_config_init('/path/to/project/config', project_name='projectX');
% optionally open the GUI to review the values:
prestus_config_init('/path/to/project/config', project_name='projectX', open_gui=true);
```

This writes two files:

- `config_default.yaml`, a copy of the toolbox defaults (layer 2). It keeps your project's baseline fixed when PRESTUS changes its defaults later, which matters for reproducibility.
- `config_projectX.yaml`, an empty study config (layer 3) for your own values.

Load them with `load_parameters('config_projectX.yaml', '/path/to/project/config')`. PRESTUS looks for the project `config_default.yaml` next to the study config, else in the folder you pass, else in the current folder (with the warning `prestus:cwdDefault`).

??? note "Updating a project copy of the defaults"
    When PRESTUS later adds new settings, a project copy made earlier does not contain them. They then take their toolbox default, and `load_parameters` warns (`prestus:staleProjectDefault`). `prestus_config_upgrade(folder)` replaces the copy with the current toolbox default, keeps a dated backup, and can first write the values your copy pinned differently into the study config, so the simulation settings stay the same. The Nexus offers the same under **Tools ▸ Update project defaults…**.

## Example Study Config

A study config for head simulations typically sets the paths, the transducer, the placement and the protocol:

```yaml
# config_projectX.yaml
path:
  anat: /data/projectX/anat            # T1 (and T2/UTE) images
  t1_pattern: 'sub-%1$03d/anat/sub-%1$03d_T1w.nii.gz'   # %03d = subject_id
  seg: /data/projectX/simnibs          # SimNIBS m2m_sub-XXX folders
  sim: /data/projectX/sim              # simulation outputs
startup:
  simnibs_bin_path: /home/user/.conda/envs/simnibs_v4.6.0/bin

simulation:
  medium: layered
  code_type: matlab_gpu

transducer:
  serial: CTX_500_026
  combo: {ds_serial: SC_203_035}
  focal_distance_ep: 55
  target_isppa_wcm2: 30

placement:
  mode: mni
  mni:
    trans_pos_mm: [-52, -8, 58]
    focus_pos_mm: [-29, -6, 7]

timing:
  pd: 0.02
  pri: 0.2
  ptd: 20
  ptri: 20
  ptrd: 20
modules:
  run_heating_sims: 1
```

On a workstation without a job scheduler, also set `platform: matlab`: the default `auto` only detects [SLURM](https://slurm.schedmd.com/) or qsub and stops with an error when neither is found. The subject is usually not in the file but set per run (`parameters.subject_id = 3`, an override, or a [study manifest](doc_group.md#study-setup)).

??? note "Where to look next"
    | Topic | Keys | Page |
    |---|---|---|
    | Paths and input images | `path.*`, `startup.*` | [Parameters](doc_parameters.md#data-paths-path) |
    | Medium and tissue properties | `simulation.medium`, `layers`, `medium_properties`, `pct.*` | [Medium Setup](doc_medium.md), [pseudo-CT](doc_pseudoCT.md) |
    | Transducer | `transducer.*` | [Transducer](doc_transducer.md), [Transducer Library](doc_transducer_library.md) |
    | Position and target | `placement.*`, `transducer.trans_pos`/`focus_pos` | [Placement](doc_placement.md) |
    | Protocol and heating | `timing.*`, `thermal.*` | [Thermal](doc_simulations-thermal.md) |
    | Which stages run | `modules.*` | [Modules](doc_modules.md) |
    | Where and how it runs | `platform`, `hpc.*`, `simulation.code_type`, `acoustic_backend` | [HPC](doc_hpc.md), [Backends](doc_backend.md), [k-Wave](doc_kwave.md) |
    | Output files | `io.*` | [Outputs](doc_outputs.md) |

## Configs From Earlier PRESTUS Versions

`load_parameters` reads configs written for PRESTUS v0.5.0 or later. It converts each renamed key to its current name, in every config layer and in the overrides, and lists the conversions in one warning (`prestus:renamedConfigKey`); keys without a replacement are dropped. Old keys without an exact counterpart stop the run (see below), as does a layer that sets an old key and its current name to different values (`prestus:configKeyConflict`). Update the config with the names from that list; the [Nexus](doc_nexus.md) marks the same keys when you open the file.

??? note "Old keys without a current name"
    A few old keys have no exact counterpart. These stop the run (`prestus:legacyConfigKey`) with a message naming the key to set instead:

    | Old key | Set instead |
    |---|---|
    | `smooth_window` [voxels] | `headmodel.smooth_fwhm_mm` [mm] |
    | `expected_focal_distance_mm` | `transducer.focal_distance_ep` (from the exit plane) or `transducer.focal_distance_bowl` (from the bowl apex) |
    | `pseudoCT_variant: marquet` | `pct.mapping_density`, `pct.mapping_soundspeed`, `pct.mapping_attenuation` |
    | `calibration.weights: 0` | `calibration.opt_weights` (now 0 = uniform weighting) |
    | `subject_subfolder: 0`, custom `seg_labels` | not supported: outputs always go to `sub-XXX` folders, and the SimNIBS charm labels are fixed |

    A config that sets an old key and its current name to different values also stops the run. Conversion covers names only: defaults that changed since your config was written (for example `transducer.target_isppa_wcm2`) apply unless the config sets them, as listed in the [changelog](CHANGELOG.md). `calibration.skip_front_peak_mm` kept its name but is now measured from the exit plane rather than the bowl apex.

Two front ends write the same YAML files:
- the [PRESTUS Nexus](doc_nexus.md) (experimental), a desktop app that also checks the config as you type and shows a placement viewer, and
- the [MATLAB GUI](doc_gui.md) (`prestus_gui`, deprecated).
## Editing Configs With a GUI

Both write only the values that differ from the defaults, so their files load with `load_parameters` like hand-written ones.

## Prompts for Missing Tool Paths

Missing tool paths (`startup.simnibs_bin_path`, `hpc.ld_library_path`) only prompt for input in an interactive MATLAB desktop session with `simulation.interactive: 1`. Otherwise PRESTUS warns and continues, so batch and cluster runs never wait on a dialog.
