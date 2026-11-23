# Segmentation

![Pipeline steps; this page belongs to 1 Segmentation](assets/figures/fig_stage-strip-segmentation_diagram.svg){ .stage-strip }

The skull reflects, attenuates, and refracts ultrasound, so the simulation needs to know where skull, brain, and soft tissue lie. [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) `charm` labels these tissues in the MRI; PRESTUS then builds the [head model](doc_head-model.md) from the labels.

On this page:

- [Running the Segmentation](#running-the-segmentation): running SimNIBS `charm` and reusing its result.
- [SimNIBS Version and charm Settings](#simnibs-version-and-charm-settings): choosing the install and changing charm's settings.
- [Defaced Copy for Sharing](#defaced-copy-for-sharing): a copy of the segmentation without the face.

!!! note "At a glance"
    - **When you need it:** Every `layered` simulation. Free-water simulations skip it.
    - **Before you start:** A T1w image (optionally a T2w) and a SimNIBS install ([Installation](doc_installation.md#install-simnibs)).

## Running the Segmentation

For a `layered` medium, the first run segments the T1 (and the T2, when `path.t2_pattern` is set) with SimNIBS `charm` into `<path.seg>/m2m_sub-XXX/`. This takes an hour or more. Later runs reuse `final_tissues.nii.gz` from that folder unless `io.overwrite_simnibs: 1`. The [air-cavity check](doc_head-model.md#air-cavity-check) also reads charm's detailed labelling, `segmentation/labeling.nii.gz`, because `final_tissues.nii.gz` has no air label: charm merges the sinuses and other internal air into the background (0). Keep the `segmentation/` folder when you tidy up m2m folders.

- `modules.segmentation_only: 1` segments and stops, e.g. to segment all subjects before simulating them.
- `segmentation.use_qform: 1` when `charm` reports a qform/sform mismatch; `segmentation.debug: 1` passes `--debug` to `charm`.
- `segmentation.noneck: 1` passes `--noneck` for a scan without a neck.

```yaml
startup:
  simnibs_bin_path: /path/to/SimNIBS-4.x/bin   # the SimNIBS install to use
modules:
  segmentation_only: 1                         # segment, then stop
segmentation:
  noneck: 0                                    # 1 for a scan without a neck
```

## SimNIBS Version and charm Settings

`startup.simnibs_bin_path` selects the SimNIBS install, and with it the version; PRESTUS needs SimNIBS 4 or later. While they are empty, `hpc.ld_library_path` (the install's libraries for SimNIBS tools) and `placement.plantus.env_path` (the Python PlanTUS runs in) follow it; a value set there overrides it. To change charm's settings, copy the install's `charm.ini` (in its `simnibs` package folder), edit the copy, and set `segmentation.settings_file` to it; keep it outside the m2m folder. Report the settings you changed.

An existing segmentation is reused even when another version or other settings are selected. PRESTUS then warns, comparing the version in `charm_log.html` and the settings in `settings.ini` of the m2m folder with the selection. Set `io.overwrite_simnibs: 1` to segment again, which deletes the m2m folder first, or point `path.seg` to a separate folder per version to keep both. What was built from the old segmentation is detected on reuse: grid, medium, acoustic, and head-model caches, placement runs, and pseudo-CTs record the segmentation's checksum ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)), and a defaced copy is made again.

In the Nexus, the Segmentation card (Medium tab) lists the SimNIBS installs on this computer, shows the version and the existing segmentation, opens `charm_report.html`, and has the re-run switch. When runs are submitted to a cluster over SSH (Run tab), *Scan cluster* lists the installs there: conda environments, `~/SimNIBS*`, `/opt/simnibs*`, and the bin folders of the `simnibs` environment modules. The Nexus keeps the result per host, so it can check a cluster path that does not exist on this computer.

## Defaced Copy for Sharing

`segmentation.deface: 1` writes a copy of the m2m folder, `m2m_sub-XXX_defaced/`, with the face removed from the anatomical images by [pydeface](https://github.com/poldracklab/pydeface) ([install it](doc_installation.md#mri-defacing) in the SimNIBS environment). Tissue labels, surfaces, warps, and meshes are copied unchanged, and the original folder is not modified. The copy is made once and again after the segmentation changes. Figures are defaced separately, with `io.deface_plots` (on by default).

Installing SimNIBS: [Installation](doc_installation.md#install-simnibs). Errors: [Troubleshooting](doc_troubleshooting.md#segmentation-simnibs-charm).
