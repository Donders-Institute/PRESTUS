# Tutorial: A Head Simulation

!!! warning "Draft"
    This tutorial has not yet been run end to end on the current PRESTUS version. Runtimes and the values under *What you should see* are placeholders (marked **TODO**) until a test run fills them in. It stays in `examples/` until then; move it to `doc/` and add it to the nav in `mkdocs.yaml` once it has been run.

This tutorial runs a complete head simulation on *Ernie*, the example head that comes with SimNIBS: segmentation, placement from MNI coordinates, the acoustic simulation, the thermal simulation, and the report. It builds on [Your first simulation](doc_getting-started.md); terms are explained in the [Glossary](doc_glossary.md).

!!! note "At a glance"
    - **What:** one subject, one target, a transducer from the library, acoustic and thermal simulation.
    - **Time:** segmentation **TODO** h; acoustic **TODO** min on a GPU; thermal **TODO** min.
    - **Before you start:** PRESTUS with k-Wave, SimNIBS 4 ([Installation](doc_installation.md)), and a GPU or a cluster ([Running on a cluster](doc_hpc.md)).

The files used here are in the toolbox: `examples/tutorial_ernie.yaml` (the study config) and `examples/tutorial_ernie.m` (the three steps).

---

## 0. Get the Ernie data

Download the SimNIBS example dataset (see the [SimNIBS documentation](https://simnibs.github.io/simnibs/build/html/dataset.html)) and unpack it. You need the T1 and T2 images, `ernie_T1.nii.gz` and `ernie_T2.nii.gz`.

Open `examples/tutorial_ernie.yaml` and set the four paths:

| Key | Set to |
|---|---|
| `path.anat` | the folder with `ernie_T1.nii.gz` and `ernie_T2.nii.gz` |
| `path.seg` | an empty folder for the segmentation |
| `path.sim` | an empty folder for the simulation outputs |
| `startup.simnibs_bin_path` | the `bin` folder of your SimNIBS environment |

The rest of the file describes the study:

- **Transducer:** a Sonic Concepts CTX-500 with its driving system, from the [transducer library](doc_transducer_library.md), set to focus 50 mm from its exit plane. Results are scaled to a free-water intensity of 30 W/cm².
- **Placement:** `mode: mni`. The transducer is placed on the left side of the head and aimed at a target given in MNI coordinates ([MNI placement](doc_placement_mni.md)). PRESTUS converts both points to Ernie's head and moves the transducer onto the scalp.
- **Protocol:** 20 ms pulses at 5 Hz for 20 s, followed by 20 s of cooling ([Thermal](doc_simulations-thermal.md)).

## 1. Segmentation

Run the first section of `examples/tutorial_ernie.m`. It runs SimNIBS `charm` on the T1 and T2 and stops (`modules.segmentation_only: 1`).

**What you should see:** a folder `m2m_sub-001` in `path.seg`. Its `final_tissues.nii.gz` labels skin, skull, brain and the other tissues. Open it over the T1 in a NIfTI viewer and check that the skull is complete and the brain is not cut off. **TODO:** runtime.

The segmentation is reused by all later runs; it is only redone with `io.overwrite_simnibs: 1`.

## 2. Placement and acoustic simulation

Run the second section. PRESTUS places the transducer, builds the simulation grid from the segmentation, runs a short free-water baseline of the transducer, and then the acoustic simulation through the head. Thermal simulation is still off (the default), so you can check the acoustic result first.

**What you should see** in `path.sim/sub-001/`:

- `sub-001_placement_mni.html`: the transducer on the left scalp, its beam pointing at the target. **TODO:** screenshot.
- `sub-001_layered_report.html`: the report. The intensity map should show the focus near the target.
- `sub-001_layered.csv`: **TODO:** typical values of `Isppa_brain`, `Ipa_target`, `real_focal_distance_mm` and `MI_tc` from the test run.

Go through steps 1 and 2 of [Reading your results](doc_results.md): is the focus at the target, and is the intensity in the brain plausible?

## 3. Thermal simulation

Run the third section. It switches on `modules.run_heating_sims` and runs the pipeline again, now simulating heating over the protocol.

**What you should see:** the report now contains temperature and thermal-dose maps, and the safety summary includes temperature rise, maximum temperature and CEM43. **TODO:** values from the test run, and whether any card is amber or red.

Read step 3 of [Reading your results](doc_results.md) for how to interpret the safety summary.

## Next steps

- Your own data: copy `tutorial_ernie.yaml` into your project and change the paths and `t1_pattern` ([Configuring a study](doc_configuration.md)).
- Another target or transducer position: [Placement](doc_placement.md).
- Subject-specific skull properties from a UTE scan: [pseudo-CT](doc_pseudoCT.md).
- Many subjects: [Studies & group analysis](doc_group.md).
