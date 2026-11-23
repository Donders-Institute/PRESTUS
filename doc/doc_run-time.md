# Run Time and Resources

A simulation run can take a few minutes on a GPU or many hours on a CPU, and its memory decides which GPUs can run it. Knowing this before a run helps to choose the resolution, solver, and device, to request job resources, and to plan a study. PRESTUS predicts the time and memory of a run from its config, and records what each run used.

The calculator below shows the expected time of each stage of a typical run, from the same table the Nexus uses. A head uses the reference head, or with a [limited field of view](doc_partial-fov.md) a cylinder around the beam; a free-water grid is sized to the transducer aperture and focal distance. A partial FOV shortens the acoustic simulation from minutes to seconds, but the head is still preprocessed in full, so the CPU stages then take most of the time. The Nexus computes the same prediction from a full config and the subject's own grid.

<div id="load-predictor"><noscript>The load calculator needs JavaScript.</noscript></div>

On this page:

- [Predicted Load](#predicted-load): how the time of each stage and the GPU memory are predicted, before any run.
- [From Earlier Runs](#from-earlier-runs): time and memory of this project's finished runs.

!!! note "At a glance"
    - **When you need it:** Before a study or a change of resolution, solver, or device; when setting job resources on a cluster.
    - **Where:** In the Nexus at the bottom of the Transducer, Medium, Sonication, Solver, and Run tabs; without the Nexus with `python -m prestus_nexus.workload` and `python -m prestus_nexus.estimate`.

## Predicted Load

Before any run, the Nexus predicts the time of each stage the current config runs, the total, and the GPU memory. The prediction appears at the bottom of the Transducer, Medium, Sonication, Solver, and Run tabs and in the Run tab's *Plan*. Without the Nexus: `python -m prestus_nexus.workload path/to/study.yaml`.

- **Grid**: sized to the transducer for free water (`grid.default_dims: 'auto'`), `grid.default_dims`, the acoustic FOV, or for a head the grid of this project's earlier runs. Without earlier runs, a reference head (the SimNIBS example head, 9–14 L depending on orientation) stands in; the subject's own grid replaces it after its first run.
- **Time steps**: as the source setup derives them from the resolution, frequency, sound speeds of the modelled layers, `grid.source_cfl`, and `grid.steady_state_transits`.
- **Acoustic simulation and [free-water baseline](doc_glossary.md#free-water-baseline)**: a fixed start-up time (about 5 s on a GPU) plus grid points with [PML](doc_glossary.md#pml) × time steps × a rate for the device and solver. The baseline runs on a box around the beam with the same time steps. The GPU is set by `hpc.gpu`. For a run on this computer (`platform: matlab`) without `hpc.gpu`, the Nexus reads the local NVIDIA GPU with `nvidia-smi` and marks it as detected. Otherwise the prediction spans all GPU types with a rate. A GPU without its own rate is scaled from GPUs of the same generation by memory bandwidth. A solver without its own rate on that GPU is scaled from `cpp_gpu` (or `cpp_cpu`): `matlab_gpu` takes about 3.5 times as long, and double precision about 1.8 times.
- **Thermal simulation**: thermal grid points × thermal steps × a rate for the device.
- **Other stages** (placement, head model, medium, source, analysis, NIfTI maps): CPU work that takes the same time on any GPU. Each takes a time that depends on the medium and the resolution; head preprocessing is the largest.
- **Segmentation, pseudo-CT, and the HTML report**: timed like the other stages once runs have recorded them. Segmentation (SimNIBS charm) runs once per subject and pseudo-CT generation only for UTE-based skulls, so the calculator lists them apart from the total.
- **Not predicted**: queue wait on a cluster.

Rates and stage times come from this project's runs when they recorded them, else from `nexus/prestus_nexus/throughput.yaml`. Each estimate shows a range from the 10th to the 90th percentile of earlier runs. Across studies, measured acoustic times fall within 0.7–1.4 times the prediction for 80% of runs.

The table was measured on the Donders cluster. On other computers, the prediction is less certain:

- **GPUs without measurements** (for example desktop RTX cards) are scaled by memory bandwidth from measured GPUs of the same generation, which held to about ±25% for the cluster GPUs. A GPU of a generation without any measured GPU has no rate.
- **CPU runs** use the cluster CPU rate; how the time changes with the number of cores is not measured yet.
- **Stages that do not simulate** run on the CPU and take the times measured on cluster nodes; a faster or slower CPU changes them.

Runs on any computer record their times in `run_stats.jsonl` and, with consent, in [telemetry](doc_telemetry.md), from which the table is refitted.

??? developer "Refitting `throughput.yaml`"
    `python -m prestus_nexus.workload fit [--update] <records> ...` fits the table from `run_stats.jsonl` files or telemetry exports (`run_end` events as JSON lines or a JSON array; `sim_stats` holds the sizes and stage times). A fitted entry replaces the entry with the same device and solver (or stage and medium), and other entries stay. The `bandwidth` list (GPU memory bandwidth from NVIDIA specifications, and the generation) is entered by hand and kept. Acoustic entries hold `t0_s` and `ns`, fitted on relative error; a fixed start-up time is fitted only from 6 or more runs spanning a 3-fold range of work, else 5 s. Runs more than 2.5 times off the fit are dropped. Stage entries hold `s_1mm` and `exponent`: seconds = `s_1mm` × (1 mm / resolution)^(3 × `exponent`). `workload.validate` checks the acoustic model leave-one-study-out. After a refit, `python tests/docs/gen_load_data.py` updates the calculator data; CI checks that it matches the table.

## From Earlier Runs

Each run that simulates appends one line to `<path.sim>/run_stats.jsonl`. The Nexus shows an estimate from these lines on the Run tab (*Show deployment, save & resource details*). Without the Nexus:

```
cd <PRESTUS>/nexus
python -m prestus_nexus.estimate path/to/study.yaml [--stats path/to/run_stats.jsonl]
```

The estimate uses earlier runs with the same `simulation.code_type`, `acoustic_backend`, `simulation.medium`, and `grid.resolution_mm`, on the same `hpc.partition` if there are any, else on any partition. It reports the median and maximum of:

- **Acoustic simulation time**, as recorded.
- **Thermal simulation time**: seconds per thermal step × the steps of the current `timing.*` protocol.
- **Peak RAM** of the job.

Segmentation, pseudo-CT, [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) preprocessing and queue wait are not estimated. Estimates appear once runs of the project have finished, and get better as more runs are recorded.

??? note "What `run_stats.jsonl` records"
    Each line holds the acoustic grid the solver ran on (voxels, dimensions, PML, and whether it was axisymmetric) and its time steps, thermal voxels and steps, the acoustic and thermal simulation durations, the duration of each timed stage (`stage_<name>_s`, as in the run log's `⏱` lines), the modelled skull properties (sound speed, density, and attenuation at the source frequency: 10th, 50th, and 90th percentile over all skull voxels, and the median and path length through skull along each transducer's beam to the focus), the run total, solver (`simulation.code_type`, `acoustic_backend`), medium, resolution, `hpc.partition`, node, CPU allocation, GPU model, and peak RAM. Runs that reuse a cached acoustic result, dry runs, and runs with the stage off record nothing for that stage. Peak RAM is the MATLAB process high-water mark, or under SLURM the job cgroup peak, which includes the k-Wave C++ binary and BabelBrain. It is recorded on Linux only. Nothing is sent anywhere; this is separate from [telemetry](doc_telemetry.md).

What a finished cluster job used, including its cost and CO₂e, is in the [post-job resource report](doc_hpc.md#post-job-resource-report).
