# Running on a Cluster (HPC)

!!! note "At a glance"
    - **What:** Running PRESTUS as jobs on a computing cluster.
    - **When you need it:** For head simulations without a local GPU, and whenever you run many simulations.
    - **Before you start:** Access to a cluster with SLURM or PBS/qsub, and PRESTUS and its dependencies installed there.

PRESTUS submits cluster jobs itself. You call the same `prestus_pipeline_start(parameters)` as for a local run, and `platform` decides where it runs:

```yaml
platform: slurm          # auto | matlab | slurm | qsub
hpc:
  timelimit: '04:00:00'
  memorylimit: 20        # GB
```

The most efficient way to use PRESTUS is to run many simulations (subjects, targets, settings) as parallel cluster jobs. Many subjects or conditions are submitted together from a [study manifest](doc_group.md#study-setup) with `prestus_study_start`. The [PRESTUS Nexus](doc_nexus.md) can submit to a cluster, also over SSH from a laptop, and follows the jobs and their logs.

## Where a run executes

| `platform` | Effect |
|---|---|
| `auto` (default) | SLURM if `sbatch` is on the PATH, else qsub if `qsub` is; **an error if neither is found**. Set `matlab` on a workstation |
| `matlab` | Run in the current MATLAB session (recommended for debugging) |
| `slurm` | Submit a SLURM job (`sbatch`) |
| `qsub` | Submit a PBS/Torque job (`qsub`) |

Job logs go to `log_hpc/` in the subject's output folder (`sub-XXX_slurm_output_<id>.log`). MATLAB returns as soon as the job is submitted, unless `hpc.wait_for_job` is set. A small dependent job then writes a [resource report](#post-job-resource-report).

Some modes split one run into several chained jobs that pass results through caches: async multi-transducer runs, multi-intensity sweeps and uncertainty runs ([How PRESTUS works](doc_architecture.md#where-a-run-executes)). Post-hoc free-water simulations (`modules.run_posthoc_water_sims`) are submitted as separate jobs.

??? note "Details: how a job is submitted"
    PRESTUS saves the parameters to a temporary file, writes a small launcher that adds the PRESTUS paths and calls `prestus_pipeline`, and submits it.

## Job resources

| Key | Default | Meaning |
|---|---|---|
| `hpc.timelimit` | `'04:00:00'` | Wall time per job |
| `hpc.memorylimit` | `20` | Memory per job [GB] |
| `hpc.gpu`, `hpc.partition`, `hpc.reservation` | `''` | GPU type, partition, reservation (see [GPU support](#gpu-support)) |
| `hpc.posthoc_water_timelimit`, `hpc.posthoc_water_memorylimit` | `'05:00:00'`, `40` | Resources of the post-hoc free-water jobs |
| `hpc.job_prefix` | `'PS'` | Prefix of job names |
| `hpc.wait_for_job`, `hpc.max_wait_checks` | `false`, `540` | Block MATLAB until the job ends (one check every 20 s) |
| `hpc.ld_library_path` | `''` | Extra shared-library path for SimNIBS, if `charm` reports `undefined symbol` errors |
| `hpc.name` | `'default'` | `'snellius'` selects the Snellius partition presets (`snellius.*`) |

All keys are listed in [Parameters](doc_parameters.md#high-performance-computing-hpc). Which simulations need a GPU and which settings are fastest is described in [k-Wave](doc_kwave.md). The time and memory earlier runs needed are summarised under [Resource estimates](#resource-estimates-from-earlier-runs).

## GPU support

TUS simulations are accelerated by GPUs, but requesting GPUs can lead to longer wait times (on the Donders cluster, the concurrent GPU limit per user is 4). To reduce wait times, run acoustic simulations first (to confirm targeting; they need less RAM) and then run thermal simulations with the final protocol. Avoid blanket simulations (e.g. circling through all participants with all permutations), especially for thermal simulations.

PRESTUS is a MATLAB toolbox and only supports *Nvidia GPUs*. When Nvidia GPUs are digitally partitioned, MATLAB R2024+ can fail to identify the assigned GPU. For SLURM jobs, MATLAB R2023b is deployed by default.

| Field | Default | Explanation |
|---|---|---|
| `hpc.gpu` | `''` | Request a specific GPU type (e.g. `"nvidia_a100-sxm4-40gb:1"`); usually not needed. `scontrol show nodes \| egrep -o gres/gpu:.*=[0-9] \| egrep -o 'nvidia_.*=' \| sort \| uniq \| sed 's/=//'` lists the available types. |
| `hpc.partition` | `''` | Queue partition. The Donders cluster offers `gpu40g` (GPUs with > 40 GB vRAM), which should be used for most thermal simulations. |
| `hpc.reservation` | `''` | Reserved queue, if you have one. |

??? note "Details: benchmark data, Nvidia GPUs"
    These are potentially unrepresentative benchmarks run on a 256 × 216 × 192 mm grid, with minor variations depending on transducer placement.

    **Acoustic simulations**

    | GPU | Memory Used | Duration | Notes |
    |---|---|---|---|
    | A100 80 GB | 12 GB | 14 mins | Slightly smaller grid size |
    | A100 80 GB (partitioned 2×40 GB) | 12 GB | 25 mins | |
    | A100 40 GB | 12 GB | 19 mins | |
    | A16 16 GB | 12 GB | 145 mins | |
    | P100 16 GB | 12 GB | 43 mins | Compiled, ~L40s |
    | L40S 47 GB | 14 GB | 28 mins | Compiled |

    **Heating simulations**

    | GPU | Memory Used | Duration | Notes |
    |---|---|---|---|
    | A100 80 GB | ?? GB | ~12 s/trial (400 trials: ~90 mins) | Slightly smaller grid size |
    | A100 40 GB | ?? GB | ~17 s/trial (400 trials: ~115 mins) | |
    | A16 16 GB | Out of RAM | ??? | |
    | P100 16 GB | Out of RAM | ??? | |
    | L40S 47 GB | ?? GB | ~4 s/trial (400 trials: ~45 mins) | |

## Resource estimates from earlier runs

Each run that simulates appends one line to `<path.sim>/run_stats.jsonl`. The Nexus shows an estimate from these lines on the Run tab (*Show deployment, save & resource details*). Without the Nexus:

```
cd <PRESTUS>/nexus
python -m prestus_nexus.estimate path/to/study.yaml [--stats path/to/run_stats.jsonl]
```

The estimate uses earlier runs with the same `simulation.code_type`, `acoustic_backend`, `simulation.medium` and `grid.resolution_mm`, on the same `hpc.partition` if there are any, else on any partition. It reports the median and maximum of:

- **Acoustic simulation time**, as recorded.
- **Thermal simulation time**: seconds per thermal step × the steps of the current `timing.*` protocol.
- **Peak RAM** of the job.

Segmentation, pseudo-CT, BabelBrain preprocessing and queue wait are not estimated. Estimates appear once runs of the project have finished, and get better as more runs are recorded.

??? note "Details: what `run_stats.jsonl` records"
    Each line holds the acoustic grid voxels and time steps, thermal voxels and steps, the acoustic and thermal simulation durations, the run total, solver (`simulation.code_type`, `acoustic_backend`), medium, resolution, `hpc.partition`, node, GPU model and peak RAM. Runs that reuse a cached acoustic result, dry runs, and runs with the stage off record nothing for that stage. Peak RAM is the MATLAB process high-water mark, or under SLURM the job cgroup peak, which includes the k-Wave C++ binary and BabelBrain. It is recorded on Linux only. Nothing is sent anywhere; this is separate from [telemetry](doc_telemetry.md).

## Post-job resource report

After a SLURM job finishes, `hpc_job_report` writes `log_hpc/sub-XXX_job<ID>_report.txt` with wall time, CPU/GPU allocation, peak memory, a CO₂e estimate and, where billing rates are known, the cost. The report runs from `hpc_wait_for_completion` when `hpc.wait_for_job` is set, and otherwise as a small dependent job. `memory_efficiency` and the `--mem` suggestion use `peak_mem_gb`.

??? note "Details: where the peak memory comes from"
    Peak memory comes from two sources:

    - `max_rss_gb`: SLURM's `sacct MaxRSS`, the maximum over job steps. SLURM samples it (every 30 s by default), so it can miss short spikes.
    - The job logs: on Linux, BabelBrain runs print memory lines to the console. MATLAB writes `[memory]` lines (from `log_memory_usage`) before and after the BabelBrain launch. The BabelBrain Python process writes `[memory_bb]` lines whenever its peak grows by 0.5 GB or more, plus one at exit. The report reads these lines from the pipeline log (`io.log_file`, fixed at submission) and from `log_hpc/<subj>_slurm_output_<ID>.log`:
        - `log_cgroup_peak_gb`: peak usage of the job's memory cgroup (MATLAB plus BabelBrain). This is the figure the out-of-memory killer acts on, and it includes page cache.
        - `bb_peak_gb` and `matlab_peak_gb`: the peak resident memory of each process.

    `peak_mem_gb` is the larger of `log_cgroup_peak_gb` and `max_rss_gb`, and `peak_mem_source` names the source it came from. Both are lower bounds: the log lines stop at BabelBrain exit, so a later stage can go above them.

    For local runs, each BabelBrain run also prints one summary line to the pipeline log, e.g.
    `[memory_summary] BabelBrain: BabelBrain peak=7.2 GB | MATLAB peak=3.1 GB | cgroup peak=10.5 GB, limit=64.0 GB`.

## Donders HPC

The Donders cluster runs SLURM. Current resource usage is shown at https://grafana.dccn.nl/ (*intranet required*); the [HPC wiki](https://hpc.dccn.nl/) (*intranet required*) has the full documentation. The workflows below start an interactive session from which you submit PRESTUS jobs.

??? note "Details: Donders HPC, SLURM workflow"
    1. Log in to an access node (`mentat005`, `mentat006`, `mentat007`): `ssh abcxyz@mentat005.dccn.nl`
    2. Start the VNC manager (`vncmanager - 2`) and connect via [VNC](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC). TigerVNC: enter ip [e.g. `mentat007.dccn.nl:12`]
    3. Load SLURM: `module load slurm`
    4. Start an interactive job:
        - Without MATLAB GUI: `srun --mem=8gb --time=01:00:00 --x11 -p interactive --pty bash -i`
        - With MATLAB GUI: `srun --partition=gpu --gres=gpu:1 --mem=8G --time=01:00:00 --x11 --pty /bin/bash -i`
    5. Start MATLAB: `module load matlab/R2024a` then `matlab`. SLURM supports CUDA 12.2; recent MATLAB versions (up to R2024) should be supported, see [this issue](https://github.com/Donders-Institute/PRESTUS/issues/50).
    6. Run PRESTUS with `platform: slurm` (or `auto`). Check job status with `squeue`. For PBS→SLURM command migration see [this documentation](https://hpc.dccn.nl/docs/cluster_howto/compute_slurm.html#migrating-from-torque-pbs-to-slurm).

??? note "Details: Donders HPC, PBS workflow"
    1. Log in to an access node (`mentat001`, `mentat002`, `mentat003`, `mentat004`): `ssh abcxyz@mentat004.dccn.nl`
    2. Start the VNC manager (`vncmanager - 2`) and connect via [VNC](https://intranet.donders.ru.nl/index.php?id=vnc00&no_cache=1&sword_list%5B%5D=VNC). TigerVNC: enter ip [e.g. `mentat004.dccn.nl:12`]
    3. Load QSUB: `module load qsub`
    4. Start an interactive job:

        ```
        qsub -I -l 'nodes=1:gpus=1,feature=cuda,walltime=05:00:00,mem=24gb,reqattr=cudacap>=5.0'
        ```

    5. Start MATLAB: `module load matlab/R2022b` then `matlab`. PBS only supports CUDA 11.2, which is dropped from R2023 on; see [this issue](https://github.com/Donders-Institute/PRESTUS/issues/50).
    6. Run PRESTUS with `platform: qsub` (or `auto`). Check job status with `qstat`.
