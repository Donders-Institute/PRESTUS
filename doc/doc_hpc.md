# Running on a Cluster (HPC)

A head simulation takes minutes to hours on a GPU, and a study runs many of them (subjects, targets, settings). The most efficient way to use PRESTUS is to run these simulations as parallel cluster jobs.

On this page:

- [Where a Run Executes](#where-a-run-executes): what each `platform` value does, and where job logs go.
- [Job Resources](#job-resources): wall time, memory, and the other `hpc.*` keys.
- [GPU Support](#gpu-support): requesting GPUs and partitions, with benchmark timings.
- [Post-Job Resource Report](#post-job-resource-report): what a finished SLURM job used.

!!! note "At a glance"
    - **When you need it:** For head simulations without a local GPU, and whenever you run many simulations.
    - **Before you start:** Access to a cluster with [SLURM](https://slurm.schedmd.com/) or PBS/qsub, and PRESTUS and its dependencies installed there.

PRESTUS submits cluster jobs itself. You call the same `prestus_pipeline_start(parameters)` as for a local run, and `platform` decides where it runs:

```yaml
platform: slurm          # auto | matlab | slurm | qsub
hpc:
  timelimit: '04:00:00'
  memorylimit: 20        # GB
```

Many subjects or conditions are submitted together from a [study manifest](doc_group_study-setup.md#study-setup) with `prestus_study_start`. The [PRESTUS Nexus](doc_nexus.md) can submit to a cluster, also over SSH from a laptop, and follows the jobs and their logs.

On the Donders cluster, see [Donders HPC](doc_hpc_donders.md) for interactive sessions and your own job scripts.

## Where a Run Executes

| `platform` | Effect |
|---|---|
| `auto` (default) | SLURM if `sbatch` is on the PATH, else qsub if `qsub` is; **an error if neither is found**. Set `matlab` on a workstation |
| `matlab` | Run in the current [MATLAB](https://www.mathworks.com/products/matlab.html) session (recommended for debugging) |
| `slurm` | Submit a SLURM job (`sbatch`) |
| `qsub` | Submit a PBS/Torque job (`qsub`) |

Job logs go to `log_hpc/` in the subject's output folder (`sub-XXX_slurm_output_<id>.log`; `sub-XXX_qsub_output_<id>.log` for qsub). MATLAB returns as soon as the job is submitted, unless `hpc.wait_for_job` is set. A small dependent job then writes a [resource report](#post-job-resource-report).

Some modes split one run into several chained jobs that pass results through caches: async multi-transducer runs, multi-intensity sweeps, and uncertainty runs ([How PRESTUS works](doc_architecture.md#run-execution-and-chained-jobs)). Post-hoc free-water simulations (`modules.run_posthoc_water_sims`) are submitted as separate jobs.

To submit, PRESTUS saves the parameters to a temporary file, writes a small launcher that adds the PRESTUS paths and calls `prestus_pipeline`, and submits it.

## Job Resources

Each job requests 4 hours of wall time (`hpc.timelimit`) and 20 GB of memory (`hpc.memorylimit`); post-hoc free-water jobs request 5 hours and 40 GB. Raise these for fine grids or heating runs. GPU type, partition, and reservation are described under [GPU support](#gpu-support), node-local disk under [Node-local scratch](#node-local-scratch), and `hpc.name: 'snellius'` selects the [Snellius](https://www.surf.nl/en/services/compute/snellius-the-national-supercomputer) presets. All `hpc` keys, with their defaults, are listed under [HPC parameters](doc_parameters_run.md#high-performance-computing-hpc).

### Node-local scratch

The C++ k-Wave codes (`cpp_cpu`, `cpp_gpu`) write an input and an output `.h5` file per simulation, often 5–15 GB together. Inside a job they go to `$TMPDIR` when the scheduler sets it, else to `<output>/cache/`. To reserve node-local disk for them and choose where they go, set `hpc.tmp_gb`:

```yaml
hpc:
  tmp_gb: 50          # adds '#SBATCH --tmp=50G' to the job script; 0 = off (default)
  scratch_dir: auto   # or a path, e.g. '/scratch/${USER}/${SLURM_JOB_NAME}' (created if missing)
```

With `scratch_dir: auto`, PRESTUS takes the first existing folder of `$PRESTUS_SCRATCH`, `$SLURM_TMPDIR`, `/scratch/$USER/$SLURM_JOB_NAME`, and `/scratch/$USER/$SLURM_JOB_ID`. If none exists but `/scratch/$USER` does, it creates `/scratch/$USER/$SLURM_JOB_NAME` and removes it again at the end if it is empty. The last fallback is `$TMPDIR`. On a cluster whose scratch is not under `/scratch`, set `PRESTUS_SCRATCH_BASE` before MATLAB starts. The job log names the folder in a `Node-local scratch:` line, with the free space.

PRESTUS works in `prestus_<jobid>/sub-NNN<affix>/` inside that folder and removes it when the run ends or fails, including `.h5` files a failed k-Wave run left behind. A job killed by the scheduler leaves the folder to the cluster's own cleanup. `cache/` always stays in the output folder, because later stages and runs read it.

If no node can provide `--tmp=<tmp_gb>G`, the job stays pending: lower `hpc.tmp_gb` or set it to `0`.

All keys are listed in [Parameters](doc_parameters_run.md#high-performance-computing-hpc). Which simulations need a GPU and which settings are fastest is described in [k-Wave](doc_kwave.md). Predicted run time and memory, and the time and memory earlier runs needed, are on [Run time and resources](doc_run-time.md).

## GPU Support

TUS simulations are accelerated by GPUs, but requesting GPUs can lead to longer wait times (on the Donders cluster, the concurrent GPU limit per user is 4). To reduce wait times, run acoustic simulations first (to confirm targeting; they need less RAM) and then run thermal simulations with the final protocol. Avoid blanket simulations (e.g. circling through all participants with all permutations), especially for thermal simulations.

A job requests a GPU when a solver it runs computes on one:
- k-Wave with `simulation.code_type` `matlab_gpu` or `cpp_gpu`, for the acoustic stage, k-Wave heating, or time reversal;
- jwave, unless `jwave.use_gpu` is `false`;
- BabelBrain, unless `babelbrain.computing_backend` is `cpu`.

`hpc.partition` and `hpc.gpu` override the partition and GPU type.

PRESTUS is a MATLAB toolbox and only supports *Nvidia GPUs*. When Nvidia GPUs are digitally partitioned, MATLAB R2024+ can fail to identify the assigned GPU. Submitted SLURM jobs therefore load `matlab/R2023b`, whichever MATLAB version submits them (`MATLAB/2024b` with `hpc.name: 'snellius'`). PBS jobs load `matlab/R2022b`, because the PBS nodes support CUDA 11.2 only.

`hpc.gpu` requests a specific GPU type and is usually not needed. On the Donders cluster, set `hpc.partition: 'gpu40g'` (GPUs with more than 40 GB vRAM) for most thermal simulations; `hpc.reservation` selects a reserved queue, if you have one. See [HPC parameters](doc_parameters_run.md#high-performance-computing-hpc).

??? note "Benchmark data, Nvidia GPUs"
    These are potentially unrepresentative benchmarks run on a 256 × 216 × 192 mm grid, with minor variations depending on transducer placement.

    **Acoustic simulations**

    | GPU | Memory Used | Duration | Notes |
    |---|---|---|---|
    | A100 80 GB | 12 GB | 14 mins | Slightly smaller grid size |
    | A100 80 GB (partitioned 2×40 GB) | 12 GB | 25 mins | |
    | A100 40 GB | 12 GB | 19 mins | |
    | A16 16 GB | 12 GB | 145 mins | |
    | P100 16 GB | 12 GB | 43 mins | Compiled; ~1.5× the L40S duration |
    | L40S 47 GB | 14 GB | 28 mins | Compiled |

    **Heating simulations**

    | GPU | Memory Used | Duration | Notes |
    |---|---|---|---|
    | A100 80 GB | ?? GB | ~12 s/trial (400 trials: ~90 mins) | Slightly smaller grid size |
    | A100 40 GB | ?? GB | ~17 s/trial (400 trials: ~115 mins) | |
    | A16 16 GB | Out of RAM | ??? | |
    | P100 16 GB | Out of RAM | ??? | |
    | L40S 47 GB | ?? GB | ~4 s/trial (400 trials: ~45 mins) | |

## Post-Job Resource Report

After a SLURM job finishes, `hpc_job_report` writes `log_hpc/sub-XXX_job<ID>_report.txt` with wall time, CPU/GPU allocation, peak memory, a CO₂e estimate and, where billing rates are known, the cost. The report runs from `hpc_wait_for_completion` when `hpc.wait_for_job` is set, and otherwise as a small dependent job. `memory_efficiency` and the `--mem` suggestion use `peak_mem_gb`.

??? note "Where the peak memory comes from"
    Peak memory comes from two sources:

    - `max_rss_gb`: SLURM's `sacct MaxRSS`, the maximum over job steps. SLURM samples it (every 30 s by default), so it can miss short spikes.
    - The job logs: on Linux, BabelBrain runs print memory lines to the console. MATLAB writes `[memory]` lines (from `log_memory_usage`) before and after the BabelBrain launch. The BabelBrain Python process writes `[memory_bb]` lines whenever its peak grows by 0.5 GB or more, plus one at exit. The report reads these lines from the pipeline log (`io.log_file`, fixed at submission) and from `log_hpc/<subj>_slurm_output_<ID>.log`:
        - `log_cgroup_peak_gb`: peak usage of the job's memory cgroup (MATLAB plus BabelBrain). This is the figure the out-of-memory killer acts on, and it includes page cache.
        - `bb_peak_gb` and `matlab_peak_gb`: the peak resident memory of each process.

    `peak_mem_gb` is the larger of `log_cgroup_peak_gb` and `max_rss_gb`, and `peak_mem_source` names the source it came from. Both are lower bounds: the log lines stop at BabelBrain exit, so a later stage can go above them.

    For local runs, each BabelBrain run also prints one summary line to the pipeline log, e.g.
    `[memory_summary] BabelBrain: BabelBrain peak=7.2 GB | MATLAB peak=3.1 GB | cgroup peak=10.5 GB, limit=64.0 GB`.
