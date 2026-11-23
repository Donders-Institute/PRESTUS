# PRESTUS Telemetry

PRESTUS can send anonymous usage statistics, one record per pipeline event, so the developers can see which features and platforms are in use, improve run-time estimates, and publish simulation statistics across sites. Participation is voluntary, you choose what to share in three levels, and nothing is sent without your choice. This page shows how to make or change that choice, how batch and cluster runs handle it, and every field a record may contain.

On this page:

- [Choosing What to Share](#choosing-what-to-share): the prompt, the three levels, how to record or change your choice, the `PRESTUS_TELEMETRY` override, and batch or cluster runs.
- [Collected Data](#collected-data): every field of a record, grouped by level.

!!! note "At a glance"
    - **When you need it:** To decide what to share, to change your choice, or to silence the prompt on a cluster.
    - **Before you start:** Nothing; your choice lives in `~/.prestus/telemetry.json`.

## Choosing What to Share

Usage statistics come in three levels. You can share any combination:

| Level | Name in settings | What it adds |
|---|---|---|
| Usage and errors | `basic` | Which features, platforms, and versions are used, and whether runs succeed (error identifier and failing stage, never the message text). |
| Performance | `performance` | Hardware, requested resources, grid size, stage times, and peak memory. These calibrate the [predicted load](doc_run-time.md#predicted-load). |
| Summary outcomes (anonymous) | `outcomes` | Eight rounded summary statistics: brain intensity and its ratio to free water, mechanical index, targeting error, brain and skull temperature rise, and the duty cycle and duration they depend on; plus the median modelled skull properties. Never maps, coordinates, or target names. |

Every record also carries a small envelope whichever levels you share: a random installation ID, a random run ID, the PRESTUS version, the run type, the run status, and the levels themselves.

The prompt shows **Share all (recommended)** with the three levels below it, all ticked. Untick the levels you do not want to share; unticking all of them shares nothing. How the prompt appears depends on the session:

- **Interactive console session**: the levels are listed; press Enter to share all, type the numbers of the levels to share (e.g. `1 2`), or `0` to share none.
- **MATLAB GUIs and the [PRESTUS Nexus](doc_nexus.md)**: a dialog with the same boxes, and **Save** or **Ask me later**. In the Nexus, **Tools ▸ Usage statistics…** opens it at any time, and the status bar shows what you share.
- **Non-interactive session** (`matlab -batch`, deployed app, [SLURM](https://slurm.schedmd.com/)/PBS job, or a session with neither desktop nor display; see [below](#batch-and-cluster-runs)): no prompt is shown and the run never waits for an answer.

The prompt repeats every run until you make a choice. A choice made before the levels existed keeps sharing what it covered (usage, performance, and the skull medians); an interactive session offers the levels once. To change your choice:

- **Share levels**: run `telemetry_set_consent(true)` in [MATLAB](https://www.mathworks.com/products/matlab.html) for all levels, or name them, e.g. `telemetry_set_consent({'basic', 'performance'})`.
- **Share nothing**: run `telemetry_set_consent(false)`.
- **Undecided again**: `telemetry_setup_reset()` deletes the file, so the next run asks again.

### Per-Process Override (`PRESTUS_TELEMETRY`)

The `PRESTUS_TELEMETRY` environment variable sets the choice for one MATLAB process, without a prompt and without writing `~/.prestus/telemetry.json`:

??? note "Values and precedence"
    | Value (case-insensitive) | Effect for this process |
    |---|---|
    | `off`, `0`, `false` | Nothing sent, no prompt. |
    | `on`, `1`, `true` | All levels shared if no choice is recorded, no prompt. A recorded choice still applies. The random ID in `~/.prestus/uuid.txt` is created on the first send. |
    | levels, comma-separated (e.g. `basic,performance`) | Those levels shared if no choice is recorded, no prompt. |
    | unset or empty | No override: the recorded choice applies, or the prompt/notice described above. |
    | anything else | Warning `prestus:telemetry:invalidEnv` (once per value per MATLAB session), then treated as unset. |

    `PRESTUS_TELEMETRY=off` takes precedence over a recorded choice: it sends nothing even after you chose to share. Sharing through the variable never overrides a recorded choice; it only applies when no choice is recorded. The variable never changes the recorded choice; use `telemetry_set_consent` for a lasting one. Set the variable before MATLAB starts (`PRESTUS_TELEMETRY=off matlab -batch "..."`, or `export PRESTUS_TELEMETRY=off` in a job script) or with `setenv('PRESTUS_TELEMETRY', 'off')` in the session.

### Batch and Cluster Runs

PRESTUS never prompts, and never waits for an answer, in these sessions:

??? note "How these sessions are detected"
    - `matlab -batch` runs (detected with `batchStartupOptionUsed`, R2019a and later), whether or not a display is attached: cron jobs, CI runners, and scripts started from a terminal;
    - deployed applications and jobs under SLURM (`SLURM_JOB_ID`) or PBS/Torque (`PBS_JOBID`);
    - sessions with neither the MATLAB desktop nor a display (e.g. `matlab -nodisplay -r ...` from cron). An interactive `matlab -nodisplay` terminal session is also treated this way; record a choice there with `telemetry_set_consent`.

The full notice is printed to stdout (visible in job logs) on every run, but **no choice is recorded and no data is sent** unless a choice was already recorded or `PRESTUS_TELEMETRY` is set. The notice reappears on every run until you record a choice or set `PRESTUS_TELEMETRY`.

To choose from an HPC environment, either:

1. Run PRESTUS once interactively (e.g. on a login node) and answer the prompt,
2. Run `telemetry_set_consent(true)`, a list of levels, or `false` in MATLAB,
3. Set `PRESTUS_TELEMETRY` in the job script, for that job only, or
4. Create `~/.prestus/telemetry.json` by hand:

    ```json
    {"opt_in": true, "levels": ["basic", "performance", "outcomes"], "decided_on": "YYYY-MM-DD", "prestus_ver": "", "notice": 3}
    ```

### Technical, Development, and Deployment Runs

Test and CI runs that exercise `prestus_pipeline` are tagged `run_type: technical` automatically, so their pass/fail rate can be tracked separately from real deployment usage on the [usage statistics](doc_usage_statistics.md) page. Set `PRESTUS_RUN_TYPE=technical` yourself when running PRESTUS purely to validate an install or a pipeline change, and `PRESTUS_RUN_TYPE=development` (in your shell or `startup.m`) on a machine where you develop PRESTUS; leave it unset for real runs (default: `deployment`). Jobs that PRESTUS submits to a cluster report the run type of the session that submitted them, including one set with `setenv` in MATLAB. Independently of the run type, `prestus_branch` and `prestus_dirty` mark runs from a checkout off `main` or with local changes.

## Collected Data

Each pipeline run emits a `run_start` record and then `run_end` (on clean exit) or `run_error` (if the pipeline throws an unhandled exception), with the same `run_id`. A run with a start but no end record was stopped from outside, e.g. by a scheduler's time or memory limit. The tables below list every field that may be sent, by level. A field listed under several levels is sent when any of them is shared.

**Never stored in the usage statistics:** subject IDs, file paths, directory names, coordinates, target names, free text, IP addresses, host names, or any map or image of the simulation. Free-text settings are sent as classes only (`transducer_model`, `hpc_partition`). Like any web service, the hosting provider ([Supabase](https://supabase.com/docs/guides/observability/logs)) keeps a request log with the sender's IP address and approximate location for 1 day; it is not linked to the statistics.

### Every Record

??? note "Fields"
    | Field | Description |
    |---|---|
    | `event` | Event type: `run_start`, `run_end`, or `run_error` |
    | `telemetry_schema` | Version of the record format (currently `3`) |
    | `consent_levels`, `consent_notice` | The shared levels, and the version of the notice the choice was made under |
    | `run_type` | `technical` (test/CI run), `development`, or `deployment` (default). See [Technical, Development, and Deployment Runs](#technical-development-and-deployment-runs). |
    | `uuid` | Random UUID generated locally on first send; stored in `~/.prestus/uuid.txt`. Not linked to your identity or machine. |
    | `run_id` | Random UUID generated per pipeline run, so the records of one run can be told apart |
    | `prestus_ver`, `prestus_hash` | Version string (from `git describe` or the `VERSION` file, e.g. `v0.5.1`) and short [git](https://git-scm.com/) commit hash of the running PRESTUS |
    | `status` | `"success"` or `"error"` (on `run_end` and `run_error`) |

### Usage and Errors (`basic`)

??? note "Fields"
    | Field | Description |
    |---|---|
    | `prestus_branch` | Branch class of the PRESTUS checkout: `main`, `development`, `other`, or `unknown` (never the branch name) |
    | `prestus_dirty` | `yes` when tracked PRESTUS files differ from the commit, else `no` or `unknown` |
    | `matlab_ver`, `platform`, `kwave_ver` | MATLAB release (e.g. `R2024a`), CPU architecture from `computer('arch')`, and the [k-Wave](https://github.com/ucl-bug/k-wave) version from `getComputerInfo` |
    | `sim_platform` | Execution platform: `matlab`, `slurm`, `qsub`, etc. |
    | `code_type`, `acoustic_backend`, `thermal_backend` | k-Wave code variant (`matlab_cpu`, `matlab_gpu`, `cpp_cpu`, `cpp_gpu`), acoustic solver (`kwave`, `kwave_elastic`, `jwave`, or `babelbrain`), and thermal solver |
    | `medium`, `pct_enabled` | Medium type (`water`, `layered`, or `phantom`) and whether pseudo-CT is enabled |
    | `pct_density`, `pct_speed`, `pct_atten` | Names of the pseudo-CT mappings (names, not values) |
    | `modules_enabled`, `pipeline_mode` | Names of the `modules.*` flags set to 1, and the pipeline mode: `single`, `async`, `multi_isppa`, or `sequential` |
    | `acoustic_fov`, `headmodel_reuse` | Booleans: a limited acoustic field of view is set; `headmodel.reuse_across_placements` |
    | `transducer_type`, `freq_hz`, `n_transducer_elements` | Geometry type (e.g. `focused_bowl`), nominal frequency [Hz], and element count of the first transducer |
    | `transducer_model` | Equipment model of the first transducer when its `serial` or `name` matches a shipped entry (e.g. `CTX_500`); `custom` otherwise. Never the name or serial. |
    | `n_transducers`, `transducer_coupling` | Number of transducers and how they are combined (`simulation.transducer_coupling`) |
    | `duration_s` | Total pipeline wall-clock time in seconds |
    | `error_id`, `failed_stage` | MATLAB error identifier (e.g. `MATLAB:badsubscript`; **no message text**) and the stage the run was in when it failed (e.g. `acoustic`) |

### Performance (`performance`)

??? note "Fields"
    | Field | Description |
    |---|---|
    | `cpu_model`, `ram_gb`, `gpu_model` | CPU model name, total RAM [GB], and GPU device name (`no_pct` without the [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html), `none` without a GPU) |
    | `precision` | Floating-point precision: `single` or `double` |
    | `hpc_timelimit_h`, `hpc_memory_gb` | Requested wall time [h] and memory [GB] of the job (`hpc.timelimit`, `hpc.memorylimit`) |
    | `hpc_partition` | Class of `hpc.partition`: `gpu` (the name contains "gpu"), `other`, or empty when unset; never the partition name |
    | `n_cpu_alloc` | CPU cores of the SLURM allocation (`SLURM_CPUS_PER_TASK`); `null` outside SLURM |
    | `resolution_mm`, `fit_alpha_power` | Grid resolution [mm], and whether attenuation is refitted to one power law |
    | `source_cfl`, `source_ppw`, `steady_state_transits`, `sim_time_us` | Time-step and duration settings of the acoustic simulation (`grid.*`; `null` when unset) |
    | `sim_stats.acoustic_voxels`, `sim_stats.acoustic_dims`, `sim_stats.acoustic_pml` | Grid points and dimensions of the acoustic grid the solver ran on, and its PML thickness in grid points |
    | `sim_stats.axisymmetric` | 1 for an axisymmetric (2D) acoustic grid |
    | `sim_stats.acoustic_nt`, `sim_stats.acoustic_s` | Acoustic time steps and simulation time in seconds |
    | `sim_stats.thermal_voxels`, `sim_stats.thermal_steps`, `sim_stats.thermal_s` | Thermal grid points, time steps, and simulation time in seconds |
    | `sim_stats.stage_<name>_s` | Duration in seconds of each timed pipeline stage (e.g. `stage_preproc_s`) |
    | `sim_stats.peak_ram_gb`, `sim_stats.n_cpu`, `sim_stats.gpu_mem_gb` | Peak RAM of the job [GB] (Linux only), CPU cores of the allocation or machine, and GPU memory [GB] (a MIG slice reports its own) |
    | `sim_stats.queue_wait_s` | Seconds the SLURM job waited in the queue |

    A stage that did not simulate (reused cache, stage off) leaves its fields out. A failed run sends the stages it finished.

### Summary Outcomes (`outcomes`)

Rounded to two significant figures (temperatures to 0.1 °C, skull medians to three significant figures), for the first transducer and the first field a run analyses. Values a run does not compute, such as brain values in a water run, are left out. These records also carry `medium`, `pct_enabled`, the solvers, the transducer model, type and frequency, and `resolution_mm`, so the numbers can be read without the other levels.

??? note "Fields"
    | Field | Description |
    |---|---|
    | `outcomes.isppa_brain_wcm2` | Peak pulse-average intensity in the brain [W/cm²] |
    | `outcomes.brain_transmission` | Brain I<sub>SPPA</sub> divided by the free-water I<sub>SPPA</sub> of the same drive: the share of intensity that reaches the brain |
    | `outcomes.mi_brain` | Peak mechanical index in the brain |
    | `outcomes.focal_error_mm` | Distance from the intended focus to the centre of the −3 dB focal region [mm] |
    | `outcomes.riset_brain_c`, `outcomes.riset_skull_c` | Peak temperature rise above baseline in brain and skull [°C] |
    | `outcomes.duty_cycle`, `outcomes.sonication_s` | Overall duty cycle (pulse duty × train duty, 0–1) and total sonication time [s] from `timing.*` |
    | `sim_stats.skull_c_median`, `sim_stats.skull_rho_median`, `sim_stats.skull_atten_median` | Median sound speed [m/s], density [kg/m³], and attenuation at the source frequency [dB/cm] of the modelled skull, over all skull voxels |
    | `sim_stats.skull_path_c_median`, `sim_stats.skull_path_rho_median`, `sim_stats.skull_path_atten_median` | The same medians over the skull on the beam path |

## Usage Statistics

See the [Usage Statistics](doc_usage_statistics.md) page for live charts populated from the anonymised telemetry database.
