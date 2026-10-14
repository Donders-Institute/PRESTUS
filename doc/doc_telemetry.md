# PRESTUS Telemetry

PRESTUS can optionally send anonymous usage statistics to help the developers
understand which features and platforms are in active use. Participation is
entirely voluntary. No data is ever sent without an explicit opt-in.

## Opt-in / opt-out

On each pipeline run PRESTUS prints a short notice and asks whether to enable
telemetry. How that prompt is shown depends on the session:

- **Interactive console session**: a `y`/`n` prompt that blocks until you answer.
- **GUI session** (`prestus_gui`): a modal **Allow / Decline / Ask me later**
  dialog shown once at startup, before the "Run Simulation" button (which may
  dispatch to a background worker with no console) can be pressed.
- **`matlab -batch` run locally with a display attached**: the same
  Allow / Decline / Ask me later dialog as the GUI, since there is no console
  to answer `y`/`n` on.
- **Headless session** (deployed app, or a SLURM/PBS job — see
  [below](#non-interactive-environments-hpc-batch-jobs)): no prompt of any
  kind can be shown.

The prompt repeats every run until you give an explicit answer. Your decision
is stored in:

```
~/.prestus/telemetry.json
```

You can change your decision at any time:

- **Opt out**: set `"opt_in": false` in that file, or run
  `telemetry_setup_reset()` in MATLAB to be asked again on the next run.
- **Opt in**: set `"opt_in": true` in the same file.

### Non-interactive environments (HPC batch jobs)

Deployed applications and jobs that run under SLURM (`SLURM_JOB_ID`) or
PBS/Torque (`PBS_JOBID`) are always headless: no prompt of any kind is
possible, even if a figure could technically be drawn. `matlab -batch` run
locally (not under SLURM/PBS) falls back to the same silent notice only when
no display is attached; if a display is attached it gets the dialog described
above. The full notice is printed to stdout (visible in job logs) on every
run, but **no decision is recorded and no data is sent**. The message will
reappear on every run until you make an explicit decision.

To opt in from an HPC environment, either:

1. Run PRESTUS once interactively (e.g. on a login node) and answer `y`,
2. Run `telemetry_set_consent(true)` (or `false` to opt out) in MATLAB, or
3. Create `~/.prestus/telemetry.json` manually:
   ```json
   {"opt_in": true, "decided_on": "YYYY-MM-DD", "prestus_ver": ""}
   ```

### Technical vs. deployment runs

Test and CI runs that exercise `prestus_pipeline` (e.g. via
`run_all_tests('water'|'head'|'acoustic'|'all'|'physics'|...)`) are tagged
`run_type: technical` automatically, so their pass/fail rate can be tracked
separately from real deployment usage on the
[usage statistics](doc_usage_statistics.md) page. Set
`PRESTUS_RUN_TYPE=technical` yourself when running PRESTUS purely to validate
an install or a pipeline change; leave it unset for real runs (default:
`deployment`).

## What is collected

Each pipeline run emits one event: `run_end` (on clean exit) or `run_error`
(if the pipeline throws an unhandled exception). The table below lists every
field that may be sent.

### System resources

| Field | Description |
|---|---|
| `os_version` | OS version string (macOS: `sw_vers -productVersion`; Linux: kernel release; Windows: `ver` output). Returns `'unknown'` on error. |
| `cpu_model` | CPU model name (e.g. `Apple M2 Pro`). Returns `'unknown'` on error. |
| `n_cpu_cores` | Number of logical CPU cores via `feature('numcores')`. |
| `ram_gb` | Total physical RAM in GB (rounded). Returns `NaN` on error. |
| `gpu_model` | GPU device name from `gpuDevice()`. Returns `'no_pct'` when the Parallel Computing Toolbox is not licensed; `'none'` when no GPU is found; otherwise the device name string. |
| `toolboxes` | Cell array of installed MATLAB toolbox names (MATLAB itself excluded). |

All system resource fields are collected inside `try/catch` blocks and fall back to `'unknown'` or `NaN` on any error, so they never cause the telemetry payload to fail.

### Identity

| Field | Description |
|---|---|
| `event` | Event type: `run_end` or `run_error` |
| `run_type` | `technical` (test/CI run) or `deployment` (default). See [Technical vs. deployment runs](#technical-vs-deployment-runs). |
| `timestamp_utc` | POSIX timestamp (seconds since epoch, UTC) |
| `uuid` | Random UUID generated locally on first run; stored in `~/.prestus/uuid.txt`. Not linked to your identity or machine. |
| `run_id` | Random UUID generated per pipeline run, so the events of one run can be told apart. Not linked to your identity or machine. |
| `prestus_ver` | Version string of the running PRESTUS checkout (from `git describe` or the `VERSION` file, e.g. `v0.5.1`) |
| `prestus_hash` | Short git commit hash of the running PRESTUS checkout (e.g. `628edd0`) |
| `matlab_ver` | MATLAB release string (e.g. `R2024a`) |
| `kwave_ver` | k-Wave toolbox version string as reported by `getComputerInfo` |
| `platform` | CPU architecture string from `computer('arch')` (e.g. `maci64`, `glnxa64`) |

### Execution environment

| Field | Description |
|---|---|
| `sim_platform` | Execution platform: `matlab`, `slurm`, `qsub`, etc. |
| `code_type` | k-Wave code variant: `matlab_cpu`, `matlab_gpu`, `cpp_cpu`, `cpp_gpu` |
| `precision` | Floating-point precision: `single` or `double` |
| `hpc_name` | HPC cluster name if running on a known cluster; `unknown` otherwise |
| `use_gpu` | Boolean — whether a GPU device is configured |

### Simulation configuration

| Field | Description |
|---|---|
| `medium` | Medium type: `water` or `layered` |
| `layers` | List of tissue layer names (e.g. `["skin","skull","brain"]`) |
| `pct_enabled` | Boolean — whether pseudo-CT is enabled |
| `pct_density` | Name of the density mapping used for pCT (string, not values) |
| `pct_speed` | Name of the sound-speed mapping used for pCT |
| `pct_atten` | Name of the attenuation mapping used for pCT |
| `modules_enabled` | Names of the `modules.*` flags set to 1 (e.g. `["run_acoustic_sims","run_heating_sims"]`) |
| `pipeline_mode` | Pipeline mode of the run: `single`, `async` (asynchronous multi-transducer), `multi_isppa` or `sequential` |

### Transducer

| Field | Description |
|---|---|
| `transducer_type` | Transducer geometry type (e.g. `focused_bowl`) |
| `transducer_name` | Transducer model name from the config |
| `freq_hz` | Nominal operating frequency in Hz |
| `n_transducer_elements` | Number of transducer elements (`elem_n` from the transducer config) |

### Outcome (run_end / run_error only)

| Field | Description |
|---|---|
| `duration_s` | Total pipeline wall-clock time in seconds |
| `status` | `"success"` or `"error"` |
| `error_id` | MATLAB error identifier string (e.g. `MATLAB:badsubscript`) — **no message text** is transmitted |

Subject IDs, file paths, or directory names that carries a risk of containing identifiable information is not collected.

## Usage statistics

See the [Usage Statistics](doc_usage_statistics.md) page for live charts populated from the anonymised telemetry database.
