# PRESTUS Anonymous Usage Statistics

PRESTUS can optionally send anonymous usage statistics to help the developers
understand which features and platforms are in active use. Participation is
entirely voluntary and is requested once on the first pipeline run.

## Opt-in / opt-out

On first run, PRESTUS prints a short prompt and asks whether to enable
telemetry. The answer is stored in:

```
~/.prestus/telemetry.json
```

You can change your decision at any time:

- **Opt out**: set `"opt_in": false` in that file, or run
  `telemetry_setup_reset()` in MATLAB to be asked again on the next run.
- **Opt in**: set `"opt_in": true` in the same file.

In non-interactive environments (HPC batch jobs) telemetry is silently
disabled by default.

## What is collected

Each pipeline run emits up to two events: `run_start` (before any module
runs) and `run_end` (on clean exit) or `run_error` (if the pipeline throws
an unhandled exception). The table below lists every field that may be sent.

### Identity

| Field | Description |
|---|---|
| `event` | Event type: `run_start`, `run_end`, or `run_error` |
| `timestamp_utc` | POSIX timestamp (seconds since epoch, UTC) |
| `uuid` | Random UUID generated locally on first run; stored in `~/.prestus/uuid.txt`. Not linked to your identity or machine. |
| `prestus_ver` | Short git commit hash of the running PRESTUS checkout |
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
| `modules_enabled` | List of module names that are switched on (e.g. `["acoustic_sim","thermal_sim"]`) |

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

## What is never collected

- Subject IDs, file paths, or directory names
- Coordinates, focus targets, or geometry parameters
- Acoustic pressure values or simulation outputs
- Tissue property values (only whether they differ from defaults)
- Error message text or stack traces
- IP addresses or hostnames
- Any free-text field from the configuration

## Where data goes

Events are sent via HTTPS POST to a Supabase Edge Function operated by the
Donders Institute and stored in a private PostgreSQL database. Transmission
uses a fire-and-forget pattern with a 4-second timeout; a failed request is
silently discarded and never retries.

## Data retention

Data is retained for a minimum of two years to support longitudinal version
adoption analysis. No personally identifiable information is stored; records
cannot be linked back to individuals or institutions.

## Public aggregate statistics

A read-only view of aggregate statistics is published at:

> https://lruelxwhkjibezkgipme.supabase.co/functions/v1/public-stats

This endpoint returns a JSON summary of platform distributions, MATLAB version
adoption, module usage rates, and run counts over time. It is updated
continuously from the same database that receives your events.
