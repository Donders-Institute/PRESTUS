# PRESTUS Telemetry

PRESTUS can optionally send anonymous usage statistics to help the developers
understand which features and platforms are in active use. Participation is
entirely voluntary. No data is ever sent without an explicit opt-in.

## Opt-in / opt-out

On each pipeline run PRESTUS prints a short notice and asks whether to enable
telemetry. The prompt repeats every run until you give an explicit answer.
Your decision is stored in:

```
~/.prestus/telemetry.json
```

You can change your decision at any time:

- **Opt out**: set `"opt_in": false` in that file, or run
  `telemetry_setup_reset()` in MATLAB to be asked again on the next run.
- **Opt in**: set `"opt_in": true` in the same file.

### Non-interactive environments (HPC batch jobs)

In non-interactive sessions PRESTUS cannot prompt for input. The full notice
is printed to stdout (visible in job logs) on every run, but **no decision is
recorded and no data is sent**. The message will reappear on every run until
you make an explicit decision.

To opt in from an HPC environment, either:

1. Run PRESTUS once interactively (e.g. on a login node) and answer `y`, or
2. Create `~/.prestus/telemetry.json` manually:
   ```json
   {"opt_in": true, "decided_on": "YYYY-MM-DD", "prestus_ver": ""}
   ```

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

Subject IDs, file paths, or directory names that carries a risk of containing identifiable information is not collected.

## Usage statistics

The charts below are populated live from the anonymised telemetry database.

<style>
.telem-stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:.75rem; margin:1.2rem 0; }
.telem-stat  { background:var(--md-code-bg-color,#f5f5f5); border-radius:8px; padding:1rem; text-align:center; }
.telem-stat .tv { font-size:1.8rem; font-weight:700; line-height:1.1; }
.telem-stat .tl { font-size:.72rem; color:var(--md-default-fg-color--light,#666); text-transform:uppercase; letter-spacing:.04em; margin-top:.2rem; }
.telem-charts { display:grid; grid-template-columns:repeat(auto-fit,minmax(340px,1fr)); gap:1rem; margin:1rem 0; }
.telem-chart  { background:var(--md-code-bg-color,#f5f5f5); border-radius:8px; padding:1rem; }
.telem-chart h4 { font-size:.75rem; font-weight:600; text-transform:uppercase; letter-spacing:.04em; color:var(--md-default-fg-color--light,#666); margin:0 0 .75rem; }
.telem-chart canvas { max-height:200px; }
#telem-error { display:none; color:#c0392b; font-size:.85rem; padding:.5rem 0; }
#telem-loading { color:var(--md-default-fg-color--light,#666); font-size:.85rem; }
</style>

<div id="telem-loading">Loading usage data…</div>
<div id="telem-error"></div>

<div class="telem-stats" id="telem-stats" style="display:none"></div>

<div class="telem-charts" id="telem-charts" style="display:none">
  <div class="telem-chart" style="grid-column:1/-1"><h4>Simulations per month</h4><canvas id="tc-timeline"></canvas></div>
  <div class="telem-chart"><h4>Status</h4><canvas id="tc-status"></canvas></div>
  <div class="telem-chart"><h4>Execution backend</h4><canvas id="tc-codetype"></canvas></div>
</div>

<script>
(function () {
  const URL = 'https://lruelxwhkjibezkgipme.supabase.co';
  const KEY = 'sb_publishable_p18tafIOjIJuGwXb4szLSg_84vF3FL2';
  const PAL = ['#6366f1','#22d3ee','#f59e0b','#10b981','#f43f5e','#a855f7','#14b8a6','#fb923c','#84cc16','#e879f9'];

  function fg()   { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color')        || '#333'; }
  function fgL()  { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--light')  || '#999'; }
  function gridC(){ return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--lightest')|| '#e0e0e0'; }

  function countBy(arr, key) {
    return arr.reduce((a, r) => { const v = r[key] ?? '(unknown)'; a[v] = (a[v]||0)+1; return a; }, {});
  }

  function doughnutChart(id, labels, data) {
    new Chart(document.getElementById(id), {
      type: 'doughnut',
      data: { labels, datasets: [{ data, backgroundColor: PAL.slice(0, labels.length), borderWidth: 2 }] },
      options: { plugins: { legend: { position: 'right', labels: { color: fg(), boxWidth: 12, font: { size: 11 } } } }, cutout: '60%' }
    });
  }

  async function init() {
    const client = supabase.createClient(URL, KEY);
    const { data, error } = await client
      .from('events')
      .select('received_at,status,prestus_ver,prestus_hash,duration_s,uuid,code_type')
      .order('received_at', { ascending: false })
      .limit(5000);

    document.getElementById('telem-loading').style.display = 'none';

    if (error) {
      const el = document.getElementById('telem-error');
      el.textContent = 'Could not load telemetry data: ' + error.message;
      el.style.display = 'block';
      return;
    }

    // stat cards
    const total    = data.length;
    const durs     = data.map(r => r.duration_s).filter(Boolean);
    const avgDur   = durs.length ? (durs.reduce((a,b)=>a+b,0)/durs.length/60).toFixed(1)+'min' : '—';
    const versions = [...new Set(data.map(r => r.prestus_ver).filter(Boolean))];
    const users    = new Set(data.map(r => r.uuid).filter(Boolean)).size;

    const statsEl = document.getElementById('telem-stats');
    statsEl.style.display = 'grid';
    [['Run count', total], ['Unique users', users], ['Avg duration', avgDur], ['Versions', versions.length]]
      .forEach(([label, val]) => {
        statsEl.insertAdjacentHTML('beforeend',
          `<div class="telem-stat"><div class="tv">${val}</div><div class="tl">${label}</div></div>`);
      });

    document.getElementById('telem-charts').style.display = 'grid';

    // ── Simulations per month, stacked by PRESTUS version ──────────
    const allVersions = [...new Set(data.map(r => r.prestus_ver).filter(Boolean))].sort();
    const monthSet = new Set();
    data.forEach(r => { if (r.received_at) monthSet.add(r.received_at.slice(0, 7)); });
    const months = [...monthSet].sort();

    // count per (month, version)
    const mvCount = {};
    data.forEach(r => {
      if (!r.received_at || !r.prestus_ver) return;
      const key = r.received_at.slice(0, 7) + '|' + r.prestus_ver;
      mvCount[key] = (mvCount[key] || 0) + 1;
    });

    new Chart(document.getElementById('tc-timeline'), {
      type: 'bar',
      data: {
        labels: months,
        datasets: allVersions.map((ver, i) => ({
          label: ver,
          data: months.map(m => mvCount[m + '|' + ver] || 0),
          backgroundColor: PAL[i % PAL.length],
          borderRadius: 2,
          borderSkipped: false,
        }))
      },
      options: {
        plugins: { legend: { position: 'right', labels: { color: fg(), boxWidth: 12, font: { size: 11 } } } },
        scales: {
          x: { stacked: true, ticks: { color: fgL(), maxRotation: 45, font: { size: 10 } }, grid: { display: false } },
          y: { stacked: true, ticks: { color: fgL() }, grid: { color: gridC() }, beginAtZero: true }
        }
      }
    });

    // ── Status donut ───────────────────────────────────────────────
    const sc = countBy(data, 'status');
    doughnutChart('tc-status', Object.keys(sc), Object.values(sc));

    // ── Execution backend (code_type) ──────────────────────────────
    const cc = countBy(data, 'code_type');
    const ck = ['matlab_cpu','matlab_gpu','cpp_cpu','cpp_gpu'].filter(k => cc[k]);
    doughnutChart('tc-codetype', ck, ck.map(k => cc[k]));

  }

  if (typeof Chart !== 'undefined' && typeof supabase !== 'undefined') {
    init();
  } else {
    window.addEventListener('load', init);
  }
})();
</script>
