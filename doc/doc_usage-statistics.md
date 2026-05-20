# Usage Statistics

The charts below are populated live from the anonymised telemetry database. Data is only collected from users who have opted in. See [Telemetry](doc_telemetry.md) for details on what is collected and how to opt in or out.

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
  <div class="telem-chart" style="grid-column:1/-1"><h4>Completed simulations per month</h4><canvas id="tc-timeline"></canvas></div>
  <div class="telem-chart"><h4>Execution platform</h4><canvas id="tc-platform"></canvas></div>
  <div class="telem-chart"><h4>Execution backend</h4><canvas id="tc-codetype"></canvas></div>
  <div class="telem-chart"><h4>Simulation medium</h4><canvas id="tc-medium"></canvas></div>
  <div class="telem-chart"><h4>Mean run duration by medium</h4><canvas id="tc-duration"></canvas></div>
</div>

<script>
(function () {
  const URL = 'https://lruelxwhkjibezkgipme.supabase.co';
  const KEY = 'sb_publishable_p18tafIOjIJuGwXb4szLSg_84vF3FL2';
  const PAL = ['#6366f1','#22d3ee','#f59e0b','#10b981','#f43f5e','#a855f7','#14b8a6','#fb923c','#84cc16','#e879f9'];

  // Maps known git hashes (7-char) to { tag, n } from `git describe --tags --long`.
  // Format shown: "tag - n - hash". Add new entries when hashes appear in telemetry.
  const HASH_INFO = {
    '9eb180b': { tag: 'v0.1.0', n: 0   },
    '2d7c014': { tag: 'v0.2.0', n: 0   },
    '05146df': { tag: 'v0.2.1', n: 0   },
    'd3886f7': { tag: 'v0.3.0', n: 0   },
    '24debe5': { tag: 'v0.4.0', n: 0   },
    '2f56459': { tag: 'v0.4.1', n: 0   },
    'a907aa3': { tag: 'v0.4.2', n: 0   },
    'bc67583': { tag: 'v0.5.0', n: 0   },
    'b7cce9f': { tag: 'v0.5.0', n: 389 },
    'a9305dc': { tag: 'v0.6.0', n: 0   },
    '75badde': { tag: 'v0.6.0', n: 15  },
    'ad10be4': { tag: 'v0.6.0', n: 16  },
  };

  function resolveVersion(row) {
    const hash = (row.prestus_hash || '').slice(0, 7);
    const ver  = row.prestus_ver  || '';
    if (hash && HASH_INFO[hash]) {
      const { tag, n } = HASH_INFO[hash];
      return `${tag} - ${n} - ${hash}`;
    }
    // Recent runs: prestus_ver is a tag but hash not yet in lookup table
    if (ver.startsWith('v')) {
      return hash ? `${ver} - ? - ${hash}` : ver;
    }
    // Older runs: prestus_ver was a bare hash
    const h = (hash || ver || '').slice(0, 7);
    return h || 'unknown';
  }

  function classifyMedium(row) {
    const med = row.medium || 'unknown';
    if (med === 'layered' && row.pct_enabled) return 'layered + pCT';
    return med;
  }

  function fg()   { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color')         || '#333'; }
  function fgL()  { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--light')  || '#999'; }
  function gridC(){ return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--lightest')|| '#e0e0e0'; }

  function countBy(arr, fn) {
    return arr.reduce((a, r) => { const v = fn(r); a[v] = (a[v]||0)+1; return a; }, {});
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
      .select('received_at,event,status,duration_s,uuid,sim_platform,code_type,medium,pct_enabled,prestus_ver,prestus_hash')
      .in('event', ['run_end', 'run_error'])
      .order('received_at', { ascending: false })
      .limit(5000);

    document.getElementById('telem-loading').style.display = 'none';

    if (error) {
      const el = document.getElementById('telem-error');
      el.textContent = 'Could not load telemetry data: ' + error.message;
      el.style.display = 'block';
      return;
    }

    // ── Stat cards ─────────────────────────────────────────────────
    const total     = data.length;
    const successes = data.filter(r => r.status === 'success');
    const successPct = total ? Math.round(successes.length / total * 100) + '%' : '—';
    const users     = new Set(data.map(r => r.uuid).filter(Boolean)).size;
    const versions  = new Set(data.map(resolveVersion).filter(v => v !== 'unknown')).size;

    const statsEl = document.getElementById('telem-stats');
    statsEl.style.display = 'grid';
    [
      ['Completed runs', total],
      ['Success rate',   successPct],
      ['Versions',       versions],
      ['Unique users',   users],
    ].forEach(([label, val]) => {
      statsEl.insertAdjacentHTML('beforeend',
        `<div class="telem-stat"><div class="tv">${val}</div><div class="tl">${label}</div></div>`);
    });

    document.getElementById('telem-charts').style.display = 'grid';

    // ── Simulations per month, stacked by resolved version tag ─────
    const allVersions = [...new Set(data.map(resolveVersion))].sort();
    const monthSet = new Set(data.map(r => r.received_at?.slice(0, 7)).filter(Boolean));
    const months   = [...monthSet].sort();

    const mvCount = {};
    data.forEach(r => {
      const key = r.received_at?.slice(0, 7) + '|' + resolveVersion(r);
      if (key) mvCount[key] = (mvCount[key] || 0) + 1;
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

    // ── Execution platform ─────────────────────────────────────────
    const pc = countBy(data, r => r.sim_platform ?? '(unknown)');
    doughnutChart('tc-platform', Object.keys(pc), Object.values(pc));

    // ── Execution backend (code_type) ──────────────────────────────
    const cc = countBy(data, r => r.code_type ?? '(unknown)');
    const ck = ['matlab_cpu','matlab_gpu','cpp_cpu','cpp_gpu'].filter(k => cc[k]);
    const unknownCC = Object.keys(cc).filter(k => !['matlab_cpu','matlab_gpu','cpp_cpu','cpp_gpu'].includes(k));
    doughnutChart('tc-codetype',
      [...ck, ...unknownCC],
      [...ck.map(k => cc[k]), ...unknownCC.map(k => cc[k])]
    );

    // ── Simulation medium (layered+pCT as its own category) ────────
    const mc = countBy(data, classifyMedium);
    doughnutChart('tc-medium', Object.keys(mc), Object.values(mc));

    // ── Mean run duration by medium (all completed runs) ─────────
    const durByMedium = {};
    data.filter(r => r.duration_s).forEach(r => {
      const med = classifyMedium(r);
      if (!durByMedium[med]) durByMedium[med] = [];
      durByMedium[med].push(r.duration_s / 60);
    });
    const durOrder  = ['water', 'phantom', 'layered', 'layered + pCT'];
    const durLabels = [...durOrder.filter(m => durByMedium[m]),
                       ...Object.keys(durByMedium).filter(m => !durOrder.includes(m))];
    const durMeans  = durLabels.map(m => {
      const vals = durByMedium[m];
      return +(vals.reduce((a, b) => a + b, 0) / vals.length).toFixed(1);
    });
    new Chart(document.getElementById('tc-duration'), {
      type: 'bar',
      data: {
        labels: durLabels,
        datasets: [{
          data: durMeans,
          backgroundColor: durLabels.map((_, i) => PAL[i % PAL.length]),
          borderRadius: 4,
          borderSkipped: false,
        }]
      },
      options: {
        plugins: { legend: { display: false } },
        scales: {
          x: { ticks: { color: fgL(), font: { size: 11 } }, grid: { display: false } },
          y: { ticks: { color: fgL(), callback: v => v + ' min' }, grid: { color: gridC() }, beginAtZero: true }
        }
      }
    });
  }

  if (typeof Chart !== 'undefined' && typeof supabase !== 'undefined') {
    init();
  } else {
    window.addEventListener('load', init);
  }
})();
</script>
