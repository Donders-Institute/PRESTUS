// Load calculator (doc_run-time.md): the expected time of each stage of a run,
// as the Nexus predicts it (nexus/prestus_nexus/workload.py). Data: load_model_data.js, generated
// from throughput.yaml by tests/docs/gen_load_data.py. Loaded on every page via extra_javascript;
// does nothing on pages without the #load-predictor container.
(function () {
  const root = document.getElementById('load-predictor');
  const M = window.PRESTUS_LOAD_MODEL;
  if (!root || !M) return;

  // ── formatting (workload._hms / _span) ──────────────────────────────────
  function hms(s) {
    const m = Math.round(s / 60);
    if (m >= 60) return Math.floor(m / 60) + 'h' + String(m % 60).padStart(2, '0');
    return s >= 90 ? m + ' min' : Math.max(Math.round(s), 1) + ' s';
  }
  function span(t) {
    const lo = hms(t[1]), hi = hms(t[2]);
    return lo === hi ? hms(t[0]) : hms(t[0]) + ' (' + lo + '–' + hi + ')';
  }
  function mega(n) { return n >= 1e6 ? (n / 1e6).toFixed(1) + ' M' : Math.round(n / 1e3) + ' k'; }
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  // ── k-Wave grid (workload.py) ───────────────────────────────────────────
  function maxPrime(n) {
    let p = n, f = 2, out = 1;
    while (f * f <= p) { while (p % f === 0) { out = f; p = p / f; } f++; }
    return p > 1 ? Math.max(out, p) : out;
  }
  function optimalPml(n) {           // getOptimalPMLSize: smallest largest prime factor in [10, 40]
    let best = 10, bestP = Infinity;
    for (let p = 10; p <= 40; p++) { const mp = maxPrime(n + 2 * p); if (mp < bestP) { bestP = mp; best = p; } }
    return best;
  }
  const withPml = (dims) => dims.reduce((a, n) => a * (n + 2 * optimalPml(n)), 1);
  const prod = (dims) => dims.reduce((a, n) => a * n, 1);

  function headDims(dx) {            // reference head (acoustic_grid)
    const scale = Math.cbrt(M.ref_head_l[1] * 1e6 / prod(M.ref_head_mm));
    return M.ref_head_mm.map((e) => Math.ceil(e * scale / dx));
  }
  function waterDims(dx, ap, fb) {   // water_grid_dims.m
    const lat = 2 * Math.ceil((ap + 2 * M.water_margin_mm) / dx / 2) + 1;
    return [lat, lat, Math.ceil((fb + M.z_beyond_focus_mm) / dx) + 2];
  }
  function fovDims(dx, d, len) {      // resolve_acoustic_fov_params.m: a cylinder around the beam
    const n = Math.ceil(d / dx);
    return [n, n, Math.ceil(len / dx)];
  }
  function baselineDims(dx, ap, fb, dims) {   // water_baseline_box.m
    const m = 2 * Math.ceil(10 / dx);
    const lat = Math.ceil((ap / dx + m) / 8) * 8, ax = Math.ceil(((fb + 40) / dx + m) / 8) * 8;
    const box = [lat, lat, ax];
    return prod(box) > 0.8 * prod(dims) ? dims : box;
  }
  function steps(dims, dx, c) {      // source_sensor_setup.m: dt = period / ceil(ppw / cfl), ppw = c_max / (f dx).
    // The frequency cancels (dt ~ cfl dx / c_max) up to the rounding of the points per period,
    // a few percent; it enters through the resolution it needs.
    const dt = M.cfl * (dx / 1e3) / c[1];
    const diag = Math.sqrt(dims.reduce((a, n) => a + (n * dx / 1e3) ** 2, 0));
    return Math.round(M.transits * diag / c[0] / dt);
  }

  // ── controls ───────────────────────────────────────────────────────────
  const SOLVERS = [['cpp_gpu', 'k-Wave C++, GPU'], ['matlab_gpu', 'k-Wave MATLAB, GPU'],
    ['cpp_cpu', 'k-Wave C++, CPU'], ['matlab_cpu', 'k-Wave MATLAB, CPU']];
  root.innerHTML =
    '<div class="lp-controls">' +
    '<label>Solver<select data-k="solver">' + SOLVERS.map((s) => '<option value="' + s[0] + '">' + s[1] +
      '</option>').join('') + '</select></label>' +
    '<label>Device<select data-k="dev"></select></label>' +
    '<label>Medium<select data-k="medium"><option value="head">head</option><option value="water">water</option></select></label>' +
    '<label>Resolution [mm]<select data-k="dx">' + [0.5, 0.6, 0.75, 1, 1.25, 1.5, 2].map((v) =>
      '<option' + (v === 0.5 ? ' selected' : '') + '>' + v + '</option>').join('') + '</select></label>' +
    '<label data-only="geom">Aperture [mm]<input data-k="ap" type="number" min="5" max="200" step="1" value="64"></label>' +
    '<label data-only="geom">Bowl to focus [mm]<input data-k="fb" type="number" min="5" max="200" step="1" value="85"></label>' +
    '<label data-only="head">Acoustic grid<select data-k="fov"><option value="full">full head</option>' +
    '<option value="fov">limited FOV</option></select></label>' +
    '<label data-only="fov">FOV diameter [mm]<input data-k="fovd" type="number" min="5" max="300" step="1" placeholder="auto: aperture"></label>' +
    '<label data-only="fov">FOV length [mm]<input data-k="fovl" type="number" min="5" max="300" step="1" placeholder="auto: focus + 40"></label>' +
    '<label data-only="head">Head model<select data-k="pre"><option value="preproc">built per run</option>' +
    '<option value="preproc_shared">shared across placements</option></select></label>' +
    '<label data-only="head" class="lp-check"><input data-k="base" type="checkbox" checked>Free-water baseline</label>' +
    '<label class="lp-check"><input data-k="heat" type="checkbox">Heating</label>' +
    '<label data-only="heat">Thermal steps<input data-k="nth" type="number" min="1" step="1" value="200"></label>' +
    '</div><div class="lp-out" aria-live="polite"></div>' +
    '<div class="lp-cmp"><div class="lp-cmp-head"><button type="button" class="lp-btn" data-act="pin">Pin setup to compare</button>' +
    '<label class="lp-check" hidden><input type="checkbox" data-act="log">Log scale (totals only)</label>' +
    '<button type="button" class="lp-btn lp-btn-quiet" data-act="clear" hidden>Clear all</button></div>' +
    '<div class="lp-cmp-out" aria-live="polite"></div></div>';
  const ctl = {};
  root.querySelectorAll('[data-k]').forEach((el) => { ctl[el.dataset.k] = el; });
  const out = root.querySelector('.lp-out');

  function fillDevices() {             // the devices with a rate for the solver (GPUs, or the CPU)
    const keep = ctl.dev.value;
    const key = ctl.solver.value + '|single';
    const opts = M.devices.map((d, i) => [d, i]).filter((x) => x[0].rates[key]);
    ctl.dev.innerHTML = opts.map((x) => '<option value="' + x[1] + '"' + (String(x[1]) === keep ? ' selected' : '') + '>' +
      esc(x[0].label) + '</option>').join('');
  }

  // ── prediction ─────────────────────────────────────────────────────────
  const STAGE_NAMES = {
    placement: 'Placement', preproc: 'Head model', preproc_shared: 'Head model (shared)', medium: 'Medium',
    nifti_medium: 'Medium maps', source: 'Source', acoustic_analysis: 'Acoustic analysis',
    nifti_acoustic: 'NIfTI maps', thermal_analysis: 'Thermal analysis', nifti_thermal: 'Thermal maps',
  };
  function stage(medium, label, dx) {
    const e = (M.stages[medium] || {})[label];
    if (!e) return null;
    const s = e.s_1mm * Math.pow(1 / dx, 3 * e.exponent);
    return [s, s * e.ratio[0], s * e.ratio[1]];
  }

  function compute() {
    const d = M.devices[+ctl.dev.value];
    const medium = ctl.medium.value, dx = +ctl.dx.value;
    const ap = +ctl.ap.value, fb = +ctl.fb.value, heat = ctl.heat.checked;
    root.querySelectorAll('[data-only="head"]').forEach((el) => { el.hidden = medium !== 'head'; });
    root.querySelectorAll('[data-only="heat"]').forEach((el) => { el.hidden = !heat; });
    const fov = medium === 'head' && ctl.fov.value === 'fov';
    root.querySelectorAll('[data-only="fov"]').forEach((el) => { el.hidden = !fov; });
    root.querySelectorAll('[data-only="geom"]').forEach((el) => { el.hidden = !(fov || medium === 'water'); });
    const r = d.rates[ctl.solver.value + '|single'];
    const dims = fov ? fovDims(dx, +ctl.fovd.value || ap, +ctl.fovl.value || fb + M.z_beyond_focus_mm)
      : medium === 'head' ? headDims(dx) : waterDims(dx, ap, fb);
    const c = M.speeds[medium], nt = steps(dims, dx, c), pts = withPml(dims);
    const rows = [];                  // [name, [s, lo, hi], kind]
    const sim = (p) => r ? [0, 1, 2].map((i) => r.t0 + [r.ns, r.lo, r.hi][i] * 1e-9 * p * nt) : null;
    rows.push(['Acoustic simulation', sim(pts), 'gpu']);
    if (medium === 'head' && ctl.base.checked) rows.push(['Free-water baseline', sim(withPml(baselineDims(dx, ap, fb, dims))), 'gpu']);
    if (heat) {
      const th = d.thermal, n = +ctl.nth.value;
      rows.push(['Thermal simulation', th ? [th.ns, th.lo, th.hi].map((ns) => ns * 1e-9 * prod(dims) * n) : null, 'heat']);
    }
    const labels = ['placement', medium === 'head' ? ctl.pre.value : 'preproc', 'medium', 'nifti_medium', 'source',
      'acoustic_analysis', 'nifti_acoustic'].concat(heat ? ['thermal_analysis', 'nifti_thermal'] : []);
    labels.forEach((lb) => { const s = stage(medium, lb, dx); if (s) rows.push([STAGE_NAMES[lb], s, 'cpu']); });
    // once per subject or only when used: shown, not added to the total
    const extra = (label, name) => { const s = stage(medium, label, dx); rows.push([name, s, s ? 'opt' : 'none']); };
    if (medium === 'head') extra('segmentation', 'Segmentation (SimNIBS charm, once per subject)');
    if (medium === 'head') extra('pct', 'Pseudo-CT (if used)');
    extra('report', 'HTML report (if on)');
    const timed = rows.filter((x) => x[1] && x[2] !== 'opt');
    const total = [0, 1, 2].map((i) => timed.reduce((a, x) => a + x[1][i], 0));
    const mem = d.kind === 'gpu' ? 40 * pts * 4 / 1e9 : null;
    render(rows, total, { dims, pts, nt, mem, d, r, medium, fov });
    // the settings that define this setup, for the comparison (label, value)
    const set = [['Solver', SOLVERS.find((x) => x[0] === ctl.solver.value)[1]], ['Device', d.label], ['Medium', medium],
      ['Resolution', dx + ' mm']];
    if (medium === 'head') set.push(['Acoustic grid', fov ? 'limited FOV' : 'full head']);
    if (fov || medium === 'water') set.push(['Aperture', ap + ' mm'], ['Bowl to focus', fb + ' mm']);
    if (fov) set.push(['FOV', (+ctl.fovd.value || ap) + ' × ' + (+ctl.fovl.value || fb + M.z_beyond_focus_mm) + ' mm']);
    if (medium === 'head') set.push(['Head model', ctl.pre.value === 'preproc' ? 'built per run' : 'shared'],
      ['Free-water baseline', ctl.base.checked ? 'on' : 'off']);
    set.push(['Heating', heat ? +ctl.nth.value + ' thermal steps' : 'off']);
    current = { set, total, mem, incomplete: rows.some((x) => !x[1] && x[2] !== 'none'),
      stages: rows.filter((x) => x[1] && x[2] !== 'opt').map((x) => [x[0], x[1][0], x[3]]) };
  }

  // ── rendering ──────────────────────────────────────────────────────────
  const COLOURS = { none: ['transparent'], opt: ['transparent'], gpu: ['#0e7490', '#38a3b5'], heat: ['#c2410c'], cpu: ['#57606a', '#8c959f', '#6e7781', '#afb8c1', '#424a53', '#9aa4ae', '#76808a', '#c3cad1', '#5c6670'] };
  function render(rows, total, info) {
    const used = { none: 0, opt: 0, gpu: 0, heat: 0, cpu: 0 };
    rows.forEach((x) => { x.push(COLOURS[x[2]][used[x[2]]++ % COLOURS[x[2]].length]); });
    const timed = rows.filter((x) => x[1] && x[2] !== 'opt');
    let xpos = 0, bars = '';
    timed.forEach((x) => {
      const w = total[0] > 0 ? 1000 * x[1][0] / total[0] : 0;
      bars += '<rect x="' + xpos.toFixed(1) + '" y="0" width="' + Math.max(w, 0).toFixed(1) + '" height="28" fill="' + x[3] +
        '"><title>' + esc(x[0]) + ': ' + span(x[1]) + '</title></rect>';
      xpos += w;
    });
    const where = { gpu: info.d.kind === 'gpu' ? 'GPU' : 'CPU', heat: info.d.kind === 'gpu' ? 'GPU' : 'CPU', cpu: 'CPU, any device',
      none: '', opt: 'CPU, not in the total' };
    const table = rows.map((x) => '<tr><td><span class="lp-swatch" style="background:' + x[3] + '"></span>' + esc(x[0]) +
      '</td><td>' + (x[1] ? hms(x[1][0]) : x[2] === 'none' ? 'not timed yet' : 'no rate for this device') + '</td><td>' +
      (x[1] ? hms(x[1][1]) + '–' + hms(x[1][2]) : '') + '</td><td>' + where[x[2]] + '</td></tr>').join('');
    const incomplete = rows.some((x) => !x[1] && x[2] !== 'none');
    out.innerHTML =
      '<div class="lp-total"><span class="lp-big">' + (timed.length ? span(total) : '–') + '</span>' +
      (incomplete ? ' <span class="lp-note">for the timed stages</span>' : '') +
      (info.mem !== null ? '<span class="lp-mem">GPU memory ~' + (info.mem >= 10 ? info.mem.toFixed(0) : info.mem.toFixed(1)) + ' GB</span>' : '') +
      '</div>' +
      '<svg class="lp-bar" viewBox="0 0 1000 28" preserveAspectRatio="none" role="img" aria-label="Time per stage">' + bars + '</svg>' +
      '<table class="lp-table"><thead><tr><th>Stage</th><th>Time</th><th>Range (p10–p90)</th><th>Runs on</th></tr></thead><tbody>' +
      table + '</tbody></table>' +
      '<p class="lp-note">Grid ' + info.dims.join(' × ') + (info.fov ? ' (limited FOV)' : info.medium === 'head' ? ' (reference head)' : ' (sized to the transducer)') +
      ', ' + mega(info.pts) + ' points with PML, ' + info.nt.toLocaleString('en') + ' time steps. Acoustic rate: ' +
      (info.r ? esc(info.r.basis) : 'none') + '.</p>';
  }

  // ── comparison: pinned setups as vertical stacked bars ─────────────────
  let current = null, pinned = [];
  const KEY = 'prestus-load-compare';
  try { pinned = JSON.parse(localStorage.getItem(KEY)) || []; } catch (e) { pinned = []; }
  const save = () => { try { localStorage.setItem(KEY, JSON.stringify(pinned)); } catch (e) { /* storage off: keep in memory */ } };
  const cmp = root.querySelector('.lp-cmp-out');
  const logBox = root.querySelector('[data-act="log"]');
  function niceTicks(max) {          // linear axis: 4–6 ticks at 1, 2, or 5 × 10^k seconds, or whole minutes / hours
    const units = [1, 2, 5, 10, 15, 30, 60, 120, 300, 600, 900, 1800, 3600, 7200, 10800, 21600, 43200, 86400];
    const step = units.find((u) => max / u <= 5) || Math.ceil(max / 5 / 86400) * 86400;
    const t = [];
    for (let v = 0; v <= max + 1e-9; v += step) t.push(v);
    if (t[t.length - 1] < max) t.push(t[t.length - 1] + step);
    return t;
  }
  function renderCompare() {
    root.querySelector('[data-act="clear"]').hidden = pinned.length === 0;
    logBox.parentElement.hidden = pinned.length < 2;
    if (!pinned.length) { cmp.innerHTML = ''; return; }
    const log = logBox.checked;
    const W = 640, H = 260, L = 56, B = 22, T = 10, bw = Math.min(70, (W - L - 10) / pinned.length * 0.6);
    const slot = (W - L - 10) / pinned.length;
    const maxT = Math.max.apply(null, pinned.map((p) => p.total[2]));
    let ticks, y;
    if (log) {
      const minT = Math.min.apply(null, pinned.map((p) => p.total[1]));
      const LOG_TICKS = [1, 10, 60, 600, 3600, 36000, 86400, 864000];   // 1 s, 10 s, 1 min, 10 min, 1 h, 10 h, 1 d, 10 d
      const i0 = Math.max(LOG_TICKS.filter((v) => v <= minT).length - 1, 0);
      const i1 = Math.max(LOG_TICKS.findIndex((v) => v >= maxT), i0 + 1);
      ticks = LOG_TICKS.slice(i0, (i1 < 0 ? LOG_TICKS.length - 1 : i1) + 1);
      const a = Math.log10(ticks[0]), b = Math.log10(ticks[ticks.length - 1]);
      y = (v) => T + (H - T - B) * (1 - (Math.log10(Math.max(v, ticks[0])) - a) / (b - a));
    } else {
      ticks = niceTicks(maxT);
      y = (v) => T + (H - T - B) * (1 - v / ticks[ticks.length - 1]);
    }
    let svg = '';
    ticks.forEach((v) => {
      svg += '<line class="lp-grid" x1="' + L + '" x2="' + W + '" y1="' + y(v).toFixed(1) + '" y2="' + y(v).toFixed(1) + '"/>' +
        '<text class="lp-axis" x="' + (L - 6) + '" y="' + (y(v) + 3).toFixed(1) + '" text-anchor="end">' + (v === 0 ? '0' : log ? ({ 1: '1 s', 10: '10 s', 60: '1 min', 600: '10 min', 3600: '1 h', 36000: '10 h', 86400: '1 d', 864000: '10 d' })[v] : hms(v)) + '</text>';
    });
    pinned.forEach((p, i) => {
      const cx = L + slot * (i + 0.5), x0 = (cx - bw / 2).toFixed(1);
      if (log) {
        svg += '<rect x="' + x0 + '" y="' + y(p.total[0]).toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' +
          (H - B - y(p.total[0])).toFixed(1) + '" class="lp-tot"><title>#' + (i + 1) + ': ' + span(p.total) + '</title></rect>';
      } else {
        let acc = 0;
        p.stages.forEach((s) => {
          const y1 = y(acc + s[1]), y0 = y(acc);
          svg += '<rect x="' + x0 + '" y="' + y1.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + Math.max(y0 - y1, 0).toFixed(1) +
            '" fill="' + s[2] + '"><title>#' + (i + 1) + ' ' + esc(s[0]) + ': ' + hms(s[1]) + '</title></rect>';
          acc += s[1];
        });
      }
      // whisker: p10–p90 of the total
      const wl = y(p.total[1]), wh = y(p.total[2]);
      svg += '<path class="lp-whisk" d="M' + cx.toFixed(1) + ' ' + wl.toFixed(1) + 'V' + wh.toFixed(1) +
        'M' + (cx - 6).toFixed(1) + ' ' + wh.toFixed(1) + 'h12M' + (cx - 6).toFixed(1) + ' ' + wl.toFixed(1) + 'h12"/>' +
        '<text class="lp-axis lp-num" x="' + cx.toFixed(1) + '" y="' + (H - 6) + '" text-anchor="middle">#' + (i + 1) + '</text>' +
        '<text class="lp-axis" x="' + cx.toFixed(1) + '" y="' + (Math.min(wh, y(p.total[0])) - 4).toFixed(1) + '" text-anchor="middle">' + hms(p.total[0]) + '</text>';
    });
    // legend: one swatch per stage name, coloured as in the setups
    const seen = {};
    pinned.forEach((p) => p.stages.forEach((s) => { if (!(s[0] in seen)) seen[s[0]] = s[2]; }));
    const legend = log ? '' : '<div class="lp-legend">' + Object.keys(seen).map((n) =>
      '<span><span class="lp-swatch" style="background:' + seen[n] + '"></span>' + esc(n) + '</span>').join('') + '</div>';
    // settings next to each number; from #2 on, what differs from #1 in bold
    const ref = {}; pinned[0].set.forEach((s) => { ref[s[0]] = s[1]; });
    const list = pinned.map((p, i) => '<li><span class="lp-num">#' + (i + 1) + '</span><span class="lp-set">' +
      p.set.map((s) => { const diff = i > 0 && ref[s[0]] !== s[1];
        return (diff ? '<b>' : '') + esc(s[0]) + ': ' + esc(s[1]) + (diff ? '</b>' : ''); }).join(' · ') +
      '</span><span class="lp-settot">' + span(p.total) + (p.incomplete ? ' (timed stages)' : '') +
      (p.mem !== null ? ' · GPU ~' + (p.mem >= 10 ? p.mem.toFixed(0) : p.mem.toFixed(1)) + ' GB' : '') + '</span>' +
      '<span class="lp-rowbtns"><button type="button" class="lp-btn lp-btn-quiet" data-mv="' + i + ':-1" aria-label="Move setup ' + (i + 1) +
      ' up"' + (i === 0 ? ' disabled' : '') + '>↑</button><button type="button" class="lp-btn lp-btn-quiet" data-mv="' + i + ':1" aria-label="Move setup ' +
      (i + 1) + ' down"' + (i === pinned.length - 1 ? ' disabled' : '') + '>↓</button>' +
      '<button type="button" class="lp-btn lp-btn-quiet" data-rm="' + i + '" aria-label="Remove setup ' + (i + 1) + '">Remove</button></span></li>').join('');
    cmp.innerHTML = '<svg class="lp-cmp-chart" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Predicted time per pinned setup">' +
      svg + '</svg>' + legend + '<ol class="lp-setups">' + list + '</ol>' +
      '<p class="lp-note">Bars stack the stages in the total; whiskers span the 10th–90th percentile of the total. ' +
      (log ? 'On the log scale, bars show the total only. ' : '') + '↑ and ↓ reorder; setup 1 is the reference for the settings in bold. Pinned setups are kept in this browser only.</p>';
    cmp.querySelectorAll('[data-rm]').forEach((b) => b.addEventListener('click', () => {
      pinned.splice(+b.dataset.rm, 1); save(); renderCompare(); }));
    cmp.querySelectorAll('[data-mv]').forEach((b) => b.addEventListener('click', () => {
      const [i, d] = b.dataset.mv.split(':').map(Number), j = i + d;
      if (j < 0 || j >= pinned.length) return;
      [pinned[i], pinned[j]] = [pinned[j], pinned[i]]; save(); renderCompare();
      const f = cmp.querySelector('[data-mv="' + j + ':' + d + '"]') || cmp.querySelector('[data-mv="' + j + ':' + (-d) + '"]');
      if (f) f.focus(); }));
  }
  root.querySelector('[data-act="pin"]').addEventListener('click', () => {
    if (current && current.stages.length) { pinned.push(JSON.parse(JSON.stringify(current))); save(); renderCompare(); }
  });
  root.querySelector('[data-act="clear"]').addEventListener('click', () => { pinned = []; save(); renderCompare(); });
  logBox.addEventListener('change', renderCompare);

  ctl.solver.addEventListener('change', () => { fillDevices(); compute(); });
  root.querySelectorAll('[data-k]').forEach((el) => { if (el !== ctl.solver) el.addEventListener('input', compute); });
  root.querySelectorAll('select[data-k]').forEach((el) => { if (el !== ctl.solver) el.addEventListener('change', compute); });
  fillDevices();
  compute();
  renderCompare();
})();
