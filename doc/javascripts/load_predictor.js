// Load calculator (doc_hpc.md, Resource Estimates): the expected time of each stage of a run,
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
    '</div><div class="lp-out" aria-live="polite"></div>';
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

  ctl.solver.addEventListener('change', () => { fillDevices(); compute(); });
  root.querySelectorAll('[data-k]').forEach((el) => { if (el !== ctl.solver) el.addEventListener('input', compute); });
  root.querySelectorAll('select[data-k]').forEach((el) => { if (el !== ctl.solver) el.addEventListener('change', compute); });
  fillDevices();
  compute();
})();
