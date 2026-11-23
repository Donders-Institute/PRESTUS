// Interactive illustration for doc_multitransducer.md ("How the Fields Combine"):
// two focused transducers on opposite sides of a 2-D domain, combined as
// coherent at one frequency, coherent at different frequencies, or async
// (firing in turn). Schematic: lossless, homogeneous, 2-D, slowed down; the
// physics of the combination is the one PRESTUS uses.
// Renders into every <div class="mt-fields"></div> on the page.
(function () {
  "use strict";

  var W = 220, H = 140;                 // field resolution [px]
  var LAMBDA1 = 22;                     // wavelength of transducer 1 [px]
  var FOCUS = [W / 2, H / 2];
  var DISPLAY_HZ = 0.6;                 // cycles of transducer 1 per second on screen
  var GATE_CYCLES = 6;                  // async: cycles per turn (pulses shortened for display)
  var running = [];

  function arcElements(cx, side) {
    // focused arc (bowl section) of radius |cx - focus|, aperture +-35 deg, facing the focus
    var R = Math.abs(FOCUS[0] - cx), els = [];
    for (var a = -35; a <= 35; a += 5) {
      var t = a * Math.PI / 180;
      els.push([FOCUS[0] - side * R * Math.cos(t), FOCUS[1] + R * Math.sin(t)]);
    }
    return els;
  }
  var TX = [arcElements(18, 1), arcElements(W - 18, -1)];
  var TXC = ["#1f6fd1", "#c8384f"];      // transducer 1, 2 (arcs and trace)

  // Complex single-frequency field of one transducer (unit drive): sum over
  // elements of exp(-i k r) / sqrt(r), normalised to 1 at the focus.
  function field(els, lambda) {
    var k = 2 * Math.PI / lambda, re = new Float32Array(W * H), im = new Float32Array(W * H);
    for (var y = 0; y < H; y++) {
      for (var x = 0; x < W; x++) {
        var sr = 0, si = 0;
        for (var e = 0; e < els.length; e++) {
          var dx = x - els[e][0], dy = y - els[e][1];
          var r = Math.max(Math.sqrt(dx * dx + dy * dy), 1.5);
          var a = 1 / Math.sqrt(r);
          sr += a * Math.cos(k * r); si -= a * Math.sin(k * r);
        }
        re[y * W + x] = sr; im[y * W + x] = si;
      }
    }
    var i0 = Math.round(FOCUS[1]) * W + Math.round(FOCUS[0]);
    var nr = re[i0], ni = im[i0], n2 = nr * nr + ni * ni;   // divide by the focal value
    for (var j = 0; j < W * H; j++) {
      var r0 = re[j], i1 = im[j];
      re[j] = (r0 * nr + i1 * ni) / n2; im[j] = (i1 * nr - r0 * ni) / n2;
    }
    return { re: re, im: im };
  }

  // Diverging blue-white-red for pressure, dark-to-yellow for intensity.
  function diverging(v) {
    v = Math.max(-1, Math.min(1, v));
    if (v >= 0) return [255, Math.round(255 * (1 - v)), Math.round(255 * (1 - v))];
    return [Math.round(255 * (1 + v)), Math.round(255 * (1 + v)), 255];
  }
  function sequential(v) {
    v = Math.max(0, Math.min(1, v));
    return [Math.round(255 * Math.min(1, 1.6 * v)), Math.round(255 * Math.max(0, 1.4 * v - 0.3)),
            Math.round(255 * Math.max(0, 0.45 - 0.9 * Math.abs(v - 0.35)))];
  }

  function cssVar(name, fallback) {
    var v = getComputedStyle(document.body).getPropertyValue(name).trim();
    return v || fallback;
  }

  var MODES = {
    single: {
      label: "Single transducer",
      text: "Transducer 1 alone: the reference for the readouts (1 = one transducer at its focus). The fringes of a single focused beam stand still and its time-averaged intensity peaks at the focus. PRESTUS simulates it in one run."
    },
    same: {
      label: "Coherent, one frequency",
      text: "Both transducers fire together at one frequency. Their relative phase is fixed, so the fields add as complex numbers: the interference fringes stand still, and the time-averaged intensity at a point can reach four times that of one transducer (or cancel). PRESTUS simulates them in one run."
    },
    diff: {
      label: "Coherent, different frequencies",
      text: "Both fire together at different frequencies. The fringes drift: the relative phase changes at the beat frequency |f1 - f2|. Peaks still add, up to the sum of the amplitudes when the tones line up; for low-ratio pairs such as 1:2 the peak depends on their starting phases (try the slider). The time-averaged intensity is the sum of the two intensities whatever the phase. PRESTUS simulates each transducer alone and combines them this way."
    },
    async: {
      label: "Async, firing in turn",
      text: "The transducers fire in separate windows (pulses shortened here; blue: transducer 1, red: transducer 2). Only one field exists at a time: the peak pressure is the larger of the two, and heating follows each transducer's intensity while it is on, here half the time each. PRESTUS simulates each transducer alone."
    }
  };

  function build(root) {
    root.innerHTML = "";
    var ui = document.createElement("div");
    ui.className = "mt-fields__controls";
    var modeSel = document.createElement("div");
    modeSel.className = "mt-fields__modes";
    modeSel.setAttribute("role", "radiogroup");
    var state = { mode: "same", ratio: 1.5, phase: 0, t: 2 * GATE_CYCLES + 1, playing: true, last: null, U: null };
    try {   // ?mt=diff (or same, async) opens a setup directly
      var q = new URLSearchParams(window.location.search).get("mt");
      if (q && MODES[q]) state.mode = q;
    } catch (e) { /* old browsers */ }

    Object.keys(MODES).forEach(function (m) {
      var b = document.createElement("button");
      b.type = "button"; b.textContent = MODES[m].label; b.dataset.mode = m;
      b.setAttribute("role", "radio");
      b.addEventListener("click", function () { state.mode = m; refresh(); });
      modeSel.appendChild(b);
    });
    ui.appendChild(modeSel);

    var opts = document.createElement("div");
    opts.className = "mt-fields__options";
    opts.innerHTML =
      '<label class="mt-fields__phase">Phase of transducer 2 at t = 0 <input type="range" min="0" max="360" step="5" value="0"> <span>0°</span></label>' +
      '<label class="mt-fields__ratio">f2 / f1 <select><option value="1.1">1.1</option><option value="1.5" selected>1.5 (2:3)</option><option value="2">2 (1:2)</option></select></label>' +
      '<button type="button" class="mt-fields__play">Pause</button>';
    ui.appendChild(opts);
    root.appendChild(ui);

    var panels = document.createElement("div");
    panels.className = "mt-fields__panels";
    function panel(title) {
      var f = document.createElement("figure");
      var c = document.createElement("canvas");
      c.width = W; c.height = H;
      var cap = document.createElement("figcaption");
      cap.textContent = title;
      f.appendChild(c); f.appendChild(cap); panels.appendChild(f);
      return c;
    }
    var cP = panel("Pressure now");
    var cI = panel("Time-averaged intensity (1 = one transducer at its focus)");
    root.appendChild(panels);

    var trace = document.createElement("canvas");
    trace.className = "mt-fields__trace";
    trace.width = 660; trace.height = 140;
    root.appendChild(trace);
    var proto = document.createElement("canvas");
    proto.className = "mt-fields__trace";
    proto.width = 660; proto.height = 96;
    root.appendChild(proto);
    var readout = document.createElement("p");
    readout.className = "mt-fields__readout";
    root.appendChild(readout);
    var text = document.createElement("p");
    text.className = "mt-fields__text";
    root.appendChild(text);

    var phaseIn = opts.querySelector(".mt-fields__phase input");
    var phaseOut = opts.querySelector(".mt-fields__phase span");
    var ratioIn = opts.querySelector(".mt-fields__ratio select");
    var playBtn = opts.querySelector(".mt-fields__play");
    phaseIn.addEventListener("input", function () { state.phase = +phaseIn.value; phaseOut.textContent = phaseIn.value + "°"; drawStatic(); });
    ratioIn.addEventListener("change", function () { state.ratio = +ratioIn.value; refresh(); });
    playBtn.addEventListener("click", function () { state.playing = !state.playing; playBtn.textContent = state.playing ? "Pause" : "Play"; });
    if (window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      state.playing = false; playBtn.textContent = "Play";
    }

    var ctxP = cP.getContext("2d"), ctxI = cI.getContext("2d"), ctxT = trace.getContext("2d");
    var ctxS = proto.getContext("2d");
    var imgP = ctxP.createImageData(W, H), imgI = ctxI.createImageData(W, H);

    function f2() { return state.mode === "diff" ? state.ratio : 1; }

    function refresh() {
      modeSel.querySelectorAll("button").forEach(function (b) {
        b.setAttribute("aria-checked", b.dataset.mode === state.mode ? "true" : "false");
      });
      opts.querySelector(".mt-fields__phase").style.display = (state.mode === "same" || state.mode === "diff") ? "" : "none";
      opts.querySelector(".mt-fields__ratio").style.display = state.mode === "diff" ? "" : "none";
      state.U = [field(TX[0], LAMBDA1), field(TX[1], LAMBDA1 / f2())];
      text.textContent = MODES[state.mode].text;
      drawStatic();
      drawProtocol();
    }

    // complex amplitude of transducer 2 including its phase (one-frequency mode)
    function u2(j) {
      var c = Math.cos(state.phase * Math.PI / 180), s = Math.sin(state.phase * Math.PI / 180);
      var r = state.U[1].re[j], i = state.U[1].im[j];
      return state.mode !== "async" ? [r * c - i * s, r * s + i * c] : [r, i];
    }

    function pressureAt(j, t) {
      // t in periods of transducer 1; returns p(t) at pixel j
      var a1 = state.U[0], w1 = 2 * Math.PI * t, w2 = w1 * f2(), b = u2(j);
      var p1 = a1.re[j] * Math.cos(w1) - a1.im[j] * Math.sin(w1);
      var p2 = b[0] * Math.cos(w2) - b[1] * Math.sin(w2);
      if (state.mode === "single") return p1;
      if (state.mode === "async") {
        return activeTx(t) === 0 ? p1 : p2;
      }
      return p1 + p2;
    }

    function intensityAt(j) {
      var a = state.U[0], b = u2(j);
      var i1 = a.re[j] * a.re[j] + a.im[j] * a.im[j], i2 = b[0] * b[0] + b[1] * b[1];
      if (state.mode === "single") return i1;
      if (state.mode === "same") {
        var r = a.re[j] + b[0], i = a.im[j] + b[1];
        return r * r + i * i;
      }
      if (state.mode === "diff") return i1 + i2;
      return (i1 + i2) / 2;   // each on half the time
    }

    function drawStatic() {
      var d = imgI.data;
      for (var j = 0; j < W * H; j++) {
        var c = sequential(intensityAt(j) / 4);    // colour scale up to 4 (two in phase)
        d[4 * j] = c[0]; d[4 * j + 1] = c[1]; d[4 * j + 2] = c[2]; d[4 * j + 3] = 255;
      }
      ctxI.putImageData(imgI, 0, 0);
      marks(ctxI);
      // readouts at the focus
      var j0 = Math.round(FOCUS[1]) * W + Math.round(FOCUS[0]);
      var pk = 0, n = 2000, span = state.mode === "diff" ? 20 : 2 * GATE_CYCLES;
      for (var k = 0; k < n; k++) pk = Math.max(pk, Math.abs(pressureAt(j0, k * span / n)));
      readout.textContent = "At the focus, relative to one transducer: peak pressure " + pk.toFixed(2) +
        " x, time-averaged intensity " + intensityAt(j0).toFixed(2) + " x.";
    }

    function activeTx(t) {   // async: which transducer is on at time t (periods of transducer 1)
      return ((Math.floor(t / GATE_CYCLES) % 2) + 2) % 2;
    }
    // On-times over a few pulse repetition intervals, one row per transducer
    // (as the firing-schedule plot of PRESTUS): together in both coherent setups,
    // offset by PRI/2 in async mode (timing.async_schedule: pulses).
    function drawProtocol() {
      var tw = proto.width, th = proto.height, n_pri = 4, pd = 0.3;
      var x0 = 92, x1 = tw - 8, pri_px = (x1 - x0) / n_pri;
      var fg = cssVar("--md-default-fg-color", "#333"), light = cssVar("--md-default-fg-color--lightest", "#ccc");
      ctxS.clearRect(0, 0, tw, th);
      ctxS.font = "12px system-ui, sans-serif";
      ctxS.fillStyle = fg;
      ctxS.fillText("On-times (schematic; the animation above shows the time within a pulse)", 6, 13);
      [0, 1].forEach(function (i) {
        var yc = 34 + i * 26;
        ctxS.fillStyle = fg;
        ctxS.fillText("Transducer " + (i + 1), 6, yc + 4);
        ctxS.strokeStyle = light; ctxS.beginPath(); ctxS.moveTo(x0, yc + 8); ctxS.lineTo(x1, yc + 8); ctxS.stroke();
        var offset = (state.mode === "async" && i === 1) ? 0.5 : 0;
        ctxS.fillStyle = TXC[i];
        if (state.mode === "single" && i === 1) return;   // transducer 2 not used
        for (var k = 0; k < n_pri; k++) {
          var xa = x0 + (k + offset) * pri_px;
          ctxS.fillRect(xa, yc - 8, pd * pri_px, 16);
        }
      });
      // PD and PRI brackets under the rows
      var yb = th - 10;
      ctxS.strokeStyle = fg; ctxS.fillStyle = fg; ctxS.beginPath();
      ctxS.moveTo(x0, yb - 4); ctxS.lineTo(x0, yb); ctxS.lineTo(x0 + pd * pri_px, yb); ctxS.lineTo(x0 + pd * pri_px, yb - 4);
      ctxS.moveTo(x0 + pri_px, yb - 4); ctxS.lineTo(x0 + pri_px, yb); ctxS.lineTo(x0 + 2 * pri_px, yb); ctxS.lineTo(x0 + 2 * pri_px, yb - 4);
      ctxS.stroke();
      ctxS.fillText("PD", x0 + pd * pri_px / 2 - 8, yb + 1 - 14 > 0 ? yb - 6 : yb);
      ctxS.fillText("PRI", x0 + 1.5 * pri_px - 10, yb - 6);
    }

    function marks(ctx, t) {
      TX.forEach(function (els, i) {
        var off = (state.mode === "single" && i === 1) ||
                  (state.mode === "async" && t !== undefined && activeTx(t) !== i);
        ctx.fillStyle = off ? "rgba(120,120,120,0.45)" : TXC[i];
        els.forEach(function (e) { ctx.fillRect(e[0] - 1.5, e[1] - 1.5, 3, 3); });
      });
      ctx.strokeStyle = "#222"; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(FOCUS[0], FOCUS[1], 3, 0, 2 * Math.PI); ctx.stroke();
    }

    function drawFrame() {
      var d = imgP.data;
      for (var j = 0; j < W * H; j++) {
        var c = diverging(pressureAt(j, state.t) / 2);
        d[4 * j] = c[0]; d[4 * j + 1] = c[1]; d[4 * j + 2] = c[2]; d[4 * j + 3] = 255;
      }
      ctxP.putImageData(imgP, 0, 0);
      marks(ctxP, state.t);
      drawProtocol();   // cheap; follows a theme switch
      // trace at the focus over a window ending now
      var j0 = Math.round(FOCUS[1]) * W + Math.round(FOCUS[0]);
      var span = state.mode === "diff" ? 20 : 2 * GATE_CYCLES, tw = trace.width, th = trace.height;
      var fg = cssVar("--md-default-fg-color", "#333"), light = cssVar("--md-default-fg-color--lightest", "#ccc");
      var top = 20, mid = top + (th - top) / 2, sc = (th - top) / 2 / 2.2;   // +-2.2 fits two in phase
      ctxT.clearRect(0, 0, tw, th);
      ctxT.strokeStyle = light; ctxT.beginPath();
      [-2, -1, 0, 1, 2].forEach(function (v) { var y = mid - v * sc; ctxT.moveTo(0, y); ctxT.lineTo(tw, y); });
      ctxT.stroke();
      ctxT.lineWidth = 1.5;
      var prev = null, col = null;
      for (var x = 0; x < tw; x++) {
        var t = state.t - span + span * x / tw;
        var y = mid - pressureAt(j0, t) * sc;
        var c = state.mode === "async" ? TXC[activeTx(t)] : fg;
        if (c !== col) {   // async: colour by the transducer that is on
          if (col !== null) ctxT.stroke();
          ctxT.strokeStyle = c; ctxT.beginPath();
          if (prev) ctxT.moveTo(prev[0], prev[1]); else ctxT.moveTo(x, y);
          col = c;
        }
        ctxT.lineTo(x, y);
        prev = [x, y];
      }
      ctxT.stroke();
      ctxT.fillStyle = fg; ctxT.font = "12px system-ui, sans-serif";
      ctxT.fillText("Pressure at the focus over time (grid lines at ±1 and ±2; 1 = one transducer)", 6, 13);
    }

    var entry = { root: root, raf: null };
    function loop(ts) {
      if (!document.body.contains(root)) { return; }   // page left (instant navigation)
      if (state.last !== null && state.playing) state.t += (ts - state.last) / 1000 * DISPLAY_HZ;
      state.last = ts;
      drawFrame();
      entry.raf = requestAnimationFrame(loop);
    }
    refresh();
    entry.raf = requestAnimationFrame(loop);
    running.push(entry);
  }

  function init() {
    running.forEach(function (e) { if (e.raf) cancelAnimationFrame(e.raf); });
    running = [];
    document.querySelectorAll(".mt-fields").forEach(build);
  }

  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(init);          // MkDocs Material, also after instant navigation
  } else if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
