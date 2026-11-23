// Interactive illustration for doc_time_reversal.md: the three steps of time
// reversal as PRESTUS computes them (time_reversal_refocus.m).
//   1. planned drive: geometric phases, aimed at the target through the skull;
//   2. reverse run: a point source at the target; by reciprocity, each element's
//      transfer to the target is G_e = sum over its points of the reverse field;
//   3. refocused drive: phases -angle(G_e), amplitudes unchanged; per element
//      (matrix), per ring (annular; here a ring is a mirrored pair of the
//      cross-section), or none (single element).
// Readouts use the code's metrics: coherence_elements = |sum C_e G_e| / sum |C_e G_e|,
// and the target pressure relative to the planned drive.
// Schematic: 2-D, straight rays with a skull phase delay and loss, slowed down.
// Renders into every <div class="tr-fields"></div> on the page.
(function () {
  "use strict";

  var W = 220, H = 160;                    // field resolution [px]
  var LAMBDA = 14;                         // wavelength in water [px]
  var K0 = 2 * Math.PI / LAMBDA;
  var C_RATIO = 1.9;                       // skull / water sound speed
  var SKULL_LOSS = 0.025;                  // amplitude loss per px of skull
  var TARGET = [110, 118];
  var R_BOWL = 100, APERTURE = 40;         // bowl radius [px] (centred on the target), half-aperture [deg]
  var N_EL = 12, N_SUB = 4;                // elements in the cross-section, points per element
  var HEAD_C = [110, 330], R_IN = 260;     // skull: inner radius about HEAD_C (between bowl and target)
  var DISPLAY_HZ = 0.6;
  var running = [];

  function thickness(x) {                  // skull thickness [px]: thicker to the right, one bump
    return 8 + 22 * Math.max(0, x - 40) / 140 + 10 * Math.exp(-Math.pow((x - 85) / 16, 2));
  }
  function inSkull(x, y) {
    var dx = x - HEAD_C[0], dy = y - HEAD_C[1], r = Math.sqrt(dx * dx + dy * dy);
    return r >= R_IN && r <= R_IN + thickness(x);
  }
  var SK = new Uint8Array(W * H);          // skull mask on the display grid
  for (var yy = 0; yy < H; yy++) for (var xx = 0; xx < W; xx++) SK[yy * W + xx] = inSkull(xx, yy) ? 1 : 0;
  function skullLength(x0, y0, x1, y1) {   // length of the straight path inside the skull
    var L = Math.hypot(x1 - x0, y1 - y0), n = Math.max(8, Math.ceil(L / 2)), hit = 0;
    for (var s = 0; s < n; s++) {
      var f = (s + 0.5) / n, x = Math.round(x0 + f * (x1 - x0)), y = Math.round(y0 + f * (y1 - y0));
      if (x >= 0 && x < W && y >= 0 && y < H && SK[y * W + x]) hit++;
    }
    return L * hit / n;
  }
  // complex propagator from (x0,y0) to (x1,y1): exp(-i phase) * amplitude
  function prop(x0, y0, x1, y1) {
    var r = Math.max(Math.hypot(x1 - x0, y1 - y0), 1.5), ls = skullLength(x0, y0, x1, y1);
    var ph = K0 * (r - ls) + (K0 / C_RATIO) * ls;          // slower phase accrual in faster skull
    var a = Math.exp(-SKULL_LOSS * ls) / Math.sqrt(r);
    return [a * Math.cos(ph), -a * Math.sin(ph)];
  }

  // element sample points on the bowl (arc about the target, facing down)
  var PTS = [], ELEM = [];
  (function () {
    var step = 2 * APERTURE / (N_EL * N_SUB);
    for (var e = 0; e < N_EL; e++) {
      ELEM.push([]);
      for (var s = 0; s < N_SUB; s++) {
        var a = (-APERTURE + step * (e * N_SUB + s + 0.5)) * Math.PI / 180;
        PTS.push([TARGET[0] + R_BOWL * Math.sin(a), TARGET[1] - R_BOWL * Math.cos(a)]);
        ELEM[e].push(PTS.length - 1);
      }
    }
  })();
  function elemCentre(e) {
    var x = 0, y = 0;
    ELEM[e].forEach(function (p) { x += PTS[p][0]; y += PTS[p][1]; });
    return [x / N_SUB, y / N_SUB];
  }

  // Propagators of every point (forward) and of the target (reverse), computed once.
  var FW = null, REV = null, GP = null;
  function precompute() {
    if (FW) return;
    FW = PTS.map(function (p) {
      var re = new Float32Array(W * H), im = new Float32Array(W * H);
      for (var y = 0; y < H; y++) for (var x = 0; x < W; x++) {
        var c = prop(p[0], p[1], x, y); re[y * W + x] = c[0]; im[y * W + x] = c[1];
      }
      return { re: re, im: im };
    });
    REV = { re: new Float32Array(W * H), im: new Float32Array(W * H) };
    for (var y = 0; y < H; y++) for (var x = 0; x < W; x++) {
      var c = prop(TARGET[0], TARGET[1], x, y); REV.re[y * W + x] = c[0]; REV.im[y * W + x] = c[1];
    }
    GP = PTS.map(function (q) { return prop(TARGET[0], TARGET[1], q[0], q[1]); });   // reverse field at each point
  }

  // The elements of each transducer type: lists of points with one drive each.
  //   matrix: every element; annular: one ring = the mirrored pair of the cross-section;
  //   single: one element over the whole aperture.
  function groups(type) {
    if (type === "matrix") return ELEM;
    if (type === "single") return [PTS.map(function (_, i) { return i; })];
    var out = [];
    for (var e = 0; e < N_EL / 2; e++) out.push(ELEM[e].concat(ELEM[N_EL - 1 - e]));
    return out;
  }
  // G_e = sum over the element's points of the reverse field there (reciprocity)
  function transfers(gr) {
    return gr.map(function (pts) {
      var r = 0, i = 0;
      pts.forEach(function (p) { r += GP[p][0]; i += GP[p][1]; });
      return [r, i];
    });
  }

  function cmul(a, b) { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }
  function cabs(a) { return Math.hypot(a[0], a[1]); }
  function unit(ph) { return [Math.cos(ph), Math.sin(ph)]; }

  // drives C_e (unit amplitude): geometric (planned) or time reversal, -angle(G_e)
  function drives(kind, gr) {
    if (kind !== "refocus") return gr.map(function () { return [1, 0]; });   // bowl centred on the target
    return transfers(gr).map(function (g) { return unit(-Math.atan2(g[1], g[0])); });
  }
  function targetPressure(gr, C) {
    var s = [0, 0];
    transfers(gr).forEach(function (g, e) { var t = cmul(C[e], g); s[0] += t[0]; s[1] += t[1]; });
    return s;
  }
  // coherence_elements: |sum C_e G_e| / sum |C_e G_e|
  function coherenceElements(gr, C) {
    var den = 0;
    transfers(gr).forEach(function (g, e) { den += cabs(cmul(C[e], g)); });
    return cabs(targetPressure(gr, C)) / den;
  }
  // coherence_aperture: the same over all points, including the spread within each element
  function coherenceAperture(gr, C) {
    var num = [0, 0], den = 0;
    gr.forEach(function (pts, e) {
      pts.forEach(function (p) { var t = cmul(C[e], GP[p]); num[0] += t[0]; num[1] += t[1]; den += cabs(t); });
    });
    return cabs(num) / den;
  }
  function forwardField(gr, C) {
    var re = new Float32Array(W * H), im = new Float32Array(W * H);
    gr.forEach(function (pts, e) {
      var c = C[e];
      pts.forEach(function (p) {
        var f = FW[p];
        for (var j = 0; j < W * H; j++) {
          re[j] += c[0] * f.re[j] - c[1] * f.im[j];
          im[j] += c[0] * f.im[j] + c[1] * f.re[j];
        }
      });
    });
    return { re: re, im: im };
  }

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
  function hue(ph) {                         // phase -> colour (cyclic)
    var h = ((ph / (2 * Math.PI)) % 1 + 1) % 1, f = function (n) {
      var k = (n + h * 6) % 6; return Math.round(255 * (1 - Math.max(0, Math.min(k, 4 - k, 1))));
    };
    return "rgb(" + f(5) + "," + f(3) + "," + f(1) + ")";
  }

  var STEPS = {
    planned: {
      label: "1  Planned drive",
      text: "The transducer fires with the geometric phases of the planned pose, aimed at the target (circle). The skull (grey band) is thicker on the right and has a bump: sound is faster in bone, so the parts of the wave front that cross more bone arrive early. The waves no longer meet in phase at the target: the focus weakens and shifts (+ marks the intensity peak)."
    },
    reverse: {
      label: "2  Reverse run",
      text: "One extra simulation: a point source at the target, on the same grid and medium. By reciprocity, the field it produces at a point on the transducer equals the field a source there produces at the target. Summed over each element's points, it gives the element's transfer G_e: with what amplitude and phase that element reaches the target through the skull. The element colours show the arrival phase; equal colours would mean no aberration."
    },
    refocus: {
      label: "3  Refocused drive",
      text: ""
    }
  };
  var TYPES = {
    matrix: { label: "Matrix array", text: "Each element gets the phase -angle(G_e) (phase conjugation), with its amplitude unchanged. All elements now arrive at the target in phase: the largest target pressure these amplitudes allow." },
    annular: { label: "Annular array", text: "One phase per ring (here a ring is the pair of elements mirrored about the axis). Rings can only correct the part of the aberration that is the same on both sides; the left-right difference remains, so the gain is smaller than with a matrix array." },
    single: { label: "Single element", text: "A single element is one transducer surface with one drive: there is no phase to set, so the refocused drive is the planned one. PRESTUS still reports the aperture coherence and the predicted focal shift, and the pose search can look for a position where the skull aberrates less." }
  };

  function build(root) {
    root.innerHTML = "";
    var state = { step: "planned", type: "matrix", t: 0, playing: true, last: null };
    try {   // ?tr=reverse (or planned, refocus) and &type=annular (or matrix, single) open a view directly
      var qs = new URLSearchParams(window.location.search);
      if (STEPS[qs.get("tr")]) state.step = qs.get("tr");
      if (TYPES[qs.get("type")]) state.type = qs.get("type");
    } catch (e) { /* old browsers */ }
    var ui = document.createElement("div");
    ui.className = "mt-fields__controls";
    var stepSel = document.createElement("div");
    stepSel.className = "mt-fields__modes"; stepSel.setAttribute("role", "radiogroup");
    stepSel.setAttribute("aria-label", "Step");
    Object.keys(STEPS).forEach(function (s) {
      var b = document.createElement("button");
      b.type = "button"; b.textContent = STEPS[s].label; b.dataset.step = s; b.setAttribute("role", "radio");
      b.addEventListener("click", function () { state.step = s; refresh(); });
      stepSel.appendChild(b);
    });
    ui.appendChild(stepSel);
    var opts = document.createElement("div");
    opts.className = "mt-fields__options";
    var typeSel = document.createElement("div");
    typeSel.className = "mt-fields__modes"; typeSel.setAttribute("role", "radiogroup");
    typeSel.setAttribute("aria-label", "Transducer type");
    Object.keys(TYPES).forEach(function (k) {
      var b = document.createElement("button");
      b.type = "button"; b.textContent = TYPES[k].label; b.dataset.type = k; b.setAttribute("role", "radio");
      b.addEventListener("click", function () { state.type = k; refresh(); });
      typeSel.appendChild(b);
    });
    opts.appendChild(typeSel);
    var playBtn = document.createElement("button");
    playBtn.type = "button"; playBtn.className = "mt-fields__play"; playBtn.textContent = "Pause";
    playBtn.addEventListener("click", function () { state.playing = !state.playing; playBtn.textContent = state.playing ? "Pause" : "Play"; });
    opts.appendChild(playBtn);
    ui.appendChild(opts);
    root.appendChild(ui);
    if (window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      state.playing = false; playBtn.textContent = "Play";
    }

    var panels = document.createElement("div");
    panels.className = "mt-fields__panels";
    function panel(title) {
      var f = document.createElement("figure"), c = document.createElement("canvas"), cap = document.createElement("figcaption");
      c.width = W; c.height = H; c.setAttribute("role", "img");
      cap.textContent = title; f.appendChild(c); f.appendChild(cap); panels.appendChild(f);
      return c;
    }
    var cP = panel("Pressure now");
    var cI = panel("Time-averaged intensity (1 = refocused matrix array at the target)");
    root.appendChild(panels);
    var readout = document.createElement("p");
    readout.className = "mt-fields__readout"; readout.setAttribute("aria-live", "polite");
    root.appendChild(readout);
    var text = document.createElement("p");
    text.className = "mt-fields__text";
    root.appendChild(text);

    var ctxP = cP.getContext("2d"), ctxI = cI.getContext("2d");
    var imgP = ctxP.createImageData(W, H), imgI = ctxI.createImageData(W, H);
    var U = null, NORM = 1, PEAK = null;

    function refresh() {
      precompute();
      stepSel.querySelectorAll("button").forEach(function (b) { b.setAttribute("aria-checked", b.dataset.step === state.step ? "true" : "false"); });
      typeSel.querySelectorAll("button").forEach(function (b) { b.setAttribute("aria-checked", b.dataset.type === state.type ? "true" : "false"); });
      var gr = groups(state.type), kind = state.step === "refocus" ? "refocus" : "planned";
      var Cn = drives(kind, gr);
      var mat = groups("matrix"), best = cabs(targetPressure(mat, drives("refocus", mat)));
      if (state.step === "reverse") {
        U = REV; NORM = cabs([REV.re[(TARGET[1] - 12) * W + TARGET[0]], REV.im[(TARGET[1] - 12) * W + TARGET[0]]]);
      } else {
        U = forwardField(gr, Cn);
        NORM = best;                         // the field at the target equals sum_e C_e G_e (reciprocity)
      }
      // intensity peak near the target (the focus)
      PEAK = null;
      if (state.step !== "reverse") {
        var pk = -1;
        for (var y = TARGET[1] - 30; y <= TARGET[1] + 30; y++) for (var x = TARGET[0] - 40; x <= TARGET[0] + 40; x++) {
          if (y < 0 || y >= H || SK[y * W + x]) continue;
          var j = y * W + x, v = U.re[j] * U.re[j] + U.im[j] * U.im[j];
          if (v > pk) { pk = v; PEAK = [x, y]; }
        }
      }
      var planned = cabs(targetPressure(gr, drives("planned", gr)));
      var nEl = gr.length, what = state.type === "annular" ? "rings" : "elements";
      if (state.step === "reverse") {
        readout.textContent = state.type === "single"
          ? "Arrival phases at the target spread over " + Math.round(spread(GP) * 180 / Math.PI) + "° across the one element's surface."
          : "Arrival phases at the target spread over " + Math.round(spread(transfers(gr)) * 180 / Math.PI) + "° across the " + nEl + " " + what + " (G_e).";
        text.textContent = STEPS.reverse.text + (state.type === "single" ? " A single element is one transducer surface: its colour varies along it, and no element phase can correct that." : "");
      } else {
        var gain = cabs(targetPressure(gr, Cn)) / planned;
        var sh = PEAK ? Math.hypot(PEAK[0] - TARGET[0], PEAK[1] - TARGET[1]) / LAMBDA : 0;
        readout.textContent = "Coherence at the target: " +
          (state.type === "single" ? "" : what + " " + coherenceElements(gr, Cn).toFixed(2) + ", ") +
          "aperture " + coherenceAperture(gr, Cn).toFixed(2) +
          " · target pressure relative to the planned drive: " + gain.toFixed(2) + "× (intensity " + (gain * gain).toFixed(2) +
          "×) · focus " + sh.toFixed(1) + " wavelengths from the target.";
        text.textContent = state.step === "planned" ? STEPS.planned.text : TYPES[state.type].text;
      }
      var di = imgI.data;
      for (var k = 0; k < W * H; k++) {
        var I = (U.re[k] * U.re[k] + U.im[k] * U.im[k]) / (NORM * NORM);
        var c = sequential(state.step === "reverse" ? Math.sqrt(I) * 0.5 : I);
        if (SK[k]) { c = [Math.round(c[0] * 0.6 + 60), Math.round(c[1] * 0.6 + 60), Math.round(c[2] * 0.6 + 60)]; }
        di[4 * k] = c[0]; di[4 * k + 1] = c[1]; di[4 * k + 2] = c[2]; di[4 * k + 3] = 255;
      }
      ctxI.putImageData(imgI, 0, 0);
      marks(ctxI);
      cP.setAttribute("aria-label", text.textContent);
    }

    function spread(list) {                   // range of arrival phases (relative to their mean)
      var m = [0, 0];
      list.forEach(function (g) { var a = Math.atan2(g[1], g[0]); m[0] += Math.cos(a); m[1] += Math.sin(a); });
      var mu = Math.atan2(m[1], m[0]), lo = 0, hi = 0;
      list.forEach(function (g) {
        var d = Math.atan2(g[1], g[0]) - mu; d = Math.atan2(Math.sin(d), Math.cos(d));
        lo = Math.min(lo, d); hi = Math.max(hi, d);
      });
      return hi - lo;
    }

    function marks(ctx) {
      var gr = groups(state.type), kind = state.step === "refocus" ? "refocus" : "planned";
      var Cn = drives(kind, gr), G = transfers(gr);
      if (state.type === "single") {
        // one bowl: a continuous arc, coloured point by point by arrival phase in the reverse run
        ctx.lineWidth = 4; ctx.lineCap = "round";
        for (var p = 0; p + 1 < PTS.length; p++) {
          ctx.strokeStyle = state.step === "reverse" ? hue(Math.atan2(GP[p][1], GP[p][0])) : "#555";
          ctx.beginPath(); ctx.moveTo(PTS[p][0], PTS[p][1]); ctx.lineTo(PTS[p + 1][0], PTS[p + 1][1]); ctx.stroke();
        }
        ctx.lineCap = "butt";
      } else {
        gr.forEach(function (pts, e) {
          var col = state.step === "reverse" ? hue(Math.atan2(G[e][1], G[e][0]))
                  : state.step === "refocus" ? hue(Math.atan2(Cn[e][1], Cn[e][0])) : "#555";
          var members = state.type === "annular" ? [e, N_EL - 1 - e] : [e];
          members.forEach(function (m) {
            var c = elemCentre(m);
            ctx.fillStyle = col; ctx.strokeStyle = "#222"; ctx.lineWidth = 0.6;
            ctx.beginPath(); ctx.arc(c[0], c[1], 3.2, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
          });
          if (state.type === "annular") {        // the two halves of one ring
            var a = elemCentre(e), b = elemCentre(N_EL - 1 - e);
            ctx.strokeStyle = "rgba(40,40,40,0.3)"; ctx.setLineDash([2, 2]); ctx.lineWidth = 0.8;
            ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
            ctx.setLineDash([]);
          }
        });
      }
      ctx.strokeStyle = "#111"; ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.arc(TARGET[0], TARGET[1], 4, 0, 2 * Math.PI); ctx.stroke();
      if (PEAK) {
        ctx.strokeStyle = "#0a7"; ctx.lineWidth = 1.5; ctx.beginPath();
        ctx.moveTo(PEAK[0] - 4, PEAK[1]); ctx.lineTo(PEAK[0] + 4, PEAK[1]);
        ctx.moveTo(PEAK[0], PEAK[1] - 4); ctx.lineTo(PEAK[0], PEAK[1] + 4); ctx.stroke();
      }
    }

    function drawFrame() {
      var d = imgP.data, w = 2 * Math.PI * state.t, cw = Math.cos(w), sw = Math.sin(w);
      for (var j = 0; j < W * H; j++) {
        var p = (U.re[j] * cw - U.im[j] * sw) / NORM;
        var c = diverging(state.step === "reverse" ? p * 0.6 : p);
        if (SK[j]) { c = [Math.round(c[0] * 0.6 + 60), Math.round(c[1] * 0.6 + 60), Math.round(c[2] * 0.6 + 60)]; }
        d[4 * j] = c[0]; d[4 * j + 1] = c[1]; d[4 * j + 2] = c[2]; d[4 * j + 3] = 255;
      }
      ctxP.putImageData(imgP, 0, 0);
      marks(ctxP);
    }

    var entry = { root: root, raf: null };
    function loop(ts) {
      if (!document.body.contains(root)) return;
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
    document.querySelectorAll(".tr-fields").forEach(build);
  }
  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(init);
  } else if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
