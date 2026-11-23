// Interactive illustration for doc_simulations-thermal.md (Pulse Model): the heat
// source and the temperature at one point for thermal.pulse_model 'pulse_resolved'
// (every pulse and gap a segment) and 'train_averaged' (each pulse train one
// segment, source x PD/PRI), with the segment counts each needs (THERMAL_SCHEDULE).
// Temperature: one point with one thermal time constant tau, dT/dt = q - T/tau,
// a schematic stand-in for diffusion and perfusion; values are relative to the
// pulse-resolved peak. Renders into every <div class="thermal-pulse-model"></div>.
(function () {
  "use strict";

  var C_RES = "#8c959f", C_AVG = "#0e7490";
  var W = 640, H = 300, ML = 46, MR = 12, MT = 14, GAP = 34;
  var PH = (H - MT - GAP - 26) / 2;          // panel height

  var DEFAULTS = { pd: 50, pri: 200, ptd: 40, ptri: 50, ptrd: 150, tau: 10, post: 1 };
  var FIELDS = [
    ["pd", "PD [ms]"], ["pri", "PRI [ms]"], ["ptd", "PTD [s]"], ["ptri", "PTRI [s]"],
    ["ptrd", "PTRD [s]"], ["tau", "tissue time constant τ [s]"], ["post", "post_pt_timestep [s]"]
  ];

  function el(tag, attrs, parent) {
    var ns = /^(svg|path|line|text|rect|g|polyline)$/.test(tag) ? "http://www.w3.org/2000/svg" : null;
    var e = ns ? document.createElementNS(ns, tag) : document.createElement(tag);
    for (var k in attrs) { if (k === "text") e.textContent = attrs[k]; else e.setAttribute(k, attrs[k]); }
    if (parent) parent.appendChild(e);
    return e;
  }

  function simulate(v) {
    var pd = v.pd / 1e3, pri = v.pri / 1e3, ptd = v.ptd, ptri = v.ptri, ptrd = v.ptrd, tau = v.tau;
    var nTrain = Math.max(1, Math.round(ptrd / ptri));
    var nPulse = Math.max(1, Math.floor(ptd / pri + 1e-9));
    var duty = pd / pri;
    var dt = Math.min(pd, pri - pd, ptd) / 6;
    var n = Math.min(400000, Math.ceil(ptrd / dt));
    dt = ptrd / n;
    var stride = Math.max(1, Math.floor(n / 1500));
    var t = [], qr = [], qa = [], Tr = [], Ta = [];
    var tr = 0, ta = 0, peakR = 0, peakA = 0;
    for (var i = 0; i <= n; i++) {
      var ti = i * dt, k = Math.floor(ti / ptri), u = ti - k * ptri;
      var inTrain = k < nTrain && u < Math.min(ptd, nPulse * pri);
      var qres = inTrain && (u % pri) < pd ? 1 : 0;
      var qavg = inTrain ? duty : 0;
      tr += dt * (qres - tr / tau);
      ta += dt * (qavg - ta / tau);
      peakR = Math.max(peakR, tr); peakA = Math.max(peakA, ta);
      if (i % stride === 0 || i === n) { t.push(ti); qr.push(qres); qa.push(qavg); Tr.push(tr); Ta.push(ta); }
    }
    var gapSeg = Math.max(0, nTrain - 1) * Math.max(1, Math.round(Math.max(0, ptri - ptd) / v.post));
    return {
      t: t, qr: qr, qa: qa, Tr: Tr, Ta: Ta, peakR: peakR, peakA: peakA, endR: tr, endA: ta, duty: duty,
      segRes: nTrain * nPulse * 2 + gapSeg,
      segAvg: nTrain * Math.max(1, Math.round(ptd / v.post)) + gapSeg
    };
  }

  function draw(svg, s, v) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var x0 = ML, x1 = W - MR, tmax = v.ptrd;
    var X = function (t) { return x0 + (x1 - x0) * t / tmax; };
    var panel = function (top, label, yMax, series) {
      var Y = function (y) { return top + PH - PH * y / yMax; };
      el("line", { x1: x0, x2: x1, y1: top + PH, y2: top + PH, class: "tpm-axis" }, svg);
      el("line", { x1: x0, x2: x0, y1: top, y2: top + PH, class: "tpm-axis" }, svg);
      el("text", { x: 6, y: top + PH / 2, class: "tpm-label", transform: "rotate(-90 12 " + (top + PH / 2) + ")", "text-anchor": "middle", text: label }, svg);
      series.forEach(function (sr) {
        var d = "";
        for (var i = 0; i < s.t.length; i++) d += (i ? "L" : "M") + X(s.t[i]).toFixed(1) + " " + Y(sr.y[i]).toFixed(1);
        el("path", { d: d, fill: "none", stroke: sr.c, "stroke-width": sr.w, opacity: sr.o }, svg);
      });
    };
    panel(MT, "heat source", 1.05, [
      { y: s.qr, c: C_RES, w: 0.6, o: 0.8 },
      { y: s.qa, c: C_AVG, w: 2, o: 1 }]);
    var top2 = MT + PH + GAP, m = Math.max(s.peakR, s.peakA) * 1.08;
    panel(top2, "temperature rise", m, [
      { y: s.Tr, c: C_RES, w: 1, o: 0.9 },
      { y: s.Ta, c: C_AVG, w: 2, o: 1 }]);
    for (var k = 0; k <= 5; k++) {
      var tt = tmax * k / 5;
      el("text", { x: X(tt), y: top2 + PH + 14, class: "tpm-tick", "text-anchor": "middle", text: (tt >= 10 ? tt.toFixed(0) : tt.toFixed(1)) + " s" }, svg);
    }
  }

  function render(root) {
    var v = {}; for (var k in DEFAULTS) v[k] = DEFAULTS[k];
    var form = el("div", { class: "tpm-controls" }, root);
    FIELDS.forEach(function (f) {
      var lab = el("label", {}, form);
      el("span", { text: f[1] }, lab);
      var inp = el("input", { type: "number", min: "0", step: "any", value: v[f[0]] }, lab);
      inp.addEventListener("input", function () {
        var x = parseFloat(inp.value);
        if (isFinite(x) && x > 0) { v[f[0]] = x; update(); }
      });
    });
    var legend = el("div", { class: "tpm-legend" }, root);
    legend.innerHTML = '<span><i style="background:' + C_RES + '"></i>pulse_resolved</span>' +
      '<span><i style="background:' + C_AVG + '"></i>train_averaged (default)</span>';
    var svg = el("svg", { viewBox: "0 0 " + W + " " + H, class: "tpm-svg", role: "img",
      "aria-label": "Heat source and temperature rise over time for the two pulse models" }, root);
    var out = el("div", { class: "tpm-out" }, root);

    function update() {
      if (v.pd >= v.pri || v.ptd > v.ptri || v.ptri > v.ptrd) {
        out.textContent = "Needs PD < PRI and PTD ≤ PTRI ≤ PTRD.";
        return;
      }
      var s = simulate(v);
      draw(svg, s, v);
      var dPeak = 100 * (s.peakA / s.peakR - 1), dEnd = 100 * (s.endA / s.endR - 1);
      out.innerHTML =
        "Duty cycle in a train <b>" + (100 * s.duty).toFixed(1) + " %</b> · peak rise, train-averaged vs pulse-resolved: <b>" +
        (dPeak >= 0 ? "+" : "") + dPeak.toFixed(1) + " %</b> (end of protocol " + (dEnd >= 0 ? "+" : "") + dEnd.toFixed(1) +
        " %) · heating segments: <b>" + s.segRes.toLocaleString() + "</b> pulse-resolved, <b>" + s.segAvg.toLocaleString() +
        "</b> train-averaged";
    }
    update();
  }

  function init() {
    var roots = document.querySelectorAll("div.thermal-pulse-model");
    for (var i = 0; i < roots.length; i++) {
      if (!roots[i].dataset.done) { roots[i].dataset.done = "1"; render(roots[i]); }
    }
  }
  if (window.document$ && window.document$.subscribe) window.document$.subscribe(init);   // mkdocs-material instant loading
  else if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
