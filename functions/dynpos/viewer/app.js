'use strict';

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
let nv            = null;   // NiiVue instance
let initData      = null;   // parsed placement_init.json
let transPosRAS   = [0, 0, 0];
let focusPosRAS   = [0, 0, 0];
let activeMarker  = 'trans';
let transducerMesh = null;  // current NiiVue mesh object for the bowl

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
function setStatus(msg) {
  document.getElementById('status').textContent = msg;
}

function dist3(a, b) {
  return Math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2);
}

function updateGeometrySummary() {
  const d = dist3(transPosRAS, focusPosRAS).toFixed(1);
  document.getElementById('g-dist').textContent = d + ' mm';
  if (initData && initData.transducer) {
    document.getElementById('g-diam').textContent =
      initData.transducer.diameter_mm.toFixed(1) + ' mm';
    document.getElementById('g-roc').textContent =
      initData.transducer.curv_radius_mm.toFixed(1) + ' mm';
  }
}

function readCoordFields() {
  transPosRAS = [
    parseFloat(document.getElementById('tx').value) || 0,
    parseFloat(document.getElementById('ty').value) || 0,
    parseFloat(document.getElementById('tz').value) || 0,
  ];
  focusPosRAS = [
    parseFloat(document.getElementById('fx').value) || 0,
    parseFloat(document.getElementById('fy').value) || 0,
    parseFloat(document.getElementById('fz').value) || 0,
  ];
}

function writeCoordFields() {
  document.getElementById('tx').value = transPosRAS[0].toFixed(2);
  document.getElementById('ty').value = transPosRAS[1].toFixed(2);
  document.getElementById('tz').value = transPosRAS[2].toFixed(2);
  document.getElementById('fx').value = focusPosRAS[0].toFixed(2);
  document.getElementById('fy').value = focusPosRAS[1].toFixed(2);
  document.getElementById('fz').value = focusPosRAS[2].toFixed(2);
}

// ---------------------------------------------------------------------------
// Transducer mesh: spherical-cap bowl + cylindrical housing
//
// PRESTUS coordinate conventions:
//   trans_pos  — apex (centre-back) of the concave bowl surface
//   focus_pos  — acoustic target; defines beam-axis direction n̂
//   R          — curv_radius_mm: radius of sphere the bowl is carved from
//   D          — diameter_mm:    outer aperture diameter
//   depth_mm   — axial depth of the housing cylinder behind the exit plane
//
// Sphere centre: C = trans_pos + R·n̂  (matches PRESTUS get_transducer_box:
//   geom_focus_pos = trans_pos − R·focal_slope,  focal_slope = −n̂)
//
// Bowl vertex at (θ, φ):
//   p = trans_pos + R·(1−cosθ)·n̂ + R·sinθ·(cosφ·û + sinφ·v̂)
//   θ ∈ [0, α_max],  α_max = arcsin(D/2R)
//
// Exit-plane offset from apex (sagitta):
//   sagitta = R − √(R²−(D/2)²)   (PRESTUS: R − dist_to_ep_mm)
//
// Housing cylinder: exit-plane ring → back ring, both at radius D/2.
//   exit-plane at trans + sagitta·n̂  (matches bowl rim)
//   back ring   at trans + (sagitta − depth_mm)·n̂
// ---------------------------------------------------------------------------

function buildTransducerSTL(trans, focus, diameter_mm, curv_radius_mm, depth_mm) {
  const dx = focus[0] - trans[0];
  const dy = focus[1] - trans[1];
  const dz = focus[2] - trans[2];
  const axisLen = Math.sqrt(dx*dx + dy*dy + dz*dz);
  if (axisLen < 1e-6 || curv_radius_mm <= 0 || diameter_mm <= 0) return null;

  const R     = curv_radius_mm;
  const halfD = diameter_mm / 2;
  if (halfD >= R) return null;  // degenerate
  const alphaMax = Math.asin(Math.min(halfD / R, 0.9998));
  // dist_to_ep_mm in PRESTUS = sqrt(R²−(D/2)²)
  // ex_plane_pos = (trans + R·n̂) − dist_to_ep·n̂ = trans + (R − dist_to_ep)·n̂
  // so the axial offset from apex to exit plane is the sagitta: R − sqrt(R²−(D/2)²)
  const dist_to_ep = Math.sqrt(R*R - halfD*halfD);
  const sagitta    = R - dist_to_ep;            // apex → exit plane (a few mm)

  // Beam-axis unit vector n̂ (trans → focus)
  const nx = dx / axisLen, ny = dy / axisLen, nz = dz / axisLen;

  // Orthonormal basis (û, v̂) perpendicular to n̂
  let ux, uy, uz;
  if (Math.abs(nx) < 0.9) { ux = 0;   uy = nz;  uz = -ny; }
  else                     { ux = -nz; uy = 0;   uz = nx;  }
  const uLen = Math.sqrt(ux*ux + uy*uy + uz*uz);
  ux /= uLen; uy /= uLen; uz /= uLen;
  const vx = ny*uz - nz*uy;
  const vy = nz*ux - nx*uz;
  const vz = nx*uy - ny*ux;

  const nRings = 12;
  const nSeg   = 24;

  // ── vertices ─────────────────────────────────────────────────────────────
  const verts = [];

  // [0] apex
  verts.push([trans[0], trans[1], trans[2]]);

  // [1 … nRings*nSeg] bowl rings
  for (let ri = 1; ri <= nRings; ri++) {
    const theta = (ri / nRings) * alphaMax;
    const sinT  = Math.sin(theta), cosT = Math.cos(theta);
    for (let si = 0; si < nSeg; si++) {
      const phi  = (si / nSeg) * 2 * Math.PI;
      const cosP = Math.cos(phi), sinP = Math.sin(phi);
      verts.push([
        trans[0] + R*(1-cosT)*nx + R*sinT*(cosP*ux + sinP*vx),
        trans[1] + R*(1-cosT)*ny + R*sinT*(cosP*uy + sinP*vy),
        trans[2] + R*(1-cosT)*nz + R*sinT*(cosP*uz + sinP*vz),
      ]);
    }
  }

  // Housing cylinder (only when depth_mm is meaningful)
  const addHousing = (depth_mm > 0);
  const iCylTop  = verts.length;    // exit-plane ring start index
  const iCylBot  = addHousing ? iCylTop + nSeg : -1;
  const iBackCtr = addHousing ? iCylTop + 2*nSeg : -1;

  if (addHousing) {
    // Exit-plane ring: at sagitta offset from apex, matching the bowl rim
    for (let si = 0; si < nSeg; si++) {
      const phi  = (si / nSeg) * 2 * Math.PI;
      const cosP = Math.cos(phi), sinP = Math.sin(phi);
      verts.push([
        trans[0] + sagitta*nx + halfD*(cosP*ux + sinP*vx),
        trans[1] + sagitta*ny + halfD*(cosP*uy + sinP*vy),
        trans[2] + sagitta*nz + halfD*(cosP*uz + sinP*vz),
      ]);
    }
    // Back ring: exit plane shifted by depth_mm away from focus (−n̂)
    const backOff = sagitta - depth_mm;
    for (let si = 0; si < nSeg; si++) {
      const phi  = (si / nSeg) * 2 * Math.PI;
      const cosP = Math.cos(phi), sinP = Math.sin(phi);
      verts.push([
        trans[0] + backOff*nx + halfD*(cosP*ux + sinP*vx),
        trans[1] + backOff*ny + halfD*(cosP*uy + sinP*vy),
        trans[2] + backOff*nz + halfD*(cosP*uz + sinP*vz),
      ]);
    }
    // Back-disk centre
    verts.push([trans[0] + backOff*nx, trans[1] + backOff*ny, trans[2] + backOff*nz]);
  }

  // ── triangles ─────────────────────────────────────────────────────────────
  const tris = [];

  // Bowl: apex fan → first ring
  for (let si = 0; si < nSeg; si++) {
    tris.push([0, 1 + si, 1 + (si + 1) % nSeg]);
  }
  // Bowl: quad strips between rings
  for (let ri = 0; ri < nRings - 1; ri++) {
    const b1 = 1 + ri * nSeg, b2 = 1 + (ri + 1) * nSeg;
    for (let si = 0; si < nSeg; si++) {
      const a = b1 + si, b = b1 + (si+1)%nSeg;
      const c = b2 + (si+1)%nSeg, d = b2 + si;
      tris.push([a, b, c]);
      tris.push([a, c, d]);
    }
  }

  if (addHousing) {
    // Cylinder side wall
    for (let si = 0; si < nSeg; si++) {
      const a = iCylTop + si,           b = iCylTop + (si+1)%nSeg;
      const c = iCylBot + (si+1)%nSeg,  d = iCylBot + si;
      tris.push([a, b, c]);
      tris.push([a, c, d]);
    }
    // Back disk (fan from centre, wound so normal faces away from focus)
    for (let si = 0; si < nSeg; si++) {
      tris.push([iBackCtr, iCylBot + (si+1)%nSeg, iCylBot + si]);
    }
  }

  // ── ASCII STL ─────────────────────────────────────────────────────────────
  const lines = ['solid transducer'];
  for (const [i0, i1, i2] of tris) {
    const [x0,y0,z0] = verts[i0], [x1,y1,z1] = verts[i1], [x2,y2,z2] = verts[i2];
    const ex1 = x1-x0, ey1 = y1-y0, ez1 = z1-z0;
    const ex2 = x2-x0, ey2 = y2-y0, ez2 = z2-z0;
    let fnx = ey1*ez2 - ez1*ey2, fny = ez1*ex2 - ex1*ez2, fnz = ex1*ey2 - ey1*ex2;
    const fnL = Math.sqrt(fnx*fnx + fny*fny + fnz*fnz) || 1;
    fnx /= fnL; fny /= fnL; fnz /= fnL;
    lines.push(`facet normal ${fnx.toFixed(6)} ${fny.toFixed(6)} ${fnz.toFixed(6)}`);
    lines.push('  outer loop');
    lines.push(`    vertex ${x0.toFixed(4)} ${y0.toFixed(4)} ${z0.toFixed(4)}`);
    lines.push(`    vertex ${x1.toFixed(4)} ${y1.toFixed(4)} ${z1.toFixed(4)}`);
    lines.push(`    vertex ${x2.toFixed(4)} ${y2.toFixed(4)} ${z2.toFixed(4)}`);
    lines.push('  endloop');
    lines.push('endfacet');
  }
  lines.push('endsolid transducer');
  return lines.join('\n');
}

async function refreshTransducerMesh() {
  if (!nv || !initData?.transducer) return;
  const { diameter_mm, curv_radius_mm, depth_mm = 0 } = initData.transducer;
  if (!(diameter_mm > 0) || !(curv_radius_mm > 0)) return;

  // Remove previous mesh
  if (transducerMesh !== null) {
    try { nv.removeMesh(transducerMesh); } catch (_) {}
    transducerMesh = null;
  }

  const stl = buildTransducerSTL(transPosRAS, focusPosRAS, diameter_mm, curv_radius_mm, depth_mm);
  if (!stl) return;

  const blob = new Blob([stl], { type: 'text/plain' });
  const url  = URL.createObjectURL(blob);
  try {
    // Pass name with .stl extension so NiiVue picks the right parser
    await nv.loadMeshes([{
      url,
      name:    'transducer.stl',
      rgba255: [220, 80, 80, 180],
      opacity: 0.65,
    }]);
    transducerMesh = nv.meshes[nv.meshes.length - 1];
  } catch (_) {
    // Mesh overlay not critical — silently skip on unsupported NiiVue builds
  } finally {
    URL.revokeObjectURL(url);
  }
}

// ---------------------------------------------------------------------------
// NiiVue annotation — two crosshairs rendered as scene annotations
// ---------------------------------------------------------------------------
function refreshAnnotations() {
  if (!nv) return;

  // NiiVue ≥ 0.43 supports drawScene annotations via nv.document.scene
  // We use addMesh-style spheres via the public "drawCrosshairs3D" path,
  // falling back to moving the crosshair cursor for older builds.

  try {
    // Clear previous custom drawings
    if (typeof nv.setDrawingEnabled === 'function') {
      nv.setDrawingEnabled(false);
    }

    // Mark the transducer position by jumping the crosshair there
    // (navigateTo expects mm in world/RAS space)
    const pos = activeMarker === 'trans' ? transPosRAS : focusPosRAS;
    if (typeof nv.scene !== 'undefined') {
      nv.scene.crosshairPos = nv.mm2frac(pos);
      nv.drawScene();
    }
  } catch (e) {
    // Silently ignore annotation errors — positions are still editable
    // via the coordinate fields.
  }

  updateGeometrySummary();
  refreshTransducerMesh(); // fire-and-forget; errors are caught internally
}

// ---------------------------------------------------------------------------
// NiiVue click handler — place the active marker at the clicked voxel
// ---------------------------------------------------------------------------
function onViewerClick(e) {
  if (!nv) return;

  const canvas = document.getElementById('gl1');
  const rect   = canvas.getBoundingClientRect();
  const x = e.clientX - rect.left;
  const y = e.clientY - rect.top;

  // NiiVue exposes the clicked mm position via the crosshair after a
  // mouse event; we can read it back from scene.crosshairPos → mm.
  // We trigger the internal pick by forwarding the raw event to NiiVue's
  // canvas handler, then read back crosshairPos.
  try {
    const frac = nv.canvasPos2frac([x, y]);
    if (!frac || frac.some(isNaN)) return;
    const mm = nv.frac2mm(frac);

    if (activeMarker === 'trans') {
      transPosRAS = Array.from(mm);
    } else {
      focusPosRAS = Array.from(mm);
    }

    writeCoordFields();
    refreshAnnotations();
  } catch (e2) {
    // canvasPos2frac may not exist in older NiiVue — ignore
  }
}

// ---------------------------------------------------------------------------
// Coordinate field change handlers
// ---------------------------------------------------------------------------
function bindCoordInputs() {
  ['tx','ty','tz','fx','fy','fz'].forEach(id => {
    document.getElementById(id).addEventListener('change', () => {
      readCoordFields();
      refreshAnnotations();
    });
  });
}

// ---------------------------------------------------------------------------
// Save / cancel
// ---------------------------------------------------------------------------
async function saveAndClose() {
  readCoordFields();
  const payload = {
    trans_pos_ras: transPosRAS,
    focus_pos_ras: focusPosRAS,
  };
  setStatus('Saving …');
  try {
    const resp = await fetch('/save_placement', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(payload),
    });
    if (resp.ok) {
      setStatus('Saved. You can close this tab.');
      document.getElementById('btn-save').disabled   = true;
      document.getElementById('btn-cancel').disabled = true;
    } else {
      setStatus('Server error — check MATLAB console.');
    }
  } catch (e) {
    setStatus('Could not reach server: ' + e.message);
  }
}

async function cancelAlignment() {
  setStatus('Cancelling …');
  try {
    await fetch('/cancel_placement', { method: 'POST' });
    setStatus('Cancelled. You can close this tab.');
  } catch (e) {
    setStatus('Could not reach server.');
  }
}

function resetToInitial() {
  if (!initData) return;
  transPosRAS = [...initData.trans_pos_ras];
  focusPosRAS = [...initData.focus_pos_ras];
  writeCoordFields();
  refreshAnnotations();
  setStatus('Reset to initial positions.');
}

// ---------------------------------------------------------------------------
// Initialise NiiVue
// ---------------------------------------------------------------------------
async function initViewer(niiFilename) {
  nv = new niivue.Niivue({
    backColor:          [0.1, 0.1, 0.1, 1],
    crosshairColor:     [1, 0.2, 0.2, 0.9],   // red → transducer side
    show3Dcrosshair:    true,
    multiplanarLayout:  niivue.MULTIPLANAR_TYPE ? niivue.MULTIPLANAR_TYPE.GRID : 2,
  });

  await nv.attachToCanvas(document.getElementById('gl1'));

  await nv.loadVolumes([{
    url:         niiFilename,
    colormap:    'gray',
    opacity:     1.0,
    visible:     true,
  }]);

  // Jump crosshair to transducer entry point
  try {
    nv.scene.crosshairPos = nv.mm2frac(transPosRAS);
    nv.drawScene();
  } catch (_) {}

  // Draw initial bowl
  await refreshTransducerMesh();

  // Wire up click-to-place
  const canvas = document.getElementById('gl1');
  canvas.addEventListener('click', onViewerClick);
}

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
async function boot() {
  setStatus('Loading …');

  // Fetch initial placement
  let resp;
  try {
    resp = await fetch('/placement_init.json');
    if (!resp.ok) throw new Error('HTTP ' + resp.status);
    initData = await resp.json();
  } catch (e) {
    setStatus('Failed to load placement_init.json: ' + e.message);
    return;
  }

  transPosRAS = [...initData.trans_pos_ras];
  focusPosRAS = [...initData.focus_pos_ras];

  writeCoordFields();
  updateGeometrySummary();

  // Active marker radio
  document.querySelectorAll('input[name=active_marker]').forEach(el => {
    el.addEventListener('change', e => {
      activeMarker = e.target.value;
      refreshAnnotations();
    });
  });

  bindCoordInputs();
  document.getElementById('btn-save').addEventListener('click', saveAndClose);
  document.getElementById('btn-cancel').addEventListener('click', cancelAlignment);
  document.getElementById('btn-reset').addEventListener('click', resetToInitial);

  // Load NIfTI then draw the bowl
  try {
    await initViewer(initData.nii_filename);
    setStatus('Ready. Click to place the active marker.');
  } catch (e) {
    setStatus('Viewer error: ' + e.message);
    console.error(e);
  }
}

document.addEventListener('DOMContentLoaded', boot);
