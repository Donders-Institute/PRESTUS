"""Analytic prediction of how several transducers' fields combine: coherent (complex pressures
add and interfere) or independent (intensities add, as async coupling assumes for heating).

A rough free-water picture to choose simulation.transducer_coupling, not a simulation:
- each transducer's focal spot is a Gaussian envelope (no side lobes, no skull) with the
  focused-piston focal widths of its aperture D and focal distance F (apex to focus):
  lateral intensity FWHM ≈ 1.02 λ F/D, axial ≈ 7 λ (F/D)²;
- inside the envelope a plane wave travels along the beam axis, so two beams crossing at
  angle θ make fringes λ / (2 sin(θ/2)) apart across their bisector;
- each spot's peak intensity is the transducer's target free-water ISPPA.

PRESTUS's coherent run starts every source at t = 0, so the relative phase of the beams at the
target follows from the path lengths and is not predicted here: the at-target numbers give the
whole range over relative phase, and the map shows in-phase arrival at the foci.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

C_WATER = 1500.0                 # m/s, free water
LATERAL_FWHM = 1.02              # × λ F/D   (intensity FWHM of a focused piston at its focus)
AXIAL_FWHM = 7.0                 # × λ (F/D)²


@dataclass
class Beam:
    apex: np.ndarray             # bowl apex, RAS mm
    focus: np.ndarray            # target, RAS mm
    aperture_mm: float
    freq_hz: float
    isppa_wcm2: float

    @property
    def axis(self) -> np.ndarray:
        d = np.asarray(self.focus, float) - np.asarray(self.apex, float)
        return d / np.linalg.norm(d)

    @property
    def focal_mm(self) -> float:
        return float(np.linalg.norm(np.asarray(self.focus, float) - np.asarray(self.apex, float)))

    def wavelength_mm(self, c: float = C_WATER) -> float:
        return c / self.freq_hz * 1e3

    def widths_mm(self, c: float = C_WATER) -> tuple[float, float]:
        """(lateral, axial) intensity FWHM of the focal spot [mm]."""
        lam, fnum = self.wavelength_mm(c), self.focal_mm / self.aperture_mm
        return LATERAL_FWHM * lam * fnum, AXIAL_FWHM * lam * fnum ** 2


def at_target(isppas: list[float]) -> dict:
    """Peak intensity where all foci coincide: independent sum, and the coherent range over
    the beams' relative phases (|Σ a_i e^{iφ_i}|², a_i = √I_i)."""
    a = np.sqrt(np.clip(np.asarray(isppas, float), 0, None))
    total = float(a.sum())
    return {"independent": float((a ** 2).sum()),
            "coherent_max": total ** 2,
            "coherent_min": max(0.0, float(a.max()) - (total - float(a.max()))) ** 2}


def beam_angle_deg(beams: list[Beam]) -> float:
    """Largest angle between two beam axes [deg]."""
    ang = 0.0
    for i in range(len(beams)):
        for j in range(i + 1, len(beams)):
            cos = float(np.clip(np.dot(beams[i].axis, beams[j].axis), -1, 1))
            ang = max(ang, math.degrees(math.acos(cos)))
    return ang


def fringe_spacing_mm(theta_deg: float, wavelength_mm: float) -> float | None:
    """Distance between interference maxima of two plane waves crossing at theta."""
    s = math.sin(math.radians(theta_deg) / 2)
    return wavelength_mm / (2 * s) if s > 1e-6 else None


def field_maps(beams: list[Beam], extent_mm: float = 20.0, n: int = 161,
               c: float = C_WATER) -> dict:
    """Intensity maps [W/cm²] in the plane through the mean focus spanned by the first two
    beam axes: 'independent' (Σ|p_i|²) and 'coherent' (|Σ p_i|², in phase at each focus),
    with in-plane coordinates u, v [mm] (u along the first beam)."""
    centre = np.mean([np.asarray(b.focus, float) for b in beams], axis=0)
    u = beams[0].axis
    other = next((b.axis for b in beams[1:] if abs(np.dot(b.axis, u)) < 0.9999), None)
    if other is None:                               # parallel beams: any perpendicular
        other = np.cross(u, [0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.cross(u, [1.0, 0.0, 0.0])
    v = other - np.dot(other, u) * u
    v /= np.linalg.norm(v)
    s = np.linspace(-extent_mm, extent_mm, n)
    uu, vv = np.meshgrid(s, s, indexing="xy")
    pts = centre + uu[..., None] * u + vv[..., None] * v
    ln2 = math.log(2)
    total = np.zeros(uu.shape, complex)
    power = np.zeros(uu.shape)
    for b in beams:
        lat_w, ax_w = b.widths_mm(c)
        rel = pts - np.asarray(b.focus, float)
        ax = rel @ b.axis
        lat = np.linalg.norm(rel - ax[..., None] * b.axis, axis=-1)
        env = np.exp(-2 * ln2 * ((lat / lat_w) ** 2 + (ax / ax_w) ** 2))   # pressure: I = env²
        p = math.sqrt(max(b.isppa_wcm2, 0.0)) * env * np.exp(1j * 2 * math.pi * ax / b.wavelength_mm(c))
        total += p
        power += np.abs(p) ** 2
    return {"u": s, "v": s, "independent": power, "coherent": np.abs(total) ** 2}
