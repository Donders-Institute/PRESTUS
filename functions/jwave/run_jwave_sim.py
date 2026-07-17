#!/usr/bin/env python3
"""
run_jwave_sim.py — PRESTUS → jwave acoustic simulation bridge

Called by jwave_wrapper.m via system(). Reads simulation inputs from an HDF5
file written by MATLAB, runs the requested jwave solver, and writes pressure
field outputs (p_max_all, p_min_all, p_final) to a second HDF5 file that
MATLAB reads back.

Usage:
    python run_jwave_sim.py \\
        --input  /path/to/sub-001_layered_jwave_input.h5 \\
        --output /path/to/sub-001_layered_jwave_output.h5 \\
        --solver pseudo_spectral

Solvers:
    pseudo_spectral  - Pseudo-spectral time-domain (PSTD) propagation via
                       jwave.acoustics.time_varying.simulate_wave_propagation.
                       Uses the Stokes (y=2) pressure-correction absorber:
                           p += c₀²·absorb_tau·ρ₀·∇·u   (absorb_tau < 0, dissipative)
                       absorb_tau is rescaled via _absorb_tau_from_power_law so
                       that att(freq_hz) matches the input (alpha_coeff,
                       alpha_power) exactly.  Off-frequency attenuation still
                       follows the Stokes f² law.

    finite_differences - Same Stokes-absorbing time-domain step as
                         pseudo_spectral but with medium fields wrapped as
                         JaxDF FiniteDifferences, so jwave's operator dispatch
                         selects local FD stencils (order --fd_accuracy) rather
                         than global FFT derivatives. Useful for benchmarking
                         numerical dispersion vs PSTD.

    helmholtz        - Frequency-domain GMRES solver via
                       jwave.acoustics.time_harmonic.helmholtz_solver.
                       Solves the heterogeneous Helmholtz equation with a
                       Fourier spectral collocation method and PML absorbing
                       boundaries. Attenuation applied via db2neper(y=2.0).
                       Requires freq_hz in the HDF5 metadata attribute.
                       alpha_coeff is floored to 0.002 dB/(cm.MHz^2) (see
                       _apply_helmholtz_attenuation_floor) to regularize GMRES
                       for lossless/near-lossless media (e.g. water benchmarks).

Exit codes:
    0  — success
    1  — error (message printed to stderr)

Diagnostics written to /diagnostics/ in the output HDF5:
    residual_norm  — exact relative GMRES residual ||A·u - b||/||b|| (Helmholtz
                     only), evaluated with jwave's own helmholtz operator and
                     scale_source_helmholtz (see compute_helmholtz_residual) —
                     not a hand-rolled approximation.
"""

import argparse
import logging
import os
import sys
import numpy as np
import h5py
import jax
jax.config.update("jax_compilation_cache_dir", os.path.join(os.path.expanduser("~"), ".cache", "jax_xla"))

# Suppress jaxdf deprecation warnings (e.g. shift_operator multi-output)
logging.getLogger('jaxdf').setLevel(logging.ERROR)


def parse_args():
    p = argparse.ArgumentParser(description='PRESTUS–jwave bridge')
    p.add_argument('--input',  required=True, help='Path to input HDF5 from MATLAB')
    p.add_argument('--output', required=True, help='Path to write output HDF5')
    p.add_argument('--solver', default='helmholtz',
                   choices=['helmholtz', 'pseudo_spectral', 'finite_differences'],
                   help='Solver backend (default: helmholtz)')
    p.add_argument('--gmres_restart', type=int, default=10,
                   help='GMRES restart parameter (default: 10, matches jwave\'s own '
                        'helmholtz_solver default; Helmholtz only)')
    p.add_argument('--gmres_maxiter', type=int, default=1000,
                   help='Max GMRES outer iterations (default: 1000, matches jwave\'s own '
                        'helmholtz_solver default; Helmholtz only)')
    p.add_argument('--gmres_tol', type=float, default=1e-3,
                   help='GMRES convergence tolerance (default: 1e-3, matches jwave\'s own '
                        'helmholtz_solver default; Helmholtz only)')
    p.add_argument('--debug', action='store_true',
                   help='Print extra diagnostics: PML/domain sizing in wavelengths, '
                        'steady-state timing margin (time-domain solvers), and a '
                        'compile-vs-solve time breakdown (Helmholtz).')
    p.add_argument('--fd_accuracy', type=int, default=8,
                   help='FD stencil accuracy order for finite_differences solver (default: 8)')
    return p.parse_args()


# ---------------------------------------------------------------------------
# HDF5 I/O helpers
# ---------------------------------------------------------------------------

def read_input(h5_path):
    """Read all simulation inputs from the MATLAB-written HDF5 file."""
    with h5py.File(h5_path, 'r') as f:
        N  = f['/grid/N'][()].flatten().astype(int)
        dx = f['/grid/dx'][()].flatten().astype(float)

        Nt = int(np.asarray(f['/time/Nt'][()]).flat[0])
        dt = float(np.asarray(f['/time/dt'][()]).flat[0])

        def read_field(key):
            v = f[key][()]
            if v.size == 1:
                return float(v)
            # MATLAB writes arrays in column-major order; HDF5/h5py reverses
            # dimension order, so a MATLAB [Nx,Ny,Nz] array appears as (Nz,Ny,Nx).
            # Transpose restores the logical layout expected by jwave.
            if v.ndim >= 2:
                return v.T
            return v

        sound_speed = read_field('/medium/sound_speed')
        density     = read_field('/medium/density')
        alpha_coeff = read_field('/medium/alpha_coeff')
        alpha_power = read_field('/medium/alpha_power')

        mask_indices        = f['/source/mask_indices'][()].flatten().astype(int)
        # MATLAB writes source.p as [N_unique x Nt] in column-major; h5py reads
        # it transposed as (Nt, N_unique). Transpose back to (N_unique, Nt).
        source_p            = f['/source/p'][()].T        # [N_unique x Nt]
        # Expand deduplicated signals: sig_map[i] is the index into source_p for voxel i.
        if '/source/sig_map' in f:
            sig_map  = f['/source/sig_map'][()].flatten().astype(int)
            source_p = source_p[sig_map, :]               # [N_src x Nt]
        record_start_index  = int(np.asarray(f['/source/record_start_index'][()]).flat[0])

        # Read as array to support per-axis vectors written by jwave_wrapper.
        # Files written by older PRESTUS versions store a scalar [1,1]; broadcast
        # that to a 1-D array of length ndims for uniform downstream handling.
        pml_raw  = np.asarray(f['/pml/size'][()]).flatten().astype(int)
        pml_size = pml_raw if pml_raw.size > 1 else pml_raw  # always ndarray

        freq_hz  = float(np.asarray(f['/'].attrs.get('freq_hz', 500e3)).flat[0])
        dims     = int(np.asarray(f['/'].attrs.get('dims', len(N))).flat[0])

    # Validate source dimensions immediately after read
    n_src_mask = len(mask_indices)
    n_src_sig  = source_p.shape[0]
    if n_src_mask != n_src_sig:
        raise ValueError(
            f'Source shape mismatch: mask has {n_src_mask} active voxels '
            f'but source.p has {n_src_sig} rows (expected {n_src_mask}). '
            'Check that source.p was built from source.p_mask with find().'
        )
    pml_str = str(pml_size.tolist()) if pml_size.size > 1 else str(int(pml_size[0]))
    print(f'[run_jwave_sim] N={N.tolist()}  Nt={Nt}  dt={dt:.3e} s  '
          f'N_src={n_src_mask}  freq={freq_hz/1e3:.1f} kHz  pml={pml_str}')

    return dict(
        N=N, dx=dx, Nt=Nt, dt=dt,
        sound_speed=sound_speed, density=density,
        alpha_coeff=alpha_coeff, alpha_power=alpha_power,
        mask_indices=mask_indices, source_p=source_p,
        record_start_index=record_start_index,
        pml_size=pml_size, freq_hz=freq_hz, dims=dims,
    )


def write_output(h5_path, p_max_all, p_min_all, p_final, diagnostics=None):
    """Write pressure field outputs for MATLAB to read.

    Arrays are transposed before writing so that MATLAB reads them with the
    correct (Nx, Ny, Nz) logical layout (HDF5/h5py reverses axis order relative
    to MATLAB's column-major convention).

    diagnostics: optional dict of scalar floats/ints written to /diagnostics/.
    """
    def to_matlab(arr):
        return np.ascontiguousarray(arr.T).astype(np.float32)

    gz = dict(compression='gzip', compression_opts=4)
    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('/pressure/p_max_all', data=to_matlab(p_max_all), **gz)
        f.create_dataset('/pressure/p_min_all', data=to_matlab(p_min_all), **gz)
        f.create_dataset('/pressure/p_final',   data=to_matlab(p_final),   **gz)
        if diagnostics:
            for key, val in diagnostics.items():
                f.create_dataset(f'/diagnostics/{key}', data=np.float64(val))


# ---------------------------------------------------------------------------
# Medium / geometry helpers
# ---------------------------------------------------------------------------

def auto_pad_for_pml(inp, stencil_clearance=0):
    """Pad the simulation grid so source elements lie outside the PML zone.

    k-Wave adds PML *outside* the specified grid, so a transducer placed at
    z=2 (1-based) is safe.  jwave applies PML *within* the grid, so the same
    transducer falls inside the PML and its energy is absorbed before it can
    propagate.  This function detects that situation and pads all medium arrays
    + shifts source positions so the source is at least pml_size voxels from
    every boundary — matching k-Wave's effective geometry.

    stencil_clearance: extra voxels beyond pml_size required for the numerical
    stencil.  For FD accuracy order N, pass N//2 so the stencil never reaches
    into the PML (e.g. accuracy=8 → stencil_clearance=4).  When > 0 this value
    is also applied as a *minimum symmetric pad on every axis*, independent of
    source position: jaxdf's FD stencil convolution pads periodically
    (mode="wrap"), so a half-width-N//2 stencil evaluated at a true array edge
    reads across the seam into the opposite face.  Guaranteeing stencil_clearance
    voxels of (edge-replicated, PML-absorbed) buffer on both sides of every axis
    keeps that seam and the outermost stencil taps out of the live medium.
    The added cost is only 2*stencil_clearance voxels per axis (e.g. 8 for
    accuracy=8) — negligible versus growing pml_size.

    Note: this guards against a real but secondary risk (stencil taps sampling
    across the periodic seam).  It was suspected as the primary cause of a
    persistent late-onset NaN on small, PML-dominated pFOV domains, but that
    NaN was actually caused by a sign error in _absorb_tau_from_power_law
    (fixed; see that function's docstring) — the Stokes absorption term was
    energy-injecting rather than dissipative, which no amount of PML/stencil
    padding or dt tuning could have fixed.  This padding is kept as legitimate
    defense-in-depth, not because it was the fix for that bug.

    pml_size may be a scalar int or a 1-D integer array of length ndims.
    When it is an array, each axis uses its own threshold so anisotropic
    PML configurations are handled correctly.

    Returns (padded_inp, pad_widths) where pad_widths is a list of (lo, hi)
    tuples per axis (in MATLAB/Fortran axis order). pad_widths is all zeros
    when no padding was needed.
    """
    pml_raw  = inp['pml_size']
    N    = inp['N'].tolist()          # MATLAB axis order [Nx, Ny, Nz]
    ndim = len(N)
    # Build per-axis safe-pml thresholds (scalar → broadcast).
    pml_arr = np.asarray(pml_raw, dtype=int).flatten()
    if pml_arr.size == 1:
        pml_arr = np.full(ndim, int(pml_arr[0]), dtype=int)
    safe_pml = pml_arr[:ndim] + stencil_clearance   # per-axis minimum clearance

    flat_idx = inp['mask_indices']   # 0-based Fortran-order linear indices

    # Decode source positions in Fortran (column-major) order
    multi_idx = np.array(np.unravel_index(flat_idx, N, order='F'))  # [ndim, N_src]

    # Per-axis padding required to move source outside PML (+ stencil clearance).
    # A minimum symmetric floor of stencil_clearance is enforced on every axis
    # (independent of source position) so the periodic-wrap stencil seam at each
    # true array edge sits in a dead buffer rather than in the live medium.
    pad = []
    for ax in range(ndim):
        lo = max(0, int(safe_pml[ax]) - int(multi_idx[ax].min()))
        hi = max(0, int(safe_pml[ax]) - (N[ax] - 1 - int(multi_idx[ax].max())))
        if stencil_clearance > 0:
            lo = max(lo, stencil_clearance)
            hi = max(hi, stencil_clearance)
        pad.append((lo, hi))

    if all(lo == 0 and hi == 0 for lo, hi in pad):
        return inp, pad   # nothing to do

    for ax, (lo, hi) in enumerate(pad):
        if lo > 0 or hi > 0:
            # Distinguish source-driven padding from the stencil-clearance floor:
            # if the source already sits >= safe_pml from both edges, the pad is
            # purely the FD stencil-seam buffer, not a source-inside-PML fix.
            src_lo_gap = int(multi_idx[ax].min())
            src_hi_gap = N[ax] - 1 - int(multi_idx[ax].max())
            source_in_pml = (src_lo_gap < safe_pml[ax]) or (src_hi_gap < safe_pml[ax])
            if source_in_pml:
                reason = 'to move source outside PML'
            else:
                reason = f'stencil-clearance buffer (fd half-width {stencil_clearance})'
            print(f'[run_jwave_sim] Padding axis {ax} '
                  f'(src pos {multi_idx[ax].min()}–{multi_idx[ax].max()}, pml={pml_arr[ax]}): '
                  f'({lo}, {hi}) voxels — {reason}')

    # Pad medium arrays by replicating boundary values
    def pad_field(arr):
        if np.isscalar(arr) or (isinstance(arr, np.ndarray) and arr.ndim == 0):
            return arr
        # Medium arrays are in MATLAB axis order (Nx, Ny, Nz)
        return np.pad(arr, pad, mode='edge')

    new_inp = dict(inp)
    for key in ('sound_speed', 'density', 'alpha_coeff', 'alpha_power'):
        new_inp[key] = pad_field(inp[key])

    # New grid dimensions
    new_N = np.array([N[ax] + pad[ax][0] + pad[ax][1] for ax in range(ndim)],
                     dtype=inp['N'].dtype)
    new_inp['N'] = new_N

    # Shift source positions and re-encode as flat Fortran indices
    shifts = np.array([p[0] for p in pad])        # [ndim]
    shifted = multi_idx + shifts[:, np.newaxis]    # [ndim, N_src]
    new_flat = np.ravel_multi_index(shifted, new_N.tolist(), order='F').astype(int)
    new_inp['mask_indices'] = new_flat

    return new_inp, pad


HELMHOLTZ_ALPHA_FLOOR_DB_CM_MHZ2 = 0.002  # physical attenuation of distilled water


def _apply_helmholtz_attenuation_floor(alpha_coeff):
    """Clamp near-zero attenuation up to a small physical floor for the Helmholtz solver.

    alpha_coeff == 0 (e.g. idealised lossless-water benchmarks) makes the Helmholtz
    wavevector term (omega/c)^2 + 2j*omega^3*alpha/c purely real, so (grad^2 + k^2)
    is exactly singular at grid-supported plane-wave modes with |k| = omega/c.
    GMRES then stalls near O(1) residual instead of converging, which is both slow
    (hits gmres_maxiter) and produces a badly under-converged (~18% low) amplitude.
    Flooring alpha_coeff to the physical attenuation of distilled water gives the
    wavevector a small but nonzero imaginary part, removing the exact degeneracy,
    with negligible effect on the field (< 0.006% attenuation over 100 mm at 500 kHz).
    """
    floor = HELMHOLTZ_ALPHA_FLOOR_DB_CM_MHZ2
    if np.isscalar(alpha_coeff):
        floored = max(float(alpha_coeff), floor)
        if floored != alpha_coeff:
            print(f'[run_jwave_sim] Helmholtz: flooring attenuation {float(alpha_coeff):.4g} -> '
                  f'{floor:.4g} dB/(cm*MHz^2) to regularize GMRES (lossless/near-lossless medium)')
        return floored
    arr = np.asarray(alpha_coeff)
    n_floored = int(np.sum(arr < floor))
    if n_floored > 0:
        print(f'[run_jwave_sim] Helmholtz: flooring attenuation in {n_floored}/{arr.size} voxels to '
              f'{floor:.4g} dB/(cm*MHz^2) to regularize GMRES (lossless/near-lossless medium)')
    return np.maximum(arr, floor)


def build_medium_and_domain(inp, use_fourier_series=False):
    """Construct jwave Domain and Medium from input dict.

    use_fourier_series: wrap spatially-varying arrays as FourierSeries fields.
    Required for the Helmholtz solver; plain JAX arrays are used for
    pseudo_spectral to avoid unnecessary FFT wrapping overhead.
    """
    from jwave.geometry import Domain, Medium
    import jax.numpy as jnp

    N  = tuple(inp['N'].tolist())
    dx = tuple(inp['dx'].tolist())
    domain = Domain(N, dx)

    def to_field(val):
        if np.isscalar(val):
            return float(val)
        arr = jnp.array(val, dtype=jnp.float32)
        if use_fourier_series:
            from jaxdf.discretization import FourierSeries
            return FourierSeries(params=arr[..., np.newaxis], domain=domain)
        return arr

    # jwave Medium accepts a scalar pml_size; use the per-axis maximum so every
    # axis is at least as well protected as its per-axis value specifies.
    pml_scalar = int(np.asarray(inp['pml_size']).max())

    # jwave Medium uses alpha in dB/(MHz^y * cm) — same convention as k-Wave.
    # The Helmholtz wavevector operator converts via db2neper with fixed y=2.0.
    medium = Medium(
        domain=domain,
        sound_speed=to_field(inp['sound_speed']),
        density=to_field(inp['density']),
        attenuation=to_field(inp['alpha_coeff']),
        pml_size=pml_scalar,
    )
    return domain, medium


def build_sources(inp, domain):
    """Reconstruct source mask and signals as jwave Sources."""
    from jwave.geometry import Sources
    import jax.numpy as jnp

    N = tuple(inp['N'].tolist())
    flat_indices = inp['mask_indices']
    # MATLAB find() uses column-major (Fortran) linear indexing.
    multi_idx = np.unravel_index(flat_indices, N, order='F')
    positions = tuple(jnp.array(ax) for ax in multi_idx)
    signals   = jnp.array(inp['source_p'], dtype=jnp.float32)  # [N_src x Nt]

    return Sources(positions=positions, signals=signals, dt=inp['dt'], domain=domain)


def build_medium_and_domain_fd(inp, accuracy=8):
    """Construct jwave Domain and Medium with FiniteDifferences field discretization.

    Identical to build_medium_and_domain(..., use_fourier_series=False) except
    medium fields are wrapped as JaxDF FiniteDifferences objects, which causes
    jwave's plum multiple-dispatch to select local FD stencils for all spatial
    derivatives rather than the global FFT-based operators used by pseudo_spectral.

    accuracy: stencil accuracy order (must be even: 2, 4, 6, 8, …).
    """
    from jwave.geometry import Domain, Medium
    from jaxdf.discretization import FiniteDifferences
    import jax.numpy as jnp

    N  = tuple(inp['N'].tolist())
    dx = tuple(inp['dx'].tolist())
    domain = Domain(N, dx)

    def to_fd_field(val):
        if np.isscalar(val):
            return float(val)
        arr = jnp.array(val, dtype=jnp.float32)
        return FiniteDifferences(params=arr[..., np.newaxis], domain=domain, accuracy=accuracy)

    pml_scalar = int(np.asarray(inp['pml_size']).max())
    medium = Medium(
        domain=domain,
        sound_speed=to_fd_field(inp['sound_speed']),
        density=to_fd_field(inp['density']),
        attenuation=to_fd_field(inp['alpha_coeff']),
        pml_size=pml_scalar,
    )
    return domain, medium


# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------
#
# History note: an earlier version of this diagnostic hand-rolled a NumPy
# approximation — ||(bare periodic Laplacian + (ω/c)²)P - source|| / ||source||
# — that omitted three things jwave's own solver applies internally: (1) the
# source scaling `source *= 2/(dx*c_min)` done by
# jwave.acoustics.operators.scale_source_helmholtz, (2) the u = P/(-iω)
# relationship between the field GMRES actually solves for and the returned
# pressure, and (3) the attenuation/PML terms in the real wavevector/Laplacian
# operators. Omitting (1)+(2) alone introduces a missing multiplicative factor
# of |−iω·2/(dx·c_min)| — ~3.35e6 for a 300 kHz / 0.75 mm / c_min=1500 m/s
# benchmark grid — which silently dominated the reported residual. This was
# caught because two runs with very different precision/GMRES settings
# (single/restart=60/maxiter=600 vs double/restart=10/maxiter=1000) produced
# near-identical bogus residuals (~3.6e6) and near-identical p_max, which is
# what you'd expect from a diagnostic bug, not from a GMRES convergence
# problem (a genuine non-convergence is sensitive to precision and restart).
# compute_helmholtz_residual below instead reuses jwave's own `helmholtz`
# operator and `scale_source_helmholtz` directly, so it evaluates exactly the
# linear residual ||A·u - b|| / ||b|| that GMRES targets — no hand-rolled
# approximation, no missing terms.

def compute_helmholtz_residual(P_raw, jw_source, medium, omega, op_params):
    """Exact relative GMRES residual for the Helmholtz solve.

    P_raw, jw_source, medium, op_params are the same objects used to run
    helmholtz_solver (see run_helmholtz). Recovers u = P/(-iω) (the field
    GMRES actually solved for) and evaluates jwave's own `helmholtz` operator
    on it against the correctly-scaled source, so PML, attenuation, and the
    internal source scaling are all handled exactly as jwave itself does —
    unlike a hand-rolled NumPy Laplacian approximation (see history note above).

    Returns the relative residual norm (scalar, dimensionless, python float).
    """
    from jwave.acoustics.operators import helmholtz, scale_source_helmholtz

    u_solved     = P_raw / (-1j * omega)
    scaled_src   = scale_source_helmholtz(jw_source, medium)
    residual_fld = helmholtz(u_solved, medium, omega=omega, params=op_params) - scaled_src

    def _to_array(field):
        return np.array(field.params[..., 0]) if hasattr(field, 'params') else np.array(field)

    residual_arr = _to_array(residual_fld)
    src_arr      = _to_array(scaled_src)
    src_norm     = np.linalg.norm(src_arr)
    if src_norm == 0:
        return float('nan')
    return float(np.linalg.norm(residual_arr) / src_norm)


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

def _absorb_tau_from_power_law(alpha_coeff, alpha_power, freq_hz, sound_speed):
    """Compute Stokes absorb_tau matched to power-law attenuation at freq_hz.

    The Stokes pressure correction p += c0^2 * absorb_tau * rho0 * div(u) has
    the plane-wave dispersion relation (derived from the linearised momentum/
    mass/pressure update, small-absorb_tau expansion):
        k ~= omega/c0 - i * omega^2 * absorb_tau / (2*c0)
    A physically decaying wave (Im(k) > 0 for exp(i(kx-wt))) therefore requires
    absorb_tau < 0 -- absorb_tau must be NEGATIVE, or the term amplifies the
    field instead of damping it (this sign was missing in an earlier version
    and caused persistent late-onset NaN regardless of dt/CFL margin, since
    the instability was structural, not a matter of insufficient headroom).

    We want -omega^2*absorb_tau/(2*c0) to equal the target power-law
    attenuation at freq_hz:
        att_target = alpha_coeff [dB/(cm*MHz^y)] * (freq_hz/1e6)^y * ln(10)/20 * 100

    Solving: absorb_tau = -2 * att_target * c0 / omega_c^2

    This is algebraically identical to k-Wave's fractional-Laplacian
    coefficient absorb_tau = -2*alpha_coeff_Np*c0^(alpha_power-1) (see
    db2neper.m / kspaceFirstOrder_createAbsorptionVariables.m) for the y=2
    case jwave always uses -- k-Wave folds the omega^2 factor into a
    frequency-independent Np-per-(rad/s)^2 coefficient, this form keeps it
    explicit, but the two are the same physical quantity at the design
    frequency.  Off-frequency attenuation will still follow Stokes (f^2)
    rather than the input alpha_power, so the match is exact only at
    freq_hz; for narrow-band CW this is sufficient.

    All array arguments are treated element-wise for spatially varying media.
    """
    f_MHz   = freq_hz / 1e6
    att_Np_m = alpha_coeff * (f_MHz ** alpha_power) * (np.log(10) / 20.0) * 100.0
    omega_c  = 2.0 * np.pi * freq_hz
    return -2.0 * att_Np_m * sound_speed / (omega_c ** 2)


def _check_pstd_stability(inp, dt, absorb_tau):
    """Pre-flight stability check for the PSTD + Stokes absorption scheme.

    Two stability metrics are reported:

    1. Acoustic CFL:  CFL = c_max * dt / dx_min
       Limit for 3-D PSTD: 1 / (pi * sqrt(ndim)) ≈ 0.184

    2. Stokes absorption metric:  |tau_max| * c_max^2 * (pi/dx_min)^2 * dt
       Must be < 1 for the explicit leapfrog absorption step to remain stable
       (this uses |tau_max|, so it is a legitimate magnitude-only bound
       regardless of absorb_tau's sign).  Values above ~0.30 can cause
       late-onset NaN in single precision (slow exponential growth from the
       numerical noise floor).

    Raises RuntimeError on CFL violation or when absorption metric exceeds 0.40.
    Prints a WARNING when the absorption metric is between 0.30 and 0.40.
    The correct fix is to reduce dt (lower grid.source_cfl in the PRESTUS
    config); jwave_wrapper.m performs this adjustment automatically.

    History note: a persistent NaN was previously chased through this metric
    (assuming insufficient margin / dt too coarse) across several rounds of
    tuning, but the actual cause was unrelated to this check: a missing
    negative sign in _absorb_tau_from_power_law made the Stokes term
    energy-injecting instead of dissipative (fixed; see that function's
    docstring). This metric is a real, still-useful CFL-like bound for a
    correctly-signed dissipative term — it was just never the cause of that
    particular bug, and would not have caught it (a sign flip doesn't change
    |tau_max|).
    """
    c_max  = float(np.max(inp['sound_speed']))
    dx_min = float(np.min(inp['dx']))   # smallest spacing → tightest k_nyq
    ndim   = len(inp['N'])

    cfl       = c_max * dt / dx_min
    cfl_limit = 1.0 / (np.pi * np.sqrt(ndim))

    tau_max = float(np.max(np.abs(absorb_tau))) if absorb_tau is not None else 0.0
    k_nyq   = np.pi / dx_min
    absorb_metric = tau_max * c_max ** 2 * k_nyq ** 2 * dt

    # Compute the dt that would bring the metric to 0.25 (safe target)
    dt_safe = 0.25 / (tau_max * c_max ** 2 * k_nyq ** 2) if tau_max > 0 else dt

    print(f'[run_jwave_sim] PSTD stability — acoustic CFL={cfl:.4f} '
          f'(limit {cfl_limit:.4f}), absorption metric={absorb_metric:.4f}')

    if cfl > cfl_limit * 1.01:
        raise RuntimeError(
            f'[run_jwave_sim] PSTD CFL violation: CFL={cfl:.4f} exceeds '
            f'{ndim}-D limit {cfl_limit:.4f}. '
            'Reduce grid.source_cfl or grid.resolution_mm.'
        )

    if absorb_metric > 0.40:
        raise RuntimeError(
            f'[run_jwave_sim] Stokes absorption metric {absorb_metric:.4f} exceeds '
            f'stability limit 0.40 — simulation would produce NaN. '
            f'The correct fix is to reduce dt from {dt:.3e} to ≤ {dt_safe:.3e} s '
            f'(i.e. reduce grid.source_cfl by factor {dt_safe/dt:.3f}). '
            'jwave_wrapper.m should have applied this automatically; '
            'check that the auto-dt-reduction block ran before HDF5 was written.'
        )

    if absorb_metric > 0.30:
        print(
            f'[run_jwave_sim] WARNING: Stokes absorption metric {absorb_metric:.4f} > 0.30. '
            'High skull attenuation is tightening the PSTD stability margin. '
            f'If NaN occurs, reduce dt to ≤ {dt_safe:.3e} s (lower grid.source_cfl '
            'by factor {:.3f}).'.format(dt_safe / dt)
        )


def _report_step_growth(step_absmax_all, dt, tag):
    """Print per-step field-magnitude diagnostics for a time-domain scan.

    p_max_all/p_min_all are only updated once n >= record_start, so a NaN
    produced early in the run (e.g. step 50 of 7447) is invisible until the
    recording window is reached, and the final "produced NaN output" error
    gives no information about when instability actually began.  This tracks
    max(|p|) at EVERY step, independent of record_start, so the true onset
    step/time can be read directly off the log instead of re-running with
    ad-hoc instrumentation after the fact.
    """
    n_steps = len(step_absmax_all)
    finite  = np.isfinite(step_absmax_all)
    if not np.all(finite):
        first_bad = int(np.argmax(~finite))
        print(f'[run_jwave_sim] [{tag}] WARNING: field magnitude goes non-finite '
              f'at step {first_bad}/{n_steps} (t={first_bad*dt*1e6:.2f} us) — '
              'this is BEFORE the recording window in general; the actual '
              'instability starts here, not at the step where the final NaN '
              'check fires.')
    n_samples   = min(10, n_steps)
    sample_idx  = np.linspace(0, n_steps - 1, n_samples).astype(int)
    trend       = ', '.join(f'{step_absmax_all[i]:.3e}' for i in sample_idx)
    print(f'[run_jwave_sim] [{tag}] max|p| trend over the run (sampled): {trend}')


def run_pseudo_spectral(inp):
    """
    Run jwave pseudo-spectral time-domain (PSTD) simulation.

    Implements the same staggered-grid symplectic integrator as kspaceFirstOrder*
    with the Stokes (y=2) absorption correction applied at the pressure update:

        p = c₀² · (ρ + absorb_tau · ρ₀ · ∇·u)

    absorb_tau is derived via _absorb_tau_from_power_law so that the effective
    Stokes attenuation matches the input power-law (alpha_coeff, alpha_power) at
    freq_hz.  Spatially varying alpha_coeff is supported.

    Returns (p_max_all, p_min_all, p_final) as numpy arrays.
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import TimeAxis
    from jwave.acoustics.time_varying import (
        TimeWavePropagationSettings,
        ongrid_wave_prop_params,
        momentum_conservation_rhs,
        mass_conservation_rhs,
        diag_jacobian,
        sum_over_dims,
        pressure_from_density,
    )
    # FourierSeries fields are required for jwave's spectral operator dispatch to
    # select optimised FFT-based stencils.  Plain JAX arrays trigger a slow generic
    # fallback path (jwave warns: "All inputs are Arrays. This is not recommended").
    domain, medium = build_medium_and_domain(inp, use_fourier_series=True)
    sources = build_sources(inp, domain)

    Nt   = inp['Nt']
    dt   = inp['dt']
    time_axis    = TimeAxis(dt=dt, t_end=Nt * dt)
    record_start = inp['record_start_index']  # 0-based

    print(f'[run_jwave_sim] Pseudo-spectral solver: N={inp["N"]}, Nt={Nt}, dt={dt:.3e} s')
    print(f'[run_jwave_sim] Recording from step {record_start} to {Nt}')

    settings    = TimeWavePropagationSettings()
    wave_params = ongrid_wave_prop_params(medium, time_axis, settings=settings)
    c_ref   = wave_params['c_ref']
    pml_rho = wave_params['pml_rho']
    pml_u   = wave_params['pml_u']

    N         = tuple(inp['N'].tolist())
    shape     = tuple(list(N) + [len(N)])
    shape_one = tuple(list(N) + [1])

    u0  = pml_u.replace_params(jnp.zeros(shape))
    p0  = pml_rho.replace_params(jnp.zeros(shape_one))
    rho0 = p0.replace_params(jnp.zeros(shape))

    # Precompute Stokes absorb_tau matched to power-law attenuation at freq_hz.
    # The Stokes term gives att(omega)=absorb_tau*omega^2/(2*c0); we rescale so
    # that this equals the target power-law attenuation at the centre frequency.
    raw_tau = _absorb_tau_from_power_law(
        inp['alpha_coeff'], inp['alpha_power'], inp['freq_hz'], inp['sound_speed'])
    _check_pstd_stability(inp, dt, raw_tau)
    if np.isscalar(raw_tau) or (isinstance(raw_tau, np.ndarray) and raw_tau.size == 1):
        absorb_tau = jnp.float32(float(raw_tau))
    else:
        absorb_tau = jnp.array(raw_tau.astype(np.float64),
                               dtype=jnp.float32)[..., jnp.newaxis]

    rho0_arr = medium.density   # plain JAX array (not a Field here)
    if hasattr(rho0_arr, 'ndim') and rho0_arr.ndim == len(N):
        rho0_arr = rho0_arr[..., jnp.newaxis]

    neg_inf = jnp.full(N, -jnp.inf, dtype=jnp.float32)
    pos_inf = jnp.full(N,  jnp.inf, dtype=jnp.float32)

    import time as _time

    report_interval = max(1, Nt // 10)
    _eta_state = {}

    def _progress_cb(n):
        if n % report_interval == 0:
            print(f'[run_jwave_sim]   step {int(n):>{len(str(Nt))}}/{Nt}', flush=True)
            # k-Wave-style ETA: time the interval between the 1st and 2nd
            # checkpoint (both post-compile) and linearly extrapolate to the
            # remaining steps. The 1st checkpoint alone is skipped because its
            # elapsed time is contaminated by XLA compile time.
            if n == report_interval:
                _eta_state['t0'] = _time.time()
            elif n == 2 * report_interval and 't0' in _eta_state:
                interval_s = _time.time() - _eta_state['t0']
                remaining_steps = Nt - n
                eta_s = interval_s / report_interval * remaining_steps
                print(f'[run_jwave_sim]   estimated remaining time: {eta_s:.1f} s', flush=True)

    def step(carry, n):
        p, u, rho, p_max, p_min = carry

        jax.debug.callback(_progress_cb, n)

        # Momentum update: u ← pml·(pml·u − dt·(1/ρ₀)·∇p)
        du  = momentum_conservation_rhs(p, u, medium, c_ref=c_ref, dt=dt)
        u   = pml_u * (pml_u * u + dt * du)

        # Mass update: ρ ← pml·(pml·ρ − dt·ρ₀·∇·u) + source
        # diag_jacobian returns the per-axis velocity gradients; recompute here
        # to obtain ∇·u for the Stokes pressure correction.
        du_tensor = diag_jacobian(u, stagger=[-0.5])
        div_u     = sum_over_dims(du_tensor)           # scalar divergence field

        drho = mass_conservation_rhs(
            p, u, sources.on_grid(n), medium, c_ref=c_ref, dt=dt)
        rho  = pml_rho * (pml_rho * rho + dt * drho)

        # Pressure: lossless term + Stokes absorption correction
        c0_sq = medium.sound_speed ** 2
        p = pressure_from_density(rho, medium)         # c₀²·Σρ
        p = p + c0_sq * absorb_tau * rho0_arr * div_u  # Stokes: absorb_tau < 0 (dissipative)

        p_arr     = p.params[..., 0]
        recording = n >= record_start
        p_max = jnp.where(recording, jnp.maximum(p_max, p_arr), p_max)
        p_min = jnp.where(recording, jnp.minimum(p_min, p_arr), p_min)

        step_absmax = jnp.max(jnp.abs(p_arr))
        return (p, u, rho, p_max, p_min), step_absmax

    init_carry = (p0, u0, rho0, neg_inf, pos_inf)

    @jax.jit
    def _scan(carry, xs):
        return jax.lax.scan(step, carry, xs)

    print('[run_jwave_sim] Compiling PSTD scan (first run may be slow on GPU)...')
    t0 = _time.time()
    result = _scan(init_carry, jnp.arange(Nt))
    # block_until_ready ensures the timer covers actual execution, not just dispatch
    jax.block_until_ready(result)
    print(f'[run_jwave_sim] PSTD scan done in {_time.time() - t0:.1f} s')
    (p_final_field, _, _, p_max_all, p_min_all), step_absmax_all = result
    _report_step_growth(np.array(step_absmax_all), dt, 'pseudo_spectral')

    p_final = p_final_field.params[..., 0]
    return np.array(p_max_all), np.array(p_min_all), np.array(p_final)


def run_finite_differences(inp, accuracy=8):
    """
    Run jwave finite-difference time-domain (FDTD) simulation with Stokes absorption.

    Same update equations and Stokes (y=2) absorption correction as
    run_pseudo_spectral, but medium fields are wrapped as JaxDF FiniteDifferences
    so that jwave's operator dispatch selects local FD stencils (order `accuracy`)
    rather than global FFT derivatives.  Useful for benchmarking numerical
    dispersion characteristics vs pseudo_spectral.

    Returns (p_max_all, p_min_all, p_final) as numpy arrays.
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import TimeAxis
    from jwave.acoustics.time_varying import (
        TimeWavePropagationSettings,
        ongrid_wave_prop_params,
        momentum_conservation_rhs,
        mass_conservation_rhs,
        diag_jacobian,
        sum_over_dims,
        pressure_from_density,
    )
    domain, medium = build_medium_and_domain_fd(inp, accuracy=accuracy)
    sources = build_sources(inp, domain)

    Nt   = inp['Nt']
    dt   = inp['dt']
    time_axis    = TimeAxis(dt=dt, t_end=Nt * dt)
    record_start = inp['record_start_index']

    print(f'[run_jwave_sim] FD time-domain solver: N={inp["N"]}, Nt={Nt}, '
          f'dt={dt:.3e} s, FD accuracy order={accuracy}')
    print(f'[run_jwave_sim] Recording from step {record_start} to {Nt}')

    settings    = TimeWavePropagationSettings()
    wave_params = ongrid_wave_prop_params(medium, time_axis, settings=settings)
    c_ref   = wave_params['c_ref']
    pml_rho = wave_params['pml_rho']
    pml_u   = wave_params['pml_u']

    N         = tuple(inp['N'].tolist())
    shape     = tuple(list(N) + [len(N)])
    shape_one = tuple(list(N) + [1])

    u0   = pml_u.replace_params(jnp.zeros(shape))
    p0   = pml_rho.replace_params(jnp.zeros(shape_one))
    rho0 = p0.replace_params(jnp.zeros(shape))

    raw_tau = _absorb_tau_from_power_law(
        inp['alpha_coeff'], inp['alpha_power'], inp['freq_hz'], inp['sound_speed'])
    _check_pstd_stability(inp, dt, raw_tau)
    if np.isscalar(raw_tau) or (isinstance(raw_tau, np.ndarray) and raw_tau.size == 1):
        absorb_tau = jnp.float32(float(raw_tau))
    else:
        absorb_tau = jnp.array(raw_tau.astype(np.float64),
                               dtype=jnp.float32)[..., jnp.newaxis]

    rho0_arr = medium.density
    if hasattr(rho0_arr, 'ndim') and rho0_arr.ndim == len(N):
        rho0_arr = rho0_arr[..., jnp.newaxis]

    neg_inf = jnp.full(N, -jnp.inf, dtype=jnp.float32)
    pos_inf = jnp.full(N,  jnp.inf, dtype=jnp.float32)

    import time as _time

    report_interval = max(1, Nt // 10)
    _eta_state = {}

    def _progress_cb(n):
        if n % report_interval == 0:
            print(f'[run_jwave_sim]   step {int(n):>{len(str(Nt))}}/{Nt}', flush=True)
            # k-Wave-style ETA: time the interval between the 1st and 2nd
            # checkpoint (both post-compile) and linearly extrapolate to the
            # remaining steps. The 1st checkpoint alone is skipped because its
            # elapsed time is contaminated by XLA compile time.
            if n == report_interval:
                _eta_state['t0'] = _time.time()
            elif n == 2 * report_interval and 't0' in _eta_state:
                interval_s = _time.time() - _eta_state['t0']
                remaining_steps = Nt - n
                eta_s = interval_s / report_interval * remaining_steps
                print(f'[run_jwave_sim]   estimated remaining time: {eta_s:.1f} s', flush=True)

    def step(carry, n):
        p, u, rho, p_max, p_min = carry

        jax.debug.callback(_progress_cb, n)

        du        = momentum_conservation_rhs(p, u, medium, c_ref=c_ref, dt=dt)
        u         = pml_u * (pml_u * u + dt * du)

        du_tensor = diag_jacobian(u, stagger=[-0.5])
        div_u     = sum_over_dims(du_tensor)

        drho = mass_conservation_rhs(
            p, u, sources.on_grid(n), medium, c_ref=c_ref, dt=dt)
        rho  = pml_rho * (pml_rho * rho + dt * drho)

        c0_sq = medium.sound_speed ** 2
        p = pressure_from_density(rho, medium)
        p = p + c0_sq * absorb_tau * rho0_arr * div_u

        p_arr     = p.params[..., 0]
        recording = n >= record_start
        p_max = jnp.where(recording, jnp.maximum(p_max, p_arr), p_max)
        p_min = jnp.where(recording, jnp.minimum(p_min, p_arr), p_min)

        step_absmax = jnp.max(jnp.abs(p_arr))
        return (p, u, rho, p_max, p_min), step_absmax

    init_carry = (p0, u0, rho0, neg_inf, pos_inf)

    @jax.jit
    def _scan(carry, xs):
        return jax.lax.scan(step, carry, xs)

    import time as _time
    print('[run_jwave_sim] Compiling FD scan (first run may be slow on GPU)...')
    t0 = _time.time()
    result = _scan(init_carry, jnp.arange(Nt))
    jax.block_until_ready(result)
    print(f'[run_jwave_sim] FD scan done in {_time.time() - t0:.1f} s')
    (p_final_field, _, _, p_max_all, p_min_all), step_absmax_all = result
    _report_step_growth(np.array(step_absmax_all), dt, 'finite_differences')

    p_final = p_final_field.params[..., 0]
    return np.array(p_max_all), np.array(p_min_all), np.array(p_final)


def run_helmholtz(inp, gmres_restart=10, gmres_maxiter=1000, gmres_tol=1e-3, debug=False):
    """
    Run jwave frequency-domain (Helmholtz) solver.

    Solves the heterogeneous Helmholtz equation at freq_hz using GMRES with
    Fourier spectral collocation and PML absorbing boundaries. Attenuation is
    applied via db2neper(alpha_coeff, y=2.0) inside jwave's wavevector operator.

    Lossless media (alpha_coeff = 0): the Helmholtz operator (∇² + k²) has a
    purely real wavenumber, making it near-singular in the bulk of the domain.
    This is regularized via _apply_helmholtz_attenuation_floor rather than by
    tightening the GMRES budget. Default restart/maxiter/tol match jwave's own
    helmholtz_solver defaults (10/1000/1e-3): empirically, tightening these
    (e.g. restart=40-60, maxiter=300-600, tol=1e-5, as earlier PRESTUS defaults
    and some benchmark scripts did) produced an identical converged answer but
    took longer once the attenuation floor properly regularizes the operator —
    the tighter budget was compensating for the missing floor, not adding
    accuracy. Medium parameters (alpha_coeff, c0, rho0) are used as passed from
    the upstream MATLAB parameters struct without modification (besides the
    attenuation floor).

    Returns:
        p_max_all = |P|        (peak amplitude)
        p_min_all = -|P|       (negative amplitude, for MI compatibility)
        p_final   = Re(P)      (real part)
        P_complex              (full complex pressure field)
        source_field           (complex source array in Pa/m², for provenance)
        residual_norm          (exact relative GMRES residual, see
                                compute_helmholtz_residual)
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import Domain, Medium
    from jwave.acoustics.operators import helmholtz
    from jwave.acoustics.time_harmonic import helmholtz_solver

    inp = dict(inp)
    inp['alpha_coeff'] = _apply_helmholtz_attenuation_floor(inp['alpha_coeff'])
    domain, medium = build_medium_and_domain(inp, use_fourier_series=True)

    freq_hz = inp['freq_hz']
    omega   = 2 * np.pi * freq_hz

    print(f'[run_jwave_sim] Helmholtz solver: N={inp["N"]}, freq={freq_hz/1e3:.1f} kHz, '
          f'GMRES restart={gmres_restart}, maxiter={gmres_maxiter}')

    # Extract complex source amplitude per element via DFT at freq_hz, restricted
    # to the steady-state window defined by record_start_index — the same "last
    # N periods, after the wavefront has crossed the domain" window that k-Wave
    # and jwave's own time-domain solvers (run_pseudo_spectral/run_finite_differences)
    # use for p_max_all. createCWSignals() applies a startup cosine up-ramp
    # (default 4 periods) to source_p; DFT-ing the full Nt-length signal averages
    # that lower-amplitude ramp in with the steady CW cycles and systematically
    # underestimates the drive amplitude (the time-domain solvers don't have this
    # problem because they peak-track with max/min rather than average). Slicing
    # to record_start_index also removes most DFT scalloping loss "for free":
    # time_points_to_record in source_sensor_setup.m is built from an integer
    # number of samples per period (round(num_periods * wave_period / dt)), so
    # the steady-state window is (to within one time step) an exact integer
    # number of cycles of freq_hz.
    Nt = inp['Nt']
    dt = inp['dt']
    record_start = inp['record_start_index']
    source_p = inp['source_p'][:, record_start:]
    Nt_steady = source_p.shape[1]
    if record_start > 0:
        print(f'[run_jwave_sim] Helmholtz: DFT source amplitude over steady-state window '
              f'[{record_start}:{Nt}] ({Nt_steady}/{Nt} samples) — excludes startup ramp.')
    freqs     = np.fft.rfftfreq(Nt_steady, dt)
    freq_idx  = int(np.argmin(np.abs(freqs - freq_hz)))
    freq_bin  = freqs[freq_idx]
    freq_err_pct = abs(freq_bin - freq_hz) / freq_hz * 100
    if freq_err_pct > 1.0:
        print(f'[run_jwave_sim] WARNING: requested freq {freq_hz/1e3:.3f} kHz, '
              f'nearest FFT bin is {freq_bin/1e3:.4f} kHz '
              f'(error={freq_err_pct:.2f}%). '
              'Increase Nt or adjust freq_hz for more accurate frequency matching.',
              file=sys.stderr)
    # np.fft.rfft uses the e^{-i2*pi*k*n/N} convention, so for a real signal
    # x(t) = A*cos(omega*t + phi) it yields A*e^{+i*phi}. jwave's own Green's
    # function (angular_spectrum's exp(1j*k*r)/r, an outgoing wave) and its
    # lossy wavevector (+2j*omega^3*alpha/c, decaying under e^{-i*omega*t})
    # both fix jwave's time convention as p(x,t) = Re[P(x)*e^{-i*omega*t}], for
    # which the correct phasor is P = A*e^{-i*phi} -- the CONJUGATE of the raw
    # DFT coefficient. Feeding the unconjugated value flips the sign of every
    # element's phase relative to jwave's convention: nearly harmless for a
    # near-uniform phase array (only relative phase differences steer the
    # beam), but a real steering phase gradient gets its focusing correction
    # applied backwards, displacing the focus.
    fft_coeff = np.conj(np.fft.rfft(source_p, axis=-1)[:, freq_idx]) * (2 / Nt_steady)  # [N_src] complex

    N = tuple(inp['N'].tolist())
    source_field = np.zeros(N, dtype=complex)
    flat_idx     = inp['mask_indices']
    multi_idx    = np.unravel_index(flat_idx, N, order='F')

    # jwave's helmholtz_solver internally scales the source by 2/(dx*c_min)
    # and solves (∇²+k²)φ = scaled_source, returning p = -iω·φ.
    # Working backwards, source_in must be in Pa (pressure amplitude).
    # Multiplying by k² here would inflate the result by ~k² ≈ 10⁶.
    source_field[multi_idx] = fft_coeff

    from jaxdf.discretization import FourierSeries
    jw_source = FourierSeries(
        params=jnp.array(source_field[..., np.newaxis], dtype=jnp.complex64),
        domain=domain,
    )

    # Precompute the Laplacian/PML/wavevector operator params outside jax.jit.
    # jwave's own harmonic-solver notebook (docs/notebooks/harmonic/helmholtz_problem.ipynb)
    # notes that in 3D this reduces compilation time (at a slightly higher per-call
    # runtime) versus leaving XLA to constant-fold the PML/Laplacian setup inside the
    # jitted GMRES trace — see https://github.com/google/jax/issues?q=constant+folding.
    # Grid size scales directly with grid.resolution_mm, so this matters most for the
    # fine (e.g. 0.5 mm) resolutions used by the free-water benchmark.
    #
    # NOTE: checkpoint must stay at jwave's default (True). checkpoint=False was
    # tried here and reproducibly crashes GMRES with a JAX UnexpectedTracerError
    # ("leaked intermediate value" from laplacian_with_pml's cached fft_rho0 params)
    # for ANY array-valued medium field (heterogeneous OR constant-valued arrays) —
    # i.e. every real PRESTUS run, since medium_setup.m always writes full grid
    # arrays, never bare scalars. Verified empirically before landing this change.
    op_params = helmholtz.default_params(jw_source, medium, omega=omega)

    @jax.jit
    def forward(medium, source, op_params):
        return helmholtz_solver(medium, omega, source, params=op_params,
                                restart=gmres_restart, maxiter=gmres_maxiter,
                                tol=gmres_tol)

    if debug:
        # jax.scipy.sparse.linalg.gmres runs a fixed restart*maxiter iteration
        # budget with no data-dependent early exit, so a second call on the same
        # (already-compiled) inputs is a pure re-execution — the gap between the
        # first call (compile+solve) and the second (solve only) isolates how
        # much of the wall-clock time is XLA compilation vs actual GMRES work.
        import time as _time
        t0 = _time.time()
        P_raw = jax.block_until_ready(forward(medium, jw_source, op_params))
        t_first = _time.time() - t0
        t0 = _time.time()
        P_raw = jax.block_until_ready(forward(medium, jw_source, op_params))
        t_second = _time.time() - t0
        print(f'[run_jwave_sim] [debug] compile+solve={t_first:.2f}s  '
              f'solve-only (2nd call, cache warm)={t_second:.2f}s  '
              f'estimated compile time={max(0.0, t_first - t_second):.2f}s')
    else:
        P_raw = forward(medium, jw_source, op_params)
    residual_norm = compute_helmholtz_residual(P_raw, jw_source, medium, omega, op_params)

    if hasattr(P_raw, 'params'):
        P_complex = np.array(P_raw.params[..., 0])
    else:
        P_complex = np.array(P_raw)

    p_amp     = np.abs(P_complex)
    p_max_all = p_amp
    p_min_all = -p_amp
    p_final   = np.real(P_complex)

    return p_max_all, p_min_all, p_final, P_complex, source_field, residual_norm


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _print_jax_device_info():
    """Print which JAX backend and device(s) are active."""
    try:
        import jax
        backend = jax.default_backend()
        devices = jax.devices()
        dev_strs = [str(d) for d in devices]
        print(f'[run_jwave_sim] JAX backend : {backend}')
        print(f'[run_jwave_sim] JAX devices : {", ".join(dev_strs)}')
    except Exception as e:
        print(f'[run_jwave_sim] JAX device info unavailable: {e}')


def _print_medium_stats(inp):
    """Print min/max of key medium fields."""
    for key in ('sound_speed', 'density', 'alpha_coeff'):
        val = inp[key]
        if np.isscalar(val):
            print(f'[run_jwave_sim]   {key:<14}: {val:.4g} (scalar)')
        else:
            print(f'[run_jwave_sim]   {key:<14}: [{np.min(val):.4g}, {np.max(val):.4g}]')


def _debug_print_domain_pml_info(inp):
    """Print domain/PML sizing in physical + wavelength terms (--debug only).

    jwave's Medium.pml_size is a single scalar (unlike k-Wave, which supports
    a per-axis PML) — see jwave.geometry.Medium. run_jwave_sim.py's max() over
    the per-axis grid.pml_size vector (itself usually 'auto'-derived from
    k-Wave's getOptimalPMLSize, an FFT-prime-factorization heuristic tuned for
    k-Wave's own backend, not for GMRES cost or physical absorption) is the
    only way to interface with that scalar API. This print exists so a
    disproportionate PML fraction (small free-water domains are most at risk,
    since 'auto' pulls from a fixed absolute point range regardless of how
    large the useful interior actually is) is visible rather than silently
    inflating every GMRES iteration's cost.
    """
    N = np.asarray(inp['N']).tolist()
    dx = np.asarray(inp['dx']).tolist()
    c_min = float(np.min(inp['sound_speed']))
    freq_hz = inp['freq_hz']
    wavelength_m = c_min / freq_hz
    pml_arr = np.asarray(inp['pml_size']).flatten()
    pml_scalar = int(pml_arr.max())

    ppw_per_axis = [wavelength_m / d for d in dx]
    pml_wavelengths = pml_scalar / np.mean(ppw_per_axis)

    total_voxels = int(np.prod(N))
    interior = [max(1, n - 2 * pml_scalar) for n in N]
    interior_voxels = int(np.prod(interior))
    pml_fraction = 1.0 - interior_voxels / total_voxels

    print(f'[run_jwave_sim] [debug] wavelength={wavelength_m*1e3:.3f} mm  '
          f'PPW={[round(p, 1) for p in ppw_per_axis]}')
    print(f'[run_jwave_sim] [debug] PML per-axis={pml_arr.tolist()}  '
          f'-> scalar used by jwave={pml_scalar} ({pml_wavelengths:.2f} wavelengths)')
    print(f'[run_jwave_sim] [debug] domain N={N}  interior={interior}  '
          f'PML fraction of total voxels={pml_fraction*100:.1f}%')
    if pml_fraction > 0.25:
        print(f'[run_jwave_sim] [debug] WARNING: PML consumes {pml_fraction*100:.0f}% of the '
              'domain — grid.pml_size=\'auto\' (k-Wave\'s FFT-prime-factorization heuristic, '
              'not tuned for GMRES cost) may be far larger than the '
              f'{pml_wavelengths:.1f}-wavelength physical requirement for this small a domain. '
              'Consider an explicit smaller grid.pml_size for this benchmark.')


def _debug_print_steady_state_margin(inp):
    """Print how much margin the recording window has past the transient wavefront
    arrival, for the time-domain (PSTD/FD) solvers (--debug only).

    source_sensor_setup.m sets t_end = domain_diagonal / c_min (a single one-way
    domain transit) and records only the last few cycles. jwave's own k-Wave
    cross-validation test (tests/test_kwave_helmholtz.m in the jwave repo) uses a
    much larger margin — roundtrips=20 full domain-diagonal transits — before its
    equivalent recording window, specifically to let any residual PML reflection
    decay before treating the field as "steady state". A tight margin here means
    parts of the domain (particularly far corners of a large free-water grid)
    may not have reached genuine periodic steady state by the time p_max_all is
    recorded, which would show up as a real (non-GMRES-related) discrepancy
    against the Helmholtz solver's exact steady-state solution.
    """
    N = np.asarray(inp['N']).tolist()
    dx = np.asarray(inp['dx']).tolist()
    c_min = float(np.min(inp['sound_speed']))
    dt = inp['dt']
    Nt = inp['Nt']
    record_start = inp['record_start_index']

    diag_m = float(np.sqrt(sum((n * d) ** 2 for n, d in zip(N, dx))))
    one_way_s = diag_m / c_min
    total_s = Nt * dt
    record_start_s = record_start * dt
    margin_s = record_start_s - one_way_s
    period_s = 1.0 / inp['freq_hz']

    print(f'[run_jwave_sim] [debug] domain diagonal={diag_m*1e3:.1f} mm  '
          f'one-way transit={one_way_s*1e6:.2f} us  total sim duration={total_s*1e6:.2f} us')
    print(f'[run_jwave_sim] [debug] record_start={record_start}/{Nt} '
          f'({record_start_s*1e6:.2f} us)  margin past one-way transit={margin_s*1e6:.2f} us '
          f'({margin_s/period_s:.2f} periods)')
    if margin_s < 3 * period_s:
        print(f'[run_jwave_sim] [debug] WARNING: recording starts only '
              f'{margin_s/period_s:.2f} periods after the wavefront reaches the far domain '
              'corner (vs. jwave\'s own roundtrips=20 validation margin). Points near the '
              'domain boundary may not be in true periodic steady state yet — if the '
              'recorded field disagrees with Helmholtz\'s exact steady-state solution, '
              'try extending grid.acoustic_fov or increasing sensor.record_start_index margin.')


def main():
    args = parse_args()

    try:
        _print_jax_device_info()
        inp = read_input(args.input)
        _print_medium_stats(inp)
        if args.debug:
            _debug_print_domain_pml_info(inp)

        diagnostics = {}

        # Pad the grid when source elements fall inside jwave's PML zone.
        # k-Wave adds PML *outside* the specified grid; jwave applies it
        # *within* the grid.  Without padding, a transducer at z=2 (standard
        # PRESTUS placement) is inside the PML and its energy is absorbed.
        orig_N = tuple(inp['N'].tolist())
        fd_clearance = args.fd_accuracy // 2 if args.solver == 'finite_differences' else 0
        inp, _pad = auto_pad_for_pml(inp, stencil_clearance=fd_clearance)

        if args.debug and args.solver in ('pseudo_spectral', 'finite_differences'):
            _debug_print_steady_state_margin(inp)

        if args.solver == 'pseudo_spectral':
            p_max_all, p_min_all, p_final = run_pseudo_spectral(inp)
            if np.any(np.isnan(np.asarray(p_max_all))):
                raise RuntimeError(
                    '[run_jwave_sim] pseudo_spectral produced NaN output. '
                    'Check the [pseudo_spectral] max|p| trend printed above to see '
                    'when growth started (a NaN from the first few steps points to a '
                    'structural bug, e.g. a sign error in the Stokes term, rather than '
                    'a tuning issue). Otherwise likely causes: (1) Stokes absorption '
                    'instability — reduce grid.source_cfl so that dt decreases '
                    '(jwave_wrapper.m should apply this automatically); '
                    '(2) acoustic CFL violation — reduce dt or grid resolution.'
                )

        elif args.solver == 'finite_differences':
            p_max_all, p_min_all, p_final = run_finite_differences(inp, accuracy=args.fd_accuracy)
            if np.any(np.isnan(np.asarray(p_max_all))):
                raise RuntimeError(
                    '[run_jwave_sim] finite_differences produced NaN output. '
                    'Check the [finite_differences] max|p| trend printed above to see '
                    'when growth started (a NaN from the first few steps points to a '
                    'structural bug, e.g. a sign error in the Stokes term, rather than '
                    'a tuning issue). Otherwise likely causes: (1) Stokes absorption '
                    'instability — reduce grid.source_cfl (jwave_wrapper.m should apply '
                    'this automatically); (2) stencil reaching into PML — try a smaller '
                    f'fd_accuracy (currently {args.fd_accuracy}) or larger pml_size; '
                    '(3) acoustic CFL violation — reduce dt or grid resolution.'
                )

        else:  # helmholtz
            p_max_all, p_min_all, p_final, P_complex, source_field, residual_norm = run_helmholtz(
                inp,
                gmres_restart=args.gmres_restart,
                gmres_maxiter=args.gmres_maxiter,
                gmres_tol=args.gmres_tol,
                debug=args.debug,
            )
            print(f'[run_jwave_sim] Helmholtz spectral residual: {residual_norm:.4e}')
            diagnostics['residual_norm'] = residual_norm
            if np.any(np.isnan(np.asarray(p_max_all))):
                raise RuntimeError(
                    '[run_jwave_sim] helmholtz produced NaN output. '
                    f'GMRES spectral residual: {residual_norm:.4e}. '
                    'Likely causes: GMRES non-convergence — increase gmres_maxiter or gmres_restart.'
                )

        # Crop padded output back to the original grid size
        if any(lo > 0 or hi > 0 for lo, hi in _pad):
            slices = tuple(
                slice(lo, lo + orig_N[ax])
                for ax, (lo, _hi) in enumerate(_pad)
            )
            p_max_all = p_max_all[slices]
            p_min_all = p_min_all[slices]
            p_final   = p_final[slices]
            print(f'[run_jwave_sim] Cropped padded output back to original grid {orig_N}')

        # Validate output shapes
        expected_shape = orig_N
        for arr, name in [(p_max_all, 'p_max_all'), (p_min_all, 'p_min_all'), (p_final, 'p_final')]:
            if arr.shape != expected_shape:
                raise ValueError(
                    f'Output {name} has shape {arr.shape} but expected grid N={expected_shape}. '
                    'Check solver output and dimension transpositions in run_jwave_sim.py.'
                )
        print(f'[run_jwave_sim] p_max={p_max_all.max():.1f} Pa  p_min={p_min_all.min():.1f} Pa')
        write_output(args.output, p_max_all, p_min_all, p_final, diagnostics=diagnostics)

    except Exception as e:
        import traceback
        print(f'[run_jwave_sim] ERROR: {e}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
