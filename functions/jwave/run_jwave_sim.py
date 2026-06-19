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
    pseudo_spectral   - Pseudo-spectral time-domain (PSTD) propagation via
                        jwave.acoustics.time_varying.simulate_wave_propagation.
                        Equivalent to kspaceFirstOrder*. Uses a semi-implicit first-order
                        Euler integrator with pseudo-spectral (global FFT) spatial derivatives.
                        Returns p_max_all, p_min_all, p_final over the recorded window.
                        NOTE: jwave ignores alpha_power (power-law exponent); attenuation
                        is applied with a fixed internal exponent regardless of the value
                        written to the HDF5. Results will differ from k-Wave when alpha_power
                        deviates from the jwave default.
    helmholtz         - Frequency-domain GMRES solver via
                        jwave.acoustics.time_harmonic.helmholtz_solver.
                        More efficient for CW simulations. Returns p_max_all
                        (amplitude), p_min_all (negative amplitude), p_final
                        (real part of complex pressure). Requires freq_hz in
                        the HDF5 metadata attribute.
    finite_differences - Finite-difference time-domain (FDTD) propagation via
                        jwave.acoustics.time_varying.simulate_wave_propagation
                        with medium fields wrapped as JaxDF FiniteDifferences
                        instead of FourierSeries. Uses the same time-stepping
                        loop as pseudo_spectral but with spatially-local FD
                        stencils rather than global FFT derivatives. Accuracy
                        order controlled by --fd_accuracy (default 8).

Exit codes:
    0  — success
    1  — error (message printed to stderr)

Diagnostics written to /diagnostics/ in the output HDF5 (helmholtz only):
    residual_norm             — relative spectral residual ||ΔP||/||b|| (scalar)
    benchmark_correlation     — spatial correlation vs time-domain on FOV subgrid (if --benchmark_fov_vox)
    benchmark_amplitude_ratio — helmholtz/pseudo_spectral peak amplitude ratio on FOV subgrid
    benchmark_peak_dist_vox   — distance between peak locations on FOV subgrid (voxels)
"""

import argparse
import os
import sys
import numpy as np
import h5py
import jax
jax.config.update("jax_compilation_cache_dir", os.path.join(os.path.expanduser("~"), ".cache", "jax_xla"))


def parse_args():
    p = argparse.ArgumentParser(description='PRESTUS–jwave bridge')
    p.add_argument('--input',  required=True, help='Path to input HDF5 from MATLAB')
    p.add_argument('--output', required=True, help='Path to write output HDF5')
    p.add_argument('--solver', default='pseudo_spectral',
                   choices=['pseudo_spectral', 'helmholtz', 'finite_differences'],
                   help='Solver: pseudo_spectral (PSTD, equivalent to kspaceFirstOrder*), '
                        'helmholtz (freq-domain GMRES), '
                        'finite_differences (FDTD with local FD stencils; default: pseudo_spectral)')
    p.add_argument('--fd_accuracy', type=int, default=8,
                   help='FD stencil accuracy order for finite_differences solver (default: 8)')
    p.add_argument('--gmres_restart', type=int, default=40,
                   help='GMRES restart parameter (default: 40; JAX default is 20)')
    p.add_argument('--gmres_maxiter', type=int, default=300,
                   help='Max GMRES outer iterations (default: 300)')
    p.add_argument('--benchmark_fov_vox', type=str, default='',
                   help='Crop box for FOV benchmark: "x0,y0,z0,x1,y1,z1" (0-based inclusive). '
                        'Helmholtz only: runs both solvers on this subgrid and compares results.')
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

        pml_size = int(np.asarray(f['/pml/size'][()]).flat[0])

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
    print(f'[run_jwave_sim] Input: N={N.tolist()}  Nt={Nt}  dt={dt:.3e} s  '
          f'N_src={n_src_mask}  freq={freq_hz/1e3:.1f} kHz')

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

def build_medium_and_domain(inp, use_fourier_series=False):
    """Construct jwave Domain and Medium from input dict.

    use_fourier_series: wrap spatially-varying medium properties as FourierSeries
    (jwave Fields) rather than plain JAX arrays. Required for all solvers to avoid
    a jwave performance warning about Array inputs.
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

    # jwave Medium uses alpha in dB/(MHz^y * cm) — same convention as k-Wave
    medium = Medium(
        domain=domain,
        sound_speed=to_field(inp['sound_speed']),
        density=to_field(inp['density']),
        attenuation=to_field(inp['alpha_coeff']),
        pml_size=inp['pml_size'],
    )
    return domain, medium


def build_sources(inp, domain):
    """Reconstruct source mask and signals as jwave Sources."""
    from jwave.geometry import Sources
    import jax.numpy as jnp

    N = tuple(inp['N'].tolist())
    flat_indices = inp['mask_indices']                   # 0-based linear indices
    # MATLAB find() uses column-major (Fortran) linear indexing; use order='F'
    # so multi-dim indices match the (Nx, Ny, Nz) logical layout.
    multi_idx = np.unravel_index(flat_indices, N, order='F')  # tuple of arrays per dim
    positions = tuple(jnp.array(ax) for ax in multi_idx)

    signals = jnp.array(inp['source_p'], dtype=jnp.float32)  # [N_src x Nt]
    dt      = inp['dt']

    return Sources(positions=positions, signals=signals, dt=dt, domain=domain)


# ---------------------------------------------------------------------------
# Convergence diagnostics
# ---------------------------------------------------------------------------

def compute_spectral_residual(P_complex, inp, omega, source_field):
    """Approximate relative residual of the Helmholtz solution using a spectral Laplacian.

    Computes ||(∇² + k²(x))P - b|| / ||b|| where ∇² is the pseudo-spectral
    (periodic) Laplacian and k(x) = ω/c(x).  The PML is not included, so this
    is an approximation, but it reliably detects gross GMRES non-convergence
    where the residual is O(1) rather than O(1e-3) or smaller.

    Returns the relative residual norm (scalar, dimensionless).
    """
    N = tuple(inp['N'].tolist())
    dx = inp['dx']

    # Build wavenumber grids for each spatial dimension
    k_axes = [np.fft.fftfreq(n, d=d) * 2 * np.pi for n, d in zip(N, dx)]
    K2 = sum(np.broadcast_to(k[reshape_idx(i, len(N))].astype(np.float64), N) ** 2
             for i, k in enumerate(k_axes))

    # Spectral Laplacian: ∇²P = IFFT(-k² * FFT(P))
    P_fft = np.fft.fftn(P_complex)
    lap_P = np.fft.ifftn(-K2 * P_fft)

    # Local wavenumber squared k²(x) = (ω/c(x))²
    c = inp['sound_speed']
    if np.isscalar(c):
        k2_local = (omega / c) ** 2 * np.ones(N)
    else:
        k2_local = (omega / c) ** 2

    residual_field = lap_P + k2_local * P_complex - source_field
    b_norm = np.linalg.norm(source_field)
    if b_norm == 0:
        return float('nan')
    return float(np.linalg.norm(residual_field) / b_norm)


def reshape_idx(axis, ndim):
    """Return an index tuple that selects along `axis` with all other dims as singletons."""
    idx = [None] * ndim
    idx[axis] = slice(None)
    return tuple(idx)


def crop_inp(inp, bbox):
    """Crop inp dict to a 3-D bounding box (x0,y0,z0,x1,y1,z1), 0-based inclusive.

    Returns a new inp dict with updated N, medium fields, and source indices.
    Source elements outside the bbox are dropped.
    """
    x0, y0, z0, x1, y1, z1 = bbox
    N_orig = tuple(inp['N'].tolist())
    N_new  = np.array([x1 - x0 + 1, y1 - y0 + 1, z1 - z0 + 1], dtype=int)

    def crop_field(arr):
        if np.isscalar(arr):
            return arr
        return arr[x0:x1 + 1, y0:y1 + 1, z0:z1 + 1]

    flat_idx = inp['mask_indices']
    xi, yi, zi = np.unravel_index(flat_idx, N_orig, order='F')
    keep = (xi >= x0) & (xi <= x1) & (yi >= y0) & (yi <= y1) & (zi >= z0) & (zi <= z1)

    if not np.any(keep):
        raise ValueError('benchmark_fov_vox contains no source elements — check the bbox.')

    new_flat_idx = np.ravel_multi_index(
        (xi[keep] - x0, yi[keep] - y0, zi[keep] - z0),
        tuple(N_new.tolist()), order='F',
    )

    inp_c = dict(inp)
    inp_c['N']           = N_new
    inp_c['dx']          = inp['dx']
    inp_c['sound_speed'] = crop_field(inp['sound_speed'])
    inp_c['density']     = crop_field(inp['density'])
    inp_c['alpha_coeff'] = crop_field(inp['alpha_coeff'])
    inp_c['alpha_power'] = crop_field(inp['alpha_power'])
    inp_c['mask_indices'] = new_flat_idx
    inp_c['source_p']    = inp['source_p'][keep, :]
    return inp_c


def run_fov_benchmark(inp, bbox_str, gmres_restart, gmres_maxiter):
    """Crop the grid to a small FOV and compare helmholtz vs pseudo_spectral on that subgrid.

    Returns a dict of scalar diagnostic values to be written to /diagnostics/.
    """
    bbox = [int(v) for v in bbox_str.split(',')]
    inp_fov = crop_inp(inp, bbox)
    N_fov = inp_fov['N']
    n_src_fov = len(inp_fov['mask_indices'])
    full_vox = int(np.prod(inp['N']))
    fov_vox  = int(np.prod(N_fov))

    print(f'[run_jwave_sim] FOV benchmark: grid {inp["N"].tolist()} → {N_fov.tolist()} '
          f'({fov_vox/full_vox*100:.1f}% of original), N_src={n_src_fov}')

    print('[run_jwave_sim] FOV benchmark: running helmholtz ...')
    p_max_h, _, _, _, _ = run_helmholtz(inp_fov, gmres_restart, gmres_maxiter)

    print('[run_jwave_sim] FOV benchmark: running pseudo_spectral ...')
    p_max_td, _, _ = run_pseudo_spectral(inp_fov)

    peak_h  = float(p_max_h.max())
    peak_td = float(p_max_td.max())
    ratio   = peak_h / peak_td if peak_td > 0 else float('nan')

    hh_loc = np.unravel_index(p_max_h.argmax(),  p_max_h.shape)
    td_loc = np.unravel_index(p_max_td.argmax(), p_max_td.shape)
    peak_dist = float(np.sqrt(sum((a - b) ** 2 for a, b in zip(hh_loc, td_loc))))

    r = float(np.corrcoef(p_max_h.ravel(), p_max_td.ravel())[0, 1])

    print(f'[run_jwave_sim] FOV benchmark results:')
    print(f'  Spatial correlation  : {r:.4f}')
    print(f'  Peak ratio (HH/PS)   : {ratio:.3f}  (HH={peak_h:.0f} Pa  PS={peak_td:.0f} Pa)')
    print(f'  Peak distance        : {peak_dist:.1f} vox  '
          f'(HH={hh_loc}  PS={td_loc})')

    verdict = 'PASS' if r > 0.95 and 0.8 < ratio < 1.2 and peak_dist < 5 else 'WARN'
    print(f'[run_jwave_sim] FOV benchmark verdict: {verdict}')

    return dict(
        benchmark_correlation=r,
        benchmark_amplitude_ratio=ratio,
        benchmark_peak_dist_vox=peak_dist,
        benchmark_fov_vox=float(fov_vox),
    )


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

def run_pseudo_spectral(inp):
    """
    Run jwave pseudo-spectral time-domain (PSTD) simulation; return (p_max_all, p_min_all, p_final).

    Equivalent to kspaceFirstOrder* — iterates a semi-implicit first-order Euler
    integrator with pseudo-spectral (global FFT) spatial derivatives, and records
    running max/min over the last `record_window` steps.

    Known mismatch vs k-Wave: jwave Medium does not expose an alpha_power parameter,
    so the power-law exponent from the HDF5 input is ignored. Attenuation is applied
    with jwave's internal fixed exponent. Results will differ from k-Wave when
    alpha_power != the jwave default.

    Memory strategy: accumulate p_max/p_min in the scan carry so only 3 volumes
    (*N each) are live at any time, matching k-Wave's on-the-fly sensor behaviour.
    The full (*N, Nt) pressure time-series is never materialised.
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import TimeAxis
    from jwave.acoustics.time_varying import (
        TimeWavePropagationSettings,
        ongrid_wave_prop_params,
        momentum_conservation_rhs,
        mass_conservation_rhs,
        pressure_from_density,
    )

    domain, medium = build_medium_and_domain(inp, use_fourier_series=True)

    Nt = inp['Nt']
    dt = inp['dt']
    time_axis = TimeAxis(dt=dt, t_end=Nt * dt)
    sources   = build_sources(inp, domain)
    record_start = inp['record_start_index']  # 0-based

    print(f'[run_jwave_sim] Pseudo-spectral simulation: N={inp["N"]}, Nt={Nt}, dt={dt:.3e} s')
    print(f'[run_jwave_sim] Recording from step {record_start} to {Nt}')

    settings   = TimeWavePropagationSettings()
    wave_params = ongrid_wave_prop_params(medium, time_axis, settings=settings)
    c_ref   = wave_params['c_ref']
    pml_rho = wave_params['pml_rho']
    pml_u   = wave_params['pml_u']

    N = tuple(inp['N'].tolist())
    shape     = tuple(list(N) + [len(N)])
    shape_one = tuple(list(N) + [1])

    u0  = pml_u.replace_params(jnp.zeros(shape))
    p0  = pml_rho.replace_params(jnp.zeros(shape_one))
    rho0 = p0.replace_params(jnp.zeros(shape))

    neg_inf = jnp.full(N, -jnp.inf, dtype=jnp.float32)
    pos_inf = jnp.full(N,  jnp.inf, dtype=jnp.float32)

    def step(carry, n):
        p, u, rho, p_max, p_min = carry

        du  = momentum_conservation_rhs(p, u, medium, c_ref=c_ref, dt=dt)
        u   = pml_u  * (pml_u  * u   + dt * du)
        drho = mass_conservation_rhs(
            p, u, sources.on_grid(n), medium, c_ref=c_ref, dt=dt)
        rho  = pml_rho * (pml_rho * rho + dt * drho)
        p    = pressure_from_density(rho, medium)

        p_arr = p.params[..., 0]
        recording = n >= record_start
        p_max = jnp.where(recording, jnp.maximum(p_max, p_arr), p_max)
        p_min = jnp.where(recording, jnp.minimum(p_min, p_arr), p_min)

        return (p, u, rho, p_max, p_min), None  # no per-step output stack

    step = jax.checkpoint(step)

    init_carry = (p0, u0, rho0, neg_inf, pos_inf)
    (p_final_field, _, _, p_max_all, p_min_all), _ = jax.lax.scan(
        step, init_carry, jnp.arange(Nt))

    p_final = p_final_field.params[..., 0]
    return np.array(p_max_all), np.array(p_min_all), np.array(p_final)


def run_helmholtz(inp, gmres_restart=40, gmres_maxiter=300):
    """
    Run jwave frequency-domain (Helmholtz) solver.

    Solves for the complex pressure amplitude at freq_hz. Returns:
        p_max_all = |P|        (peak amplitude)
        p_min_all = -|P|       (negative amplitude, for MI compatibility)
        p_final   = Re(P)      (real part)
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import Domain, Medium
    from jwave.acoustics.time_harmonic import helmholtz_solver

    domain, medium = build_medium_and_domain(inp, use_fourier_series=True)

    freq_hz = inp['freq_hz']
    omega   = 2 * np.pi * freq_hz

    print(f'[run_jwave_sim] Helmholtz solver: N={inp["N"]}, freq={freq_hz/1e3:.1f} kHz, '
          f'GMRES restart={gmres_restart}, maxiter={gmres_maxiter}')

    # Build time-harmonic source amplitude from source signals
    # Use the first time step amplitude scaled to represent CW at freq_hz
    source_p = inp['source_p']
    # Extract complex amplitude per source element via DFT at freq_hz
    Nt = inp['Nt']
    dt = inp['dt']
    freqs     = np.fft.rfftfreq(Nt, dt)
    freq_idx  = int(np.argmin(np.abs(freqs - freq_hz)))
    freq_bin  = freqs[freq_idx]
    freq_err_pct = abs(freq_bin - freq_hz) / freq_hz * 100
    if freq_err_pct > 1.0:
        print(f'[run_jwave_sim] WARNING: requested freq {freq_hz/1e3:.3f} kHz, '
              f'nearest FFT bin is {freq_bin/1e3:.4f} kHz '
              f'(error={freq_err_pct:.2f}%). '
              'Increase Nt or adjust freq_hz for more accurate frequency matching.',
              file=sys.stderr)
    else:
        print(f'[run_jwave_sim] Freq bin: requested={freq_hz:.0f} Hz, '
              f'bin={freq_bin:.1f} Hz (err={freq_err_pct:.3f}%)')
    fft_coeff = np.fft.rfft(source_p, axis=-1)[:, freq_idx] * (2 / Nt)  # [N_src] complex

    N = tuple(inp['N'].tolist())
    source_field = np.zeros(N, dtype=complex)
    flat_idx     = inp['mask_indices']
    multi_idx    = np.unravel_index(flat_idx, N, order='F')
    source_field[multi_idx] = fft_coeff

    from jaxdf.discretization import FourierSeries
    # helmholtz_solver expects source as OnGrid (FourierSeries), not TimeHarmonicSource
    jw_source = FourierSeries(
        params=jnp.array(source_field[..., np.newaxis], dtype=jnp.complex64),
        domain=domain,
    )

    @jax.jit
    def forward(medium):
        return helmholtz_solver(medium, omega, jw_source,
                                restart=gmres_restart, maxiter=gmres_maxiter)

    P_raw = forward(medium)
    # FourierSeries stores data in .params; extract and squeeze the trailing singleton
    if hasattr(P_raw, 'params'):
        P_complex = np.array(P_raw.params[..., 0])
    else:
        P_complex = np.array(P_raw)

    p_amp     = np.abs(P_complex)
    p_max_all = p_amp
    p_min_all = -p_amp
    p_final   = np.real(P_complex)

    return p_max_all, p_min_all, p_final, P_complex, source_field


def build_medium_and_domain_fd(inp, accuracy=8):
    """Construct jwave Domain and Medium with FiniteDifferences field discretization.

    Identical to build_medium_and_domain(..., use_fourier_series=False) except
    medium fields are wrapped as JaxDF FiniteDifferences objects.  This causes
    jwave's operator dispatch (via plum multiple-dispatch) to select FD stencils
    for all spatial derivatives rather than the global FFT-based Laplacian used
    by the pseudo-spectral path.

    accuracy: stencil accuracy order (must be even; 2, 4, 6, 8, …).
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

    medium = Medium(
        domain=domain,
        sound_speed=to_fd_field(inp['sound_speed']),
        density=to_fd_field(inp['density']),
        attenuation=to_fd_field(inp['alpha_coeff']),
        pml_size=inp['pml_size'],
    )
    return domain, medium


def run_finite_differences(inp, accuracy=8):
    """
    Run jwave finite-difference time-domain (FDTD) simulation.

    Uses the same simulate_wave_propagation time-stepping loop as run_pseudo_spectral
    but wraps medium fields as JaxDF FiniteDifferences so that jwave's
    multiple-dispatch selects FD stencils (local, order-`accuracy`) rather than
    pseudo-spectral (global FFT) derivatives.

    Returns (p_max_all, p_min_all, p_final) matching the pseudo_spectral output format.
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import TimeAxis
    from jwave.acoustics.time_varying import (
        TimeWavePropagationSettings,
        ongrid_wave_prop_params,
        momentum_conservation_rhs,
        mass_conservation_rhs,
        pressure_from_density,
    )

    domain, medium = build_medium_and_domain_fd(inp, accuracy=accuracy)

    Nt = inp['Nt']
    dt = inp['dt']
    time_axis = TimeAxis(dt=dt, t_end=Nt * dt)
    sources   = build_sources(inp, domain)
    record_start = inp['record_start_index']

    print(f'[run_jwave_sim] FD time-domain simulation: N={inp["N"]}, Nt={Nt}, '
          f'dt={dt:.3e} s, FD accuracy order={accuracy}')
    print(f'[run_jwave_sim] Recording from step {record_start} to {Nt}')

    settings    = TimeWavePropagationSettings()
    wave_params = ongrid_wave_prop_params(medium, time_axis, settings=settings)
    c_ref   = wave_params['c_ref']
    pml_rho = wave_params['pml_rho']
    pml_u   = wave_params['pml_u']

    N = tuple(inp['N'].tolist())
    shape     = tuple(list(N) + [len(N)])
    shape_one = tuple(list(N) + [1])

    u0   = pml_u.replace_params(jnp.zeros(shape))
    p0   = pml_rho.replace_params(jnp.zeros(shape_one))
    rho0 = p0.replace_params(jnp.zeros(shape))

    neg_inf = jnp.full(N, -jnp.inf, dtype=jnp.float32)
    pos_inf = jnp.full(N,  jnp.inf, dtype=jnp.float32)

    def step(carry, n):
        p, u, rho, p_max, p_min = carry

        du   = momentum_conservation_rhs(p, u, medium, c_ref=c_ref, dt=dt)
        u    = pml_u  * (pml_u  * u   + dt * du)
        drho = mass_conservation_rhs(
            p, u, sources.on_grid(n), medium, c_ref=c_ref, dt=dt)
        rho  = pml_rho * (pml_rho * rho + dt * drho)
        p    = pressure_from_density(rho, medium)

        p_arr = p.params[..., 0]
        recording = n >= record_start
        p_max = jnp.where(recording, jnp.maximum(p_max, p_arr), p_max)
        p_min = jnp.where(recording, jnp.minimum(p_min, p_arr), p_min)

        return (p, u, rho, p_max, p_min), None

    step = jax.checkpoint(step)

    init_carry = (p0, u0, rho0, neg_inf, pos_inf)
    (p_final_field, _, _, p_max_all, p_min_all), _ = jax.lax.scan(
        step, init_carry, jnp.arange(Nt))

    p_final = p_final_field.params[..., 0]
    return np.array(p_max_all), np.array(p_min_all), np.array(p_final)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    try:
        print(f'[run_jwave_sim] Reading input: {args.input}')
        inp = read_input(args.input)

        diagnostics = {}

        if args.solver == 'pseudo_spectral':
            p_max_all, p_min_all, p_final = run_pseudo_spectral(inp)
        elif args.solver == 'finite_differences':
            p_max_all, p_min_all, p_final = run_finite_differences(inp, accuracy=args.fd_accuracy)
        else:
            p_max_all, p_min_all, p_final, P_complex, source_field = run_helmholtz(
                inp,
                gmres_restart=args.gmres_restart,
                gmres_maxiter=args.gmres_maxiter,
            )

            # ---- Convergence diagnostics ----
            omega = 2 * np.pi * inp['freq_hz']
            residual_norm = compute_spectral_residual(P_complex, inp, omega, source_field)
            print(f'[run_jwave_sim] Helmholtz spectral residual: {residual_norm:.4e}')
            diagnostics['residual_norm'] = residual_norm

            # ---- Optional FOV benchmark ----
            if args.benchmark_fov_vox:
                bm = run_fov_benchmark(inp, args.benchmark_fov_vox,
                                       args.gmres_restart, args.gmres_maxiter)
                diagnostics.update(bm)

        # Validate output shapes match declared grid dims before writing
        expected_shape = tuple(inp['N'].tolist())
        for arr, name in [(p_max_all, 'p_max_all'), (p_min_all, 'p_min_all'), (p_final, 'p_final')]:
            if arr.shape != expected_shape:
                raise ValueError(
                    f'Output {name} has shape {arr.shape} but expected grid N={expected_shape}. '
                    'Check solver output and dimension transpositions in run_jwave_sim.py.'
                )
        print(f'[run_jwave_sim] Output shapes validated: {expected_shape}')

        print(f'[run_jwave_sim] p_max={p_max_all.max():.1f} Pa  p_min={p_min_all.min():.1f} Pa')

        print(f'[run_jwave_sim] Writing output: {args.output}')
        write_output(args.output, p_max_all, p_min_all, p_final,
                     diagnostics=diagnostics if diagnostics else None)
        print('[run_jwave_sim] Done.')

    except Exception as e:
        import traceback
        print(f'[run_jwave_sim] ERROR: {e}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
