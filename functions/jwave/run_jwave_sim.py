#!/usr/bin/env python3
"""
run_jwave_sim.py — PRESTUS → jwave acoustic simulation bridge

Called by jwave_wrapper.m via system(). Reads simulation inputs from an HDF5
file written by MATLAB, runs the jwave Helmholtz solver, and writes pressure
field outputs (p_max_all, p_min_all, p_final) to a second HDF5 file that
MATLAB reads back.

Usage:
    python run_jwave_sim.py \\
        --input  /path/to/sub-001_layered_jwave_input.h5 \\
        --output /path/to/sub-001_layered_jwave_output.h5

Solver:
    helmholtz  - Frequency-domain GMRES solver via
                 jwave.acoustics.time_harmonic.helmholtz_solver.
                 Solves the heterogeneous Helmholtz equation with a Fourier
                 spectral collocation method and PML absorbing boundaries.
                 Returns p_max_all (amplitude), p_min_all (negative amplitude),
                 p_final (real part of complex pressure).
                 Attenuation is applied via db2neper with fixed exponent y=2.0.
                 Requires freq_hz in the HDF5 metadata attribute.

Exit codes:
    0  — success
    1  — error (message printed to stderr)

Diagnostics written to /diagnostics/ in the output HDF5:
    residual_norm  — relative spectral residual ||ΔP||/||b|| (scalar)
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
    p.add_argument('--gmres_restart', type=int, default=40,
                   help='GMRES restart parameter (default: 40; JAX default is 20)')
    p.add_argument('--gmres_maxiter', type=int, default=300,
                   help='Max GMRES outer iterations (default: 300)')
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

def build_medium_and_domain(inp):
    """Construct jwave Domain and Medium from input dict."""
    from jwave.geometry import Domain, Medium
    import jax.numpy as jnp

    N  = tuple(inp['N'].tolist())
    dx = tuple(inp['dx'].tolist())
    domain = Domain(N, dx)

    def to_field(val):
        if np.isscalar(val):
            return float(val)
        arr = jnp.array(val, dtype=jnp.float32)
        from jaxdf.discretization import FourierSeries
        return FourierSeries(params=arr[..., np.newaxis], domain=domain)

    # jwave Medium uses alpha in dB/(MHz^y * cm) — same convention as k-Wave.
    # The Helmholtz wavevector operator converts via db2neper with fixed y=2.0.
    medium = Medium(
        domain=domain,
        sound_speed=to_field(inp['sound_speed']),
        density=to_field(inp['density']),
        attenuation=to_field(inp['alpha_coeff']),
        pml_size=inp['pml_size'],
    )
    return domain, medium


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


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

def run_helmholtz(inp, gmres_restart=40, gmres_maxiter=300):
    """
    Run jwave frequency-domain (Helmholtz) solver.

    Solves the heterogeneous Helmholtz equation at freq_hz using GMRES with
    Fourier spectral collocation and PML absorbing boundaries. Attenuation is
    applied via db2neper(alpha_coeff, y=2.0) inside jwave's wavevector operator.

    Returns:
        p_max_all = |P|        (peak amplitude)
        p_min_all = -|P|       (negative amplitude, for MI compatibility)
        p_final   = Re(P)      (real part)
        P_complex              (full complex pressure field)
        source_field           (complex source array, for residual diagnostics)
    """
    import jax
    import jax.numpy as jnp
    from jwave.geometry import Domain, Medium
    from jwave.acoustics.time_harmonic import helmholtz_solver

    domain, medium = build_medium_and_domain(inp)

    freq_hz = inp['freq_hz']
    omega   = 2 * np.pi * freq_hz

    print(f'[run_jwave_sim] Helmholtz solver: N={inp["N"]}, freq={freq_hz/1e3:.1f} kHz, '
          f'GMRES restart={gmres_restart}, maxiter={gmres_maxiter}')

    # Extract complex source amplitude per element via DFT at freq_hz
    source_p = inp['source_p']
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
    jw_source = FourierSeries(
        params=jnp.array(source_field[..., np.newaxis], dtype=jnp.complex64),
        domain=domain,
    )

    @jax.jit
    def forward(medium):
        return helmholtz_solver(medium, omega, jw_source,
                                restart=gmres_restart, maxiter=gmres_maxiter)

    P_raw = forward(medium)
    if hasattr(P_raw, 'params'):
        P_complex = np.array(P_raw.params[..., 0])
    else:
        P_complex = np.array(P_raw)

    p_amp     = np.abs(P_complex)
    p_max_all = p_amp
    p_min_all = -p_amp
    p_final   = np.real(P_complex)

    return p_max_all, p_min_all, p_final, P_complex, source_field


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    try:
        print(f'[run_jwave_sim] Reading input: {args.input}')
        inp = read_input(args.input)

        p_max_all, p_min_all, p_final, P_complex, source_field = run_helmholtz(
            inp,
            gmres_restart=args.gmres_restart,
            gmres_maxiter=args.gmres_maxiter,
        )

        # Convergence diagnostics
        omega = 2 * np.pi * inp['freq_hz']
        residual_norm = compute_spectral_residual(P_complex, inp, omega, source_field)
        print(f'[run_jwave_sim] Helmholtz spectral residual: {residual_norm:.4e}')
        diagnostics = {'residual_norm': residual_norm}

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
        write_output(args.output, p_max_all, p_min_all, p_final, diagnostics=diagnostics)
        print('[run_jwave_sim] Done.')

    except Exception as e:
        import traceback
        print(f'[run_jwave_sim] ERROR: {e}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
