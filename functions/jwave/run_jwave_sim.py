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
        --solver time_domain

Solvers:
    time_domain  - Pseudo-spectral time-domain propagation via
                   jwave.acoustics.time_varying.simulate_wave_propagation.
                   Equivalent to kspaceFirstOrder*. Returns p_max_all,
                   p_min_all, p_final over the recorded window.
    helmholtz    - Frequency-domain GMRES solver via
                   jwave.acoustics.time_harmonic.helmholtz_solver.
                   More efficient for CW simulations. Returns p_max_all
                   (amplitude), p_min_all (negative amplitude), p_final
                   (real part of complex pressure). Requires freq_hz in
                   the HDF5 metadata attribute.

Exit codes:
    0  — success
    1  — error (message printed to stderr)
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
    p.add_argument('--solver', default='time_domain',
                   choices=['time_domain', 'helmholtz'],
                   help='Solver backend (default: time_domain)')
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
        # MATLAB writes source.p as [N_src x Nt] in column-major; h5py reads
        # it transposed as (Nt, N_src). Transpose back to (N_src, Nt).
        source_p            = f['/source/p'][()].T        # [N_src x Nt]
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


def write_output(h5_path, p_max_all, p_min_all, p_final):
    """Write pressure field outputs for MATLAB to read.

    Arrays are transposed before writing so that MATLAB reads them with the
    correct (Nx, Ny, Nz) logical layout (HDF5/h5py reverses axis order relative
    to MATLAB's column-major convention).
    """
    def to_matlab(arr):
        return np.ascontiguousarray(arr.T).astype(np.float32)

    with h5py.File(h5_path, 'w') as f:
        f.create_dataset('/pressure/p_max_all', data=to_matlab(p_max_all))
        f.create_dataset('/pressure/p_min_all', data=to_matlab(p_min_all))
        f.create_dataset('/pressure/p_final',   data=to_matlab(p_final))


# ---------------------------------------------------------------------------
# Medium / geometry helpers
# ---------------------------------------------------------------------------

def build_medium_and_domain(inp, use_fourier_series=False):
    """Construct jwave Domain and Medium from input dict.

    use_fourier_series: wrap spatially-varying medium properties as FourierSeries
    (jwave Fields) rather than plain JAX arrays. Required for the Helmholtz solver
    to avoid a jwave performance warning; left False for time-domain until tested.
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
# Solvers
# ---------------------------------------------------------------------------

def run_time_domain(inp):
    """
    Run jwave time-domain simulation and return (p_max_all, p_min_all, p_final).

    Equivalent to kspaceFirstOrder* — iterates the pseudo-spectral symplectic
    integrator and records running max/min over the last `record_window` steps.

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

    domain, medium = build_medium_and_domain(inp)

    Nt = inp['Nt']
    dt = inp['dt']
    time_axis = TimeAxis(dt=dt, t_end=Nt * dt)
    sources   = build_sources(inp, domain)
    record_start = inp['record_start_index']  # 0-based

    print(f'[run_jwave_sim] Time-domain simulation: N={inp["N"]}, Nt={Nt}, dt={dt:.3e} s')
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

    return p_max_all, p_min_all, p_final


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    try:
        print(f'[run_jwave_sim] Reading input: {args.input}')
        inp = read_input(args.input)

        if args.solver == 'time_domain':
            p_max_all, p_min_all, p_final = run_time_domain(inp)
        else:
            p_max_all, p_min_all, p_final = run_helmholtz(
                inp,
                gmres_restart=args.gmres_restart,
                gmres_maxiter=args.gmres_maxiter,
            )

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
        write_output(args.output, p_max_all, p_min_all, p_final)
        print('[run_jwave_sim] Done.')

    except Exception as e:
        import traceback
        print(f'[run_jwave_sim] ERROR: {e}', file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
