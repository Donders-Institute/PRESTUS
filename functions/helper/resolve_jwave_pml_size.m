function pml = resolve_jwave_pml_size(pml_size, freq_hz, c_max, dx_m, pml_range)
% RESOLVE_JWAVE_PML_SIZE  Resolve grid.pml_size 'auto' to a wavelength-based scalar PML thickness for jwave
%
% k-Wave's getOptimalPMLSize (used via resolve_pml_size.m for the k-Wave
% backend) picks PML thickness purely to minimize FFT prime factors for
% k-Wave's own compiled backend, searching a FIXED absolute range
% (default [10, 40] grid points) regardless of frequency or resolution. That
% heuristic is irrelevant to jwave: JAX/XLA's FFT backend does not expose the
% same prime-factorization sensitivity, and jwave.geometry.Medium.pml_size is
% a single SCALAR applied uniformly to every axis (no per-axis PML at all).
%
% What actually matters for jwave is PML thickness IN WAVELENGTHS. jwave's own
% official example notebooks (docs/notebooks/harmonic/*.ipynb in the jwave
% repo) consistently use ~2.4-2.9 wavelengths of PML regardless of absolute
% grid size. Using a fixed absolute point range instead means 'auto' can end
% up a disproportionately large fraction of a small domain (wasting GMRES
% iterations/memory every step) or, at very fine resolutions where the point
% range clips low, insufficiently thick for the wavelength in play.
%
% c_max, not c_min: jwave's Helmholtz PML (complex_pml_on_grid) is a fixed
% function of voxel position and does not scale with omega or the local sound
% speed, so its absorption in wavelengths is lowest where it borders the
% fastest medium (e.g. skull). Sizing from c_max (the longest wavelength)
% gives every bordering medium at least TARGET_WAVELENGTHS of PML. This is
% not the same "conservative" as the PPW check, which uses c_min.
% See doc_jwave.md (PML thickness).
%
% Use as:
%   pml = resolve_jwave_pml_size(pml_size, freq_hz, c_max, dx_m)
%   pml = resolve_jwave_pml_size(pml_size, freq_hz, c_max, dx_m, pml_range)
%
% Input:
%   pml_size  - parameters.grid.pml_size (numeric, char, or string). Numeric
%               values are returned unchanged (explicit user override).
%   freq_hz   - transducer centre frequency [Hz]
%   c_max     - maximum sound speed in the medium [m/s] (longest wavelength,
%               i.e. the most conservative PML coverage requirement -- NOT
%               the same "conservative" as the PPW/spatial-sampling check,
%               which instead wants the shortest wavelength / c_min)
%   dx_m      - grid resolution [m] (assumed isotropic, matching jwave's
%               single-scalar PML thickness)
%   pml_range - [min, max] grid-point bounds (default = [10, 40], matching
%               k-Wave's getOptimalPMLSize bounds for familiarity)
%
% Output:
%   pml - scalar PML thickness in grid points
%
% See also: RESOLVE_PML_SIZE, GETOPTIMALPMLSIZE

arguments
    pml_size
    freq_hz    (1,1) double
    c_max      (1,1) double
    dx_m       (1,1) double
    pml_range  (1,2) double = [10, 40]
end

TARGET_WAVELENGTHS = 2.5;   % matches jwave's own example notebooks (~2.4-2.9 lambda)

if (ischar(pml_size) || isstring(pml_size)) && strcmp(pml_size, 'auto')
    wavelength_m = c_max / freq_hz;
    ppw = wavelength_m / dx_m;
    pml = ceil(TARGET_WAVELENGTHS * ppw);
    pml = min(max(pml, pml_range(1)), pml_range(2));
else
    pml = double(pml_size);
end

end
