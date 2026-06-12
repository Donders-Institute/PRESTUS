function sensor_data = read_babelbrain_output(h5_path, use_refocus, parameters, planimg)
% READ_BABELBRAIN_OUTPUT  Read BabelBrain DataForSim.h5 into PRESTUS sensor_data
%
% Reads the peak pressure field from a BabelBrain *DataForSim.h5 output file,
% resamples it from BabelBrain's internal simulation grid to the PRESTUS acoustic
% grid when the two resolutions differ, and packages the result into the
% sensor_data struct expected by acoustic_analysis and thermal_simulation.
%
% BabelBrain runs at a grid step derived from its PPW setting:
%   dx_bb = c / (freq x ppw)
% This may differ from parameters.grid.resolution_mm (the PRESTUS acoustic grid
% step). When they differ, p_amp is resampled onto the PRESTUS grid using the
% same affine_resample_3d approach used by read_babelbrain_thermal. When
% parameters / planimg are not supplied, no resampling is performed (legacy
% behaviour).
%
% BabelBrain HDF5 output keys (BabelIntegrationBASE.py Step 9/10):
%   /p_amp          - [Nx x Ny x Nz] float32, peak positive pressure [Pa]
%   /p_amp_refocus  - same, after refocusing back-propagation (optional)
%   /p_complex      - [Nx x Ny x Nz] complex64, steady-state Fourier pressure [Pa]
%                     stored by h5py as a compound type with fields 'r' and 'i'
%   /p_complex_refocus - same, after refocusing (optional)
%   /x_vec, /y_vec, /z_vec - voxel centre coordinates [m] (cropped simulation grid)
%
% Use as:
%   sensor_data = read_babelbrain_output(h5_path)
%   sensor_data = read_babelbrain_output(h5_path, true)
%   sensor_data = read_babelbrain_output(h5_path, use_refocus, parameters, planimg)
%
% Input:
%   h5_path      - path to *DataForSim.h5 produced by RUN_SIM().RunCases()
%   use_refocus  - (optional, default false) use /p_amp_refocus if available
%   parameters   - (optional) PRESTUS config struct; required for resampling
%                  (uses grid.resolution_mm and grid.dims)
%   planimg      - (optional) planimg struct from grid_tissue_setup; required
%                  for affine-correct resampling in layered (subject) mode
%
% Output:
%   sensor_data  - struct with fields:
%                    .p_max_all   [Nx x Ny x Nz] double [Pa], peak pressure
%                    .p_complex   [Nx x Ny x Nz] complex double [Pa], Fourier
%                                 pressure (amplitude + phase); present when
%                                 /p_complex exists in the HDF5 file
%                  Both are at parameters.grid.resolution_mm resolution when
%                  resampled. Also contains .babelbrain_source = h5_path.
%
% See also: BABELBRAIN_WRAPPER, ACOUSTIC_WRAPPER, READ_BABELBRAIN_THERMAL

arguments
    h5_path     (1,:) char
    use_refocus (1,1) logical = false
    parameters  (1,1) struct  = struct()
    planimg               = []
end

ensure_hdf5_plugin_path();

if ~isfile(h5_path)
    error('read_babelbrain_output: file not found: %s', h5_path);
end

% Inventory datasets in the file (needed for refocus fallback and p_complex check)
info     = h5info(h5_path);
ds_names = {info.Datasets.Name};

% Select amplitude pressure key
if ~any(strcmp(ds_names, 'p_amp'))
    error('read_babelbrain_output:missingDataset', ...
        '/p_amp dataset not found in BabelBrain output:\n  %s', h5_path);
end
p_key = '/p_amp';
if use_refocus
    if any(strcmp(ds_names, 'p_amp_refocus'))
        p_key = '/p_amp_refocus';
        fprintf('[read_babelbrain_output] Using refocused pressure map.\n');
    else
        warning('read_babelbrain_output: /p_amp_refocus not found; falling back to /p_amp.');
    end
end

p_amp = double(h5read(h5_path, p_key));   % [Nx x Ny x Nz], Pa

% Read spatial coordinate vectors
% If parameters are supplied (resampling may be needed), fail hard on missing coords
% rather than silently skipping the resample step with a corrupt output.
has_params = isfield(parameters, 'grid') ...
             && isfield(parameters.grid, 'resolution_mm') ...
             && isfield(parameters.grid, 'dims');

xdim = []; ydim = []; zdim = [];
coord_read_ok = false;
try
    % BabelBrain saves coordinate vectors as /x_vec /y_vec /z_vec [m]
    xdim = h5read(h5_path, '/x_vec');
    ydim = h5read(h5_path, '/y_vec');
    zdim = h5read(h5_path, '/z_vec');
    coord_read_ok = true;
    dx_bb_mm = mean(diff(xdim)) * 1e3;
    dy_bb_mm = mean(diff(ydim)) * 1e3;
    dz_bb_mm = mean(diff(zdim)) * 1e3;
    fprintf('[read_babelbrain_output] BabelBrain grid: [%d %d %d]  dx=%.4f dy=%.4f dz=%.4f mm\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3), dx_bb_mm, dy_bb_mm, dz_bb_mm);
    sensor_data.babelbrain_x_vec = xdim;
    sensor_data.babelbrain_y_vec = ydim;
    sensor_data.babelbrain_z_vec = zdim;
catch ME
    dx_bb_mm = NaN;
    p_dims = size(p_amp);
    if has_params
        expected = parameters.grid.dims(:)';
        if isequal(p_dims, expected)
            warning(['read_babelbrain_output: coordinate vectors /x_vec /y_vec /z_vec not found in:\n  %s\n' ...
                'BabelBrain output dims [%d %d %d] match the PRESTUS grid — skipping resampling.\n' ...
                'Delete the cached .h5 file and re-run to get coordinate outputs.'], ...
                h5_path, p_dims(1), p_dims(2), p_dims(3));
        else
            error(['read_babelbrain_output: coordinate vectors /x_vec /y_vec /z_vec not found in:\n  %s\n' ...
                'Error: %s\n' ...
                'BabelBrain output dims [%d %d %d] do not match PRESTUS grid [%d %d %d].\n' ...
                'Delete the cached .h5 and re-run to enable resampling.'], ...
                h5_path, ME.message, p_dims(1), p_dims(2), p_dims(3), expected(1), expected(2), expected(3));
        end
    end
    fprintf('[read_babelbrain_output] BabelBrain grid: [%d %d %d]  (no coordinate vectors; no resampling)\n', ...
        p_dims(1), p_dims(2), p_dims(3));
end

% Resample onto PRESTUS acoustic grid when resolutions differ
needs_resample = has_params && coord_read_ok && ~isnan(dx_bb_mm) && ~isempty(xdim) ...
                 && (abs(dx_bb_mm - parameters.grid.resolution_mm) > 1e-3 ...
                  || abs(dy_bb_mm - parameters.grid.resolution_mm) > 1e-3 ...
                  || abs(dz_bb_mm - parameters.grid.resolution_mm) > 1e-3);

if needs_resample
    res_mm = parameters.grid.resolution_mm;
    fprintf('[read_babelbrain_output] Resampling %.4f mm -> %.4f mm  (target dims: [%d %d %d])\n', ...
        dx_bb_mm, res_mm, parameters.grid.dims(1), parameters.grid.dims(2), parameters.grid.dims(3));
    x_mm = xdim(:)' * 1e3;
    y_mm = ydim(:)' * 1e3;
    z_mm = zdim(:)' * 1e3;
    p_amp = resample_bb_to_prestus(p_amp, x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm);
    fprintf('[read_babelbrain_output] Resampled grid: [%d %d %d]\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3));
    % Hard check: resampled dims must match the declared PRESTUS grid
    act = size(p_amp); exp = parameters.grid.dims(:)';
    if ~isequal(act, exp)
        error(['read_babelbrain_output: resampled p_amp size [%d %d %d] does not match ' ...
            'expected PRESTUS grid dims [%d %d %d].\n' ...
            'Check planimg affine transforms and BabelBrain coordinate vectors.'], ...
            act(1), act(2), act(3), exp(1), exp(2), exp(3));
    end
elseif has_params && coord_read_ok && ~isnan(dx_bb_mm)
    fprintf('[read_babelbrain_output] Grid resolutions match (%.4f mm); no resampling needed.\n', dx_bb_mm);
end

sensor_data.p_max_all         = p_amp;
sensor_data.babelbrain_source = h5_path;

p_max_pa = max(p_amp(:));
fprintf('[read_babelbrain_output] p_max = %.1f Pa (%.2f MPa)\n', p_max_pa, p_max_pa/1e6);

% Read complex Fourier pressure field when present.
% h5py stores numpy complex64 as a compound HDF5 type with fields 'r' (real)
% and 'i' (imaginary). MATLAB's h5read returns a struct with .r and .i.
cx_key = '/p_complex';
if use_refocus && any(strcmp(ds_names, 'p_complex_refocus'))
    cx_key = '/p_complex_refocus';
end
if any(strcmp(ds_names, cx_key(2:end)))
    raw = h5read(h5_path, cx_key);
    if isstruct(raw) && isfield(raw, 'r') && isfield(raw, 'i')
        p_cx = double(raw.r) + 1i * double(raw.i);
    else
        p_cx = double(raw);  % already numeric complex (uncommon)
    end
    if needs_resample
        p_cx = resample_bb_to_prestus(real(p_cx), x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm) + ...
            1i * resample_bb_to_prestus(imag(p_cx), x_mm, y_mm, z_mm, planimg, parameters.grid.dims, res_mm);
    end
    sensor_data.p_complex = p_cx;
end

end

%% Local helpers

function vol_out = resample_bb_to_prestus(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm)
% Resample a BabelBrain 3D volume (RAS mm coordinates) onto the PRESTUS
% acoustic grid via affine_resample_3d.
%
% Mirrors the same helper in read_babelbrain_thermal.m.

    dx_bb = mean(diff(x_mm));
    dy_bb = mean(diff(y_mm));
    dz_bb = mean(diff(z_mm));

    % BabelBrain RAS mm -> BabelBrain voxel (1-based)
    origin = [x_mm(1), y_mm(1), z_mm(1)];
    M_ras_to_bb = diag([1/dx_bb, 1/dy_bb, 1/dz_bb, 1]);
    M_ras_to_bb(1:3, 4) = (-[origin(1)/dx_bb, origin(2)/dy_bb, origin(3)/dz_bb] + 1)';

    if ~isempty(planimg) && isstruct(planimg) ...
            && isfield(planimg, 't1_header') && ~isempty(planimg.t1_header) ...
            && isfield(planimg, 'transf')
        % Layered / subject mode: use T1 affine chain
        T_ras_to_t1vox    = inv(planimg.t1_header.Transform.T);  % correct inverse; T' would swap rotation and translation blocks
        T_t1vox_to_acgrid = planimg.transf;
        M_bb_to_acgrid    = T_t1vox_to_acgrid * T_ras_to_t1vox * inv(M_ras_to_bb);
    else
        % Phantom / water mode: assume grids share RAS origin, scale by ratio
        M_bb_to_acgrid = diag([dx_bb/res_mm, dx_bb/res_mm, dx_bb/res_mm, 1]);
    end

    vol_out = double(affine_resample_3d(single(vol_bb), M_bb_to_acgrid, prestus_dims, 'linear', 0));
end
