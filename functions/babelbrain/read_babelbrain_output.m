function sensor_data = read_babelbrain_output(h5_path, use_refocus)
% READ_BABELBRAIN_OUTPUT  Read BabelBrain DataForSim.h5 into PRESTUS sensor_data
%
% Reads the peak pressure field from a BabelBrain *DataForSim.h5 output file
% and packages it into the sensor_data struct expected by acoustic_analysis
% and thermal_simulation.
%
% BabelBrain HDF5 output keys (BabelIntegrationBASE.py Step 9/10):
%   /p_amp          - [Nx x Ny x Nz] float32, peak positive pressure [Pa]
%   /p_amp_refocus  - same, after refocusing back-propagation (optional)
%   /XDim, /YDim, /ZDim - voxel coordinates [m]
%
% Use as:
%   sensor_data = read_babelbrain_output(h5_path)
%   sensor_data = read_babelbrain_output(h5_path, true)  % use refocused map
%
% Input:
%   h5_path      - path to *DataForSim.h5 produced by RUN_SIM().RunCases()
%   use_refocus  - (optional, default false) use /p_amp_refocus if available
%
% Output:
%   sensor_data  - struct with field p_max_all [Nx x Ny x Nz] double [Pa],
%                  matching the format returned by acoustic_wrapper.
%                  Also contains .babelbrain_source = h5_path for provenance.
%
% See also: BABELBRAIN_WRAPPER, ACOUSTIC_WRAPPER

arguments
    h5_path     (1,:) char
    use_refocus (1,1) logical = false
end

if ~isfile(h5_path)
    error('read_babelbrain_output: file not found: %s', h5_path);
end

% Select pressure key
p_key = '/p_amp';
if use_refocus
    info = h5info(h5_path);
    ds_names = {info.Datasets.Name};
    if any(strcmp(ds_names, 'p_amp_refocus'))
        p_key = '/p_amp_refocus';
        fprintf('[read_babelbrain_output] Using refocused pressure map.\n');
    else
        warning('read_babelbrain_output: /p_amp_refocus not found; falling back to /p_amp.');
    end
end

p_amp = double(h5read(h5_path, p_key));   % [Nx x Ny x Nz], Pa

% BabelBrain stores with PML stripped and zero-padded back to input mask dims.
% No further cropping needed; pass through directly.
sensor_data.p_max_all        = p_amp;
sensor_data.babelbrain_source = h5_path;

% Try to read spatial coordinates for provenance / sanity check
try
    xdim = h5read(h5_path, '/XDim');
    ydim = h5read(h5_path, '/YDim');
    zdim = h5read(h5_path, '/ZDim');
    dx_mm = mean(diff(xdim)) * 1e3;
    fprintf('[read_babelbrain_output] Grid: [%d %d %d]  dx=%.3f mm\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3), dx_mm);
    sensor_data.babelbrain_xdim = xdim;
    sensor_data.babelbrain_ydim = ydim;
    sensor_data.babelbrain_zdim = zdim;
catch
    fprintf('[read_babelbrain_output] Grid: [%d %d %d]\n', ...
        size(p_amp,1), size(p_amp,2), size(p_amp,3));
end

p_max_pa = max(p_amp(:));
fprintf('[read_babelbrain_output] p_max = %.1f Pa (%.2f MPa)\n', p_max_pa, p_max_pa/1e6);

end
