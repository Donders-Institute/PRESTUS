function [thermal_diff_obj, time_status_seq, results_heating] = thermal_simulation(...
    parameters, sensor_data, kgrid, kwave_medium, sensor, source, transf, medium_masks)
% THERMAL_SIMULATION  Simulate thermal effects of transcranial ultrasound using k-Wave
%
% Converts the acoustic pressure field (p_max_all from sensor_data) into a
% volumetric heat source Q via Q = alpha_np * p^2 / (rho * c). Runs
% kWaveDiffusion over the segments of thermal_schedule: PT ON, PT OFF,
% PTRI-OFF and post-PTRD cool-off. An async combined cache carries one scaled
% field per transducer (sensor_data.p_async_tx, [Nx Ny Nz N]); each gives its
% own Q_i, and a segment heats with the sum of the Q_i of the transducers on
% (timing.async_schedule).
% At each recorded step, temperature (T) and two CEM43 variants (k-Wave
% native; ISO 14708-3, see CEM43_ISO_INCREMENT) are accumulated into
% focal-plane snapshots, running maxima, and per-layer timeseries. GPU arrays are used when
% simulation.code_type is 'matlab_gpu' or 'cpp_gpu'.
%
% Use as:
%   [thermal_diff_obj, time_status_seq, results_heating] = ...
%       thermal_simulation(parameters, sensor_data, kgrid, kwave_medium, sensor, source, transf, medium_masks)
%
% Input:
%   parameters   - PRESTUS config; must contain thermal, timing, grid.dims,
%                  transducer(1).trans_pos, transducer(1).freq_hz [Hz],
%                  pct.enabled, simulation.code_type, simulation.precision,
%                  simulation.interactive, io.adopted_cem43 (optional).
%                  When calling with an independent thermal grid (via
%                  thermal_grid_setup), pass the parameters_th copy that
%                  has grid.dims and transducer(1).trans_pos updated to
%                  thermal grid coordinates.
%   sensor_data  - EITHER a struct from acoustic simulation containing
%                  p_max_all [Pa] (standard path; an async combined cache
%                  also holds p_async_tx [Pa], used instead), OR a char/string path
%                  to a NIfTI file (.nii or .nii.gz) whose volume will be
%                  read as p_max_all [Pa] and resampled to kgrid size via
%                  affine_resample_3d with transf.  The NIfTI must be in
%                  the same coordinate space as the T1 image (i.e. the space
%                  from which transf maps into the simulation grid).
%                  BabelBrain outputs pressure as p_amp [Pa] in a NIfTI that
%                  can be passed directly here.
%
%                  NOTE on units: the NIfTI must store peak positive pressure
%                  in Pascals (Pa), matching what PRESTUS writes to
%                  *_pressure.nii.gz.  BabelBrain also uses Pa.  If the file
%                  stores MPa, multiply after loading (not done automatically).
%   kgrid        - kWaveGrid object (may be from thermal_grid_setup for an
%                  independent thermal resolution/FOV)
%   kwave_medium - struct with fields: alpha_coeff, alpha_power,
%                  density [kg/m^3], sound_speed [m/s],
%                  thermal_conductivity [W/(m K)], specific_heat [J/(kg K)],
%                  perfusion_coeff [1/s], absorption_fraction or absorption_np
%                  [Np/m], temp_0 [°C]; temp_0 and the absorption map are the
%                  medium_plus fields
%                  reintegrated by the caller before this function is invoked.
%                  When using an independent thermal grid, pass the
%                  kwave_medium_th returned by thermal_grid_setup.
%   sensor       - not used: kWaveDiffusion gets no sensor, and all heating
%                  results come from the full temperature field. Kept for the
%                  call signature; any struct, e.g. struct('mask', ones(dims)).
%   source       - struct; Q is overwritten internally from p_max_all
%   transf       - affine forward transform (T1 → simulation grid, 4×4).
%                  Used (a) to resample NIfTI sensor_data into the grid, and
%                  (b) to resample adopted_cem43/heatmap NIfTIs.  Pass
%                  transf_th from thermal_grid_setup when using an independent
%                  thermal grid.
%   medium_masks - layer label map for per-tissue timeseries; must be sized
%                  to kgrid dimensions (use medium_masks_th from
%                  thermal_grid_setup when applicable)
%
% Output:
%   thermal_diff_obj  - kWaveDiffusion object after simulation
%   time_status_seq   - struct array with fields time [s], step, status, recorded
%   results_heating   - struct with fields:
%     maxT                 - maximum temperature over all time steps [°C]
%     focal_planeT         - temperature in focal plane [x, z, time] [°C]
%     heating_endT         - temperature at the final recorded step [°C]
%     CEM43                - running max CEM43 (k-Wave method) [min]
%     focal_planeCEM43     - CEM43 in focal plane [x, z, time] [min]
%     CEM43_end            - CEM43 at the final recorded step [min]
%     timeseries           - struct from THERMAL_UPDATE_TIMESERIES
%     CEM43_iso            - running max ISO CEM43 [min]
%     focal_planeCEM43_iso - ISO CEM43 in focal plane [x, z, time] [min]
%     CEM43_iso_end        - ISO CEM43 at the final recorded step [min]
%
% See also: THERMAL_PARAMETERS, THERMAL_PLOT_PROTOCOL, THERMAL_UPDATE_TIMESERIES,
%   THERMAL_ANALYSIS, MEDIUM_SETUP

arguments
    parameters   (1,1) struct
    sensor_data                    % struct with p_max_all [Pa], OR char/string NIfTI path
    kgrid        (1,1)
    kwave_medium (1,1) struct
    sensor       (1,1) struct
    source       (1,1) struct
    transf
    medium_masks {mustBeNumericOrLogical}
end

% A handoff that names a missing map stops the run: it would otherwise start
% from baseline temperature and zero dose and report them as cumulative.
for adopted = {'adopted_heatmap', 'adopted_cem43', 'adopted_cem43_iso'}
    if isfield(parameters.io, adopted{1}) && ~isempty(parameters.io.(adopted{1})) && ~isfile(parameters.io.(adopted{1}))
        error('thermal_simulation:adoptedMissing', ...
            'io.%s names a file that does not exist: %s (written by the NIfTI stage of the previous run).', ...
            adopted{1}, char(parameters.io.(adopted{1})));
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Notes:                                                            %
%                                                                   %
% Instead of modelling each oscillation, heating is simulated over  %
% time segments (THERMAL_SCHEDULE): ON and OFF steps within each    %
% pulse, then coarser steps between pulse trains and after PTRD.    %
% The steady-state acoustic field (what is seen in the figures and  %
% NIfTI files) is the heat source.                                  %
%                                                                   %
% For a detailed explanation on how to correctly configure your     %
% thermal parameters, see doc_simulations-thermal.md.               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Transducer position in grid coordinates
tr = parameters.transducer(1);

% [(pseudo-)CT] if 'k-plan' mapping is used for density:
% skull bone density is fixed for thermal simulation
% (note that downstream metrics such as sound speed contribute to pressure, but are not input to kWaveDiffusion)
% https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
if pct_skull_mapped(parameters) && ...
    (~isfield(parameters, 'pct') || ~isfield(parameters.pct, 'mapping_density') || strcmp(parameters.pct.mapping_density, 'k-plan'))
    skull_i = find(ismember(fieldnames(parameters.medium_properties), {'skull', 'skull_cortical', 'skull_trabecular'}));
    % [k-Plan] use fixed density for thermal simulation
    kwave_medium.density(ismember(medium_masks,skull_i)) = 1850;
    clear skull_i;
end

% Expand thermal-only fields from 2D axisymmetric to 3D when needed.
% acoustic_wrapper already expanded the acoustic fields (sound_speed, density,
% etc.); temp_0 and the absorption map were stored separately in medium_plus
% and were not expanded at that point.
if isfield(kwave_medium, 'absorption_np')
    absorption_field = 'absorption_np';         % any tissue with absorption model 'coefficient'
else
    absorption_field = 'absorption_fraction';   % all tissues with absorption model 'fraction'
end
if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
    Nlateral  = parameters.grid.axisym_bilateral_dims(2);
    ax_center = parameters.transducer(1).trans_pos_bilateral(2);
    kwave_medium.temp_0 = radialExpand2DTo3D(kwave_medium.temp_0, 'linear', Nlateral, ax_center);
    kwave_medium.(absorption_field) = radialExpand2DTo3D(kwave_medium.(absorption_field), 'linear', Nlateral, ax_center);
end

% Override initial temperature from a previous simulation's heatmap.
% Analogous to adopted_cem43 below.
if isfield(parameters.io, 'adopted_heatmap') && ~isempty(parameters.io.adopted_heatmap) && isfile(parameters.io.adopted_heatmap)
    heatmap_image = niftiread(parameters.io.adopted_heatmap);
    fprintf('\nAdopting heatmap %s from previous simulation\n', parameters.io.adopted_heatmap)
    kwave_medium.temp_0 = double(affine_resample_3d(heatmap_image, transf, ...
        size(medium_masks), 'cubic', parameters.thermal.temp_0.water));
end

if strcmp(absorption_field, 'absorption_np')
    % absorption coefficient, already in nepers/m at the source frequency (medium_setup)
    alpha_np = kwave_medium.absorption_np;
else
    % convert attenuation into absorption coefficients
    kwave_medium.alpha_coeff = ...
        kwave_medium.alpha_coeff .* kwave_medium.absorption_fraction;

    % convert the absorption coefficients to nepers/m (!)
    % see also fitPowerLawParamsMulti.m
    w = 2*pi*tr.freq_hz; % Frequency [rad/s]
    a0_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power); % [Nepers/((rad/s)^y m)]
    alpha_np = a0_np.*w.^kwave_medium.alpha_power;
    clear w a0_np;
end

% alternative simplified conversion dB to Nepers
% alpha_np = (100 * kwave_medium.alpha_coeff .* (tr.freq_hz/10^6)^kwave_medium.alpha_power)/8.686;

% save absorption coefficient [Np] for debugging
% (own name: medium_setup writes matrix_absorption_np before the acoustic FOV crop)
if contains(parameters.simulation.medium, {'layered', 'phantom'}) && debug_on(parameters, 'medium', 'volumes')
    try
        filename_absorption = fullfile(parameters.io.dir_debug_medium, 'matrix_absorption_np_thermal');
        niftiwrite(alpha_np, filename_absorption, 'Compressed',true);
        debug_register(parameters, 'medium', filename_absorption, ...
            'Absorption [Np/m] the thermal simulation uses, thermal grid', 'volume');
    catch
        warn("Error with saving absorption debug image...")
    end
end

% Get the maximum pressure (in Pa) and calculate Q, the volume rate of heat deposition
% sensor_data may be a struct (standard path) or a NIfTI file path (external acoustic)
if ischar(sensor_data) || isstring(sensor_data)
    nii_path = char(sensor_data);
    if ~isfile(nii_path)
        error('thermal_simulation:niftiNotFound', ...
            'sensor_data NIfTI not found: %s', nii_path);
    end
    fprintf('[thermal_simulation] Loading p_max_all from NIfTI: %s\n', nii_path);
    nii_vol  = single(niftiread(nii_path));   % [Pa], T1 space
    th_dims  = parameters.grid.dims;
    p = gather(max(single(affine_resample_3d(nii_vol, transf, th_dims, 'cubic', 0)), 0));
    clear nii_vol
elseif isfield(sensor_data, 'p_async_tx') && size(sensor_data.p_async_tx, 4) > 1
    % Async combined cache: one scaled field per transducer, each switched on
    % its own schedule (thermal_schedule).
    p = gather(abs(sensor_data.p_async_tx));     % stored as single
    if strcmp(char(parameters.simulation.precision), 'double')
        p = double(p);
    end
    if ~isequal(size(p, 1:3), size(alpha_np, 1:3))
        error('thermal_simulation:asyncGridMismatch', ...
            ['Per-transducer fields [%s] are not on the thermal grid [%s]. Async heating ' ...
             'needs the thermal grid to equal the acoustic grid (no thermal_resolution_mm / thermal_fov_mm).'], ...
            num2str(size(p, 1:3)), num2str(size(alpha_np, 1:3)));
    end
else
    p = gather(abs(sensor_data.p_max_all));
    if ~isscalar(alpha_np) && ~isequal(size(p), size(alpha_np))
        error('thermal_simulation:gridMismatch', ...
            ['p_max_all [%s] is not on the thermal grid [%s]: resample it with the medium ' ...
             '(thermal_stage does this for grid.thermal_resolution_mm / thermal_fov_mm).'], ...
            num2str(size(p)), num2str(size(alpha_np)));
    end
end
% Transducers at different frequencies: each field heats with the absorption
% at its own frequency, alpha(f_i) = alpha(f_1) * (f_i / f_1)^y with the tissue
% power-law exponent y (alpha_np above is at transducer(1)'s f_1). Async: one
% field per transducer; firing together: one per frequency
% (sensor_data.p_async_freqs_hz, COMBINE_ASYNC_INTENSITY).
freqs_tx = [parameters.transducer.freq_hz];
if isstruct(sensor_data) && isfield(sensor_data, 'p_async_freqs_hz') && ~isempty(sensor_data.p_async_freqs_hz)
    freqs_tx = double(sensor_data.p_async_freqs_hz(:)');
    % Firing groups: a source whose cache did not record its frequency takes its transducer's
    if isfield(sensor_data, 'p_async_source_tx') && numel(sensor_data.p_async_source_tx) == numel(freqs_tx)
        unknown = find(isnan(freqs_tx));
        for i_src = unknown
            freqs_tx(i_src) = double(parameters.transducer(sensor_data.p_async_source_tx(i_src)).freq_hz);
        end
        sensor_data.p_async_freqs_hz = freqs_tx;   % the schedule below uses the same values
    end
end
if size(p, 4) > 1 && numel(unique(freqs_tx)) > 1
    if size(p, 4) ~= numel(freqs_tx)
        error('thermal_simulation:asyncFrequencyCount', ...
            '%d heat-source fields but %d frequencies.', size(p, 4), numel(freqs_tx));
    end
    if ~isfield(kwave_medium, 'absorption_power')
        error('thermal_simulation:noAbsorptionPower', ...
            ['Transducers at different frequencies need the tissue absorption power law of the medium ' ...
             '(medium_plus.absorption_power): re-run the medium stage.']);
    end
    y = double(kwave_medium.absorption_power);
    alpha_tx = zeros(size(p), 'like', p);
    for i_tx = 1:numel(freqs_tx)
        alpha_tx(:, :, :, i_tx) = alpha_np .* (freqs_tx(i_tx) / tr.freq_hz) .^ y;
    end
    fprintf('[thermal_simulation] Absorption per field at its frequency (%s kHz).\n', ...
        num2str(freqs_tx / 1e3));
    alpha_np = alpha_tx;
    clear alpha_tx y
end
% Heat delivered to the system per source [W/m3]: Q_src(:,:,:,i)
Q_src = (alpha_np .* p.^2) ./ (kwave_medium.density .* kwave_medium.sound_speed);
n_src = size(Q_src, 4);
source.Q = sum(Q_src, 4);
source.T0 = kwave_medium.temp_0; %parameters.thermal.temp_0; % Initial temperature distribution
clear alpha_np p;

% remove unexpected fields from kWave_medium
kwave_medium = rmfield(kwave_medium, 'temp_0');
kwave_medium = rmfield(kwave_medium, intersect(fieldnames(kwave_medium), {'absorption_fraction', 'absorption_np', 'absorption_power'}));

% if perfusion was specified, implement it by specifying blood temperature
if ~all(isnan(kwave_medium.perfusion_coeff))
    kwave_medium.blood_ambient_temperature = 37; % [degC]
else
    kwave_medium = rmfield(kwave_medium, 'perfusion_coeff');
end

% create kWaveDiffusion object without a sensor: the heating results below are
% taken from the full temperature field (thermal_diff_obj.T, .cem43)
sensor = [];

% Check if DataCast is supported
try
    use_datacast = any(regexp(fileread(which('kWaveDiffusion')), '''DataCast'''));
    if ~use_datacast
        warn('MATLAB GPU support requested but kWaveDiffusion lacks ''DataCast'' (kWave ≥1.4.1 required).');
    end
catch
    warn('Cannot verify kWaveDiffusion DataCast support; assuming GPU-compatible version.');
    use_datacast = true;
end

% Set precision and enable GPU mode (if requested)
if use_datacast == true
    if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
        datacast = ['gpuArray-', char(parameters.simulation.precision)];
    else
        datacast = char(parameters.simulation.precision);
    end
end

% Build final input args
% PlotSim is always disabled here; when interactive we manage a single persistent figure
% ourselves via plotTemp after each takeTimeStep, because kWaveDiffusion.takeTimeStep
% creates a new figure handle on every call when PlotSim=true.
thermal_args = {'PlotSim', false};
if use_datacast
    thermal_args = [thermal_args, {'DataCast', datacast}];
end

log_thermal_updates = isfield(parameters, 'thermal') && isfield(parameters.thermal, 'log_thermal_updates') && parameters.thermal.log_thermal_updates;
is_first_thermal_step = true;

% Run simulation
thermal_diff_obj = kWaveDiffusion(...
    kgrid, ...
    kwave_medium, ...
    source, ...
    sensor, ...
    thermal_args{:});
source.Q = [];      % the schedule loop sets thermal_diff_obj.Q from Q_src

% In interactive mode, open a single figure that will be reused for all thermal updates
if parameters.simulation.interactive
    h_thermal_fig = figure('Name', 'kWave Thermal Simulation');
    thermal_diff_obj.plotTemp;
    drawnow;
end

% initialize field temperature
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
end
T_max = thermal_diff_obj.T;

% initialize field for cem43
if isfield(parameters.io, 'adopted_cem43') && ~isempty(parameters.io.adopted_cem43) && isfile(parameters.io.adopted_cem43)
    cumulative_heat_image = niftiread(parameters.io.adopted_cem43);
    fprintf('\nAdopting CEM43 heatmap %s from previous simulation\n', parameters.io.adopted_cem43)
    thermal_diff_obj.cem43 = double(affine_resample_3d(cumulative_heat_image, transf, ...
        size(medium_masks), 'nearest', 0));
else
    thermal_diff_obj.cem43 = zeros(size(thermal_diff_obj.T));
end
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    thermal_diff_obj.cem43 = gpuArray(thermal_diff_obj.cem43);
end
CEM43_max = thermal_diff_obj.cem43;

% initialize field for cem43 (iso variant)
if isfield(parameters.io, 'adopted_cem43_iso') && ~isempty(parameters.io.adopted_cem43_iso) && isfile(parameters.io.adopted_cem43_iso)
    cem43_iso_image = niftiread(parameters.io.adopted_cem43_iso);
    fprintf('\nAdopting ISO CEM43 heatmap %s from previous simulation\n', parameters.io.adopted_cem43_iso)
    tmp_obj.cem43_iso = double(affine_resample_3d(cem43_iso_image, transf, ...
        size(medium_masks), 'nearest', 0));
else
    tmp_obj.cem43_iso = zeros(size(thermal_diff_obj.T));
end
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    tmp_obj.cem43_iso = gpuArray(tmp_obj.cem43_iso);
end
CEM43_iso_max = tmp_obj.cem43_iso;

% Convert the thermal parameters for the simulation (validates the timing block)
params_thermal = thermal_parameters(parameters);


% Firing schedule: segments in which a fixed set of sources is on. Firing groups
% (FIRING_GROUPS): each field follows its group's windows (p_async_source_tx).
src_opts = {};
if isstruct(sensor_data) && isfield(sensor_data, 'p_async_source_tx') && ...
        numel(sensor_data.p_async_source_tx) == n_src
    src_opts = {'SourceTx', double(sensor_data.p_async_source_tx(:)'), ...
                'SourceFreq', double(sensor_data.p_async_freqs_hz(:)')};
end
sched = thermal_schedule(parameters, n_src, src_opts{:});

% Ramped pulses (timing.pulse_ramp_dur > 0): each source's heat times its ramp
% factor, from the pd of the transducer whose windows it follows (THERMAL_RAMP_FACTORS)
ramp_eta = thermal_ramp_factors(params_thermal, sched);
for i = find(ramp_eta ~= 1)
    Q_src(:,:,:,i) = Q_src(:,:,:,i) * ramp_eta(i);
    fprintf('[thermal_simulation] Pulse ramp (%s, %.3g ms per edge, PD %.3g ms): heat source %d x %.3f\n', ...
        params_thermal.ramp_shape, params_thermal.ramp_dur * 1e3, sched.source(i).pd * 1e3, i, ramp_eta(i));
end

% Visualize the thermal parameters
thermal_plot_protocol(params_thermal, parameters, 'Schedule', sched);

n_seg = numel(sched.t0);
fprintf('[thermal_simulation] ---- thermal simulation summary ----\n');
% 2D or 3D: print as many dimensions as the grid has
grid_spacing_mm = [kgrid.dx, kgrid.dy, kgrid.dz] * 1e3;
fprintf('[thermal_simulation]   grid dims      : %s\n', strjoin(string(parameters.grid.dims), ' x '));
fprintf('[thermal_simulation]   resolution     : %s mm\n', ...
    strjoin(compose('%.2f', grid_spacing_mm(1:numel(parameters.grid.dims))), ' x '));
if n_src > 1
    fprintf('[thermal_simulation]   heat sources   : %d (async schedule ''%s''; offsets %s s)\n', ...
        n_src, sched.mode, mat2str([sched.source.offset_s], 4));
end
fprintf('[thermal_simulation]   ON segments    : %d\n', nnz(any(sched.active, 2)));
fprintf('[thermal_simulation]   OFF segments   : %d\n', nnz(~any(sched.active, 2)));
fprintf('[thermal_simulation]   heating dur    : %.1f ms\n', sum(sched.dur(any(sched.active, 2)))*1e3);
fprintf('[thermal_simulation]   total sim time : %.1f ms\n', sched.t_end*1e3);
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    fprintf('[thermal_simulation]   compute mode   : GPU (%s)\n', parameters.simulation.code_type);
else
    fprintf('[thermal_simulation]   compute mode   : CPU (%s)\n', parameters.simulation.code_type);
end
fprintf('[thermal_simulation]   medium type    : %s\n', parameters.simulation.medium);
fprintf('[thermal_simulation] ------------------------------------\n');

% One recorded timepoint per segment, plus the initial state
total_timepoints = 1 + n_seg;

% Total kWaveDiffusion.takeTimeStep steps across the whole run (used for the
% ETA below). kWaveDiffusion's own built-in estimate (external/k-wave)
% requires >=10 steps within a single takeTimeStep call to fire, which the
% short per-pulse ON/OFF calls here rarely reach — so PRESTUS times the first
% call itself and extrapolates across calls instead.
total_thermal_steps = sched.total_steps;

% Focal-plane frames (x-z plane through the transducer in 3D, the whole grid in
% 2D): three arrays of one frame per timepoint. Long protocols on large planes
% keep every k-th segment (plus the initial and final state) within
% thermal.max_focal_plane_gb; the per-tissue time series and the max/end maps
% use every segment.
max_plane_gb = 2;
if isfield(parameters, 'thermal') && isfield(parameters.thermal, 'max_focal_plane_gb') ...
        && ~isempty(parameters.thermal.max_focal_plane_gb)
    max_plane_gb = parameters.thermal.max_focal_plane_gb;
end
if ndims(T_max) == 3
    plane_numel = size(T_max, 1) * size(T_max, 3);
else
    plane_numel = numel(T_max);
end
[plane_keep, plane_stride] = thermal_plane_frames(total_timepoints, plane_numel, max_plane_gb);
n_plane = nnz(plane_keep);
if plane_stride > 1
    fprintf(['[thermal_simulation]   focal plane    : every %d. segment (%d of %d frames; ' ...
        'thermal.max_focal_plane_gb = %g)\n'], plane_stride, n_plane, total_timepoints, max_plane_gb);
end

% Resize focal arrays
if ndims(T_max) == 3
    T_focal         = NaN([size(T_max,[1,3]) n_plane]);
    CEM43_focal     = NaN([size(T_max,[1,3]) n_plane]);
    CEM43_iso_focal = NaN([size(T_max,[1,3]) n_plane]);
    T_focal(:,:,1)         = squeeze(T_max(:,tr.trans_pos(2),:));
    CEM43_focal(:,:,1)     = squeeze(CEM43_max(:,tr.trans_pos(2),:));
    CEM43_iso_focal(:,:,1) = squeeze(CEM43_iso_max(:,tr.trans_pos(2),:));
elseif ndims(T_max) == 2
    T_focal         = NaN([size(T_max) n_plane]);
    CEM43_focal     = NaN([size(CEM43_max) n_plane]);
    CEM43_iso_focal = NaN([size(CEM43_iso_max) n_plane]);
    T_focal(:,:,1)         = T_max;
    CEM43_focal(:,:,1)     = CEM43_max;
    CEM43_iso_focal(:,:,1) = CEM43_iso_max;
end

% Per-layer masks for the time series, built once (on the GPU for GPU runs):
% rebuilding them at every segment cost seconds per segment on large grids
on_gpu = any(strcmp(parameters.simulation.code_type, {'matlab_gpu', 'cpp_gpu'}));
ts_masks = tissuemask_binary(parameters, medium_masks);
for f = fieldnames(ts_masks)'
    if nnz(ts_masks.(f{1})) == 0
        ts_masks = rmfield(ts_masks, f{1});
    elseif on_gpu
        ts_masks.(f{1}) = gpuArray(ts_masks.(f{1}));
    end
end
timeseries = thermal_update_timeseries(parameters, ts_masks);

% Heat source for each set of source weights that is on together (computed once).
% Weights are 1 for a source that is on; with thermal.pulse_model 'train_averaged'
% a pulse train is one segment and its source weighs PD/PRI (THERMAL_SCHEDULE).
[patterns, ~, pattern_of_seg] = unique(sched.weight, 'rows');
Q_pattern = cell(size(patterns, 1), 1);
for j = 1:size(patterns, 1)
    w = patterns(j, :);
    if any(w)
        if n_src == 1 && w == 1
            Q_pattern{j} = Q_src;       % shared copy, no extra memory
        else
            Q_pattern{j} = 0;
            for s = find(w)
                Q_pattern{j} = Q_pattern{j} + w(s) * Q_src(:,:,:,s);
            end
        end
    end
end
clear Q_src;

% Train-averaged segments last seconds: step them at the diffusion time step
% CFL * dx^2 / max diffusivity (CFL 0.1, as k-Plan's thermal solver)
if strcmp(sched.pulse_model, 'train_averaged')
    D = kwave_medium.thermal_conductivity ./ (kwave_medium.density .* kwave_medium.specific_heat);
    spacing = [kgrid.dx, kgrid.dy, kgrid.dz];
    dx = min(spacing(1:kgrid.dim));
    dt_max = 0.1 * dx^2 / double(gather(max(D(:))));
    avg = sched.fine;
    sched.n_steps(avg) = max(1, ceil(sched.dur(avg) / dt_max));
    sched.step_dur = sched.dur ./ sched.n_steps;
    sched.total_steps = sum(sched.n_steps);
    total_thermal_steps = sched.total_steps;
    fprintf('[thermal_simulation]   pulse model    : train_averaged (step %.3f s, %d steps)\n', ...
        dt_max, sched.total_steps);
end

cur_timepoint = 1;
cur_plane = 1;
seq_status = repmat({'off'}, 1, n_seg);
seq_status(any(sched.active, 2)) = {'on'};
cur_pattern = 0;
report_every = max(1, floor(n_seg / 10));

fprintf('Starting %d thermal segments (%.1f s sonication + %.1f s post)\n', ...
    n_seg, sched.t_sonication, sched.t_end - sched.t_sonication);

for seg_i = 1:n_seg
    if mod(seg_i, report_every) == 0
        fprintf('[thermal_simulation]   segment %d/%d (t = %.2f s)\n', seg_i, n_seg, sched.t0(seg_i));
    end

    % Heat source: the sum of the sources on during this segment
    j = pattern_of_seg(seg_i);
    if j ~= cur_pattern
        if isempty(Q_pattern{j})
            thermal_diff_obj.Q(:,:,:) = 0;
        else
            thermal_diff_obj.Q = Q_pattern{j};
        end
        cur_pattern = j;
    end

    n_steps_i = sched.n_steps(seg_i);
    step_dur_i = sched.step_dur(seg_i); %#ok<NASGU> (used in evalc)
    if log_thermal_updates || is_first_thermal_step
        if is_first_thermal_step
            first_step_timer = tic;
        end
        thermal_diff_obj.takeTimeStep(n_steps_i, sched.step_dur(seg_i));
        if is_first_thermal_step
            % k-Wave-style ETA: time a known number of steps, extrapolate
            % linearly to the full run (analogous to kWaveDiffusion's own
            % num_time_steps_before_simulation_time_estimate, but computed
            % across takeTimeStep calls rather than within a single one).
            per_step_s = toc(first_step_timer) / max(n_steps_i, 1);
            fprintf('[thermal_simulation]   estimated total thermal sim time: %s (%d steps @ %.3f s/step)\n', ...
                scaleTime(per_step_s * total_thermal_steps), total_thermal_steps, per_step_s);
            if ~log_thermal_updates
                fprintf('  [Subsequent k-Wave thermal outputs suppressed. Set parameters.thermal.log_thermal_updates=1 to enable.]\n');
            end
        end
        is_first_thermal_step = false;
    else
        evalc('thermal_diff_obj.takeTimeStep(n_steps_i, step_dur_i)');
    end
    if parameters.simulation.interactive && ishandle(h_thermal_fig)
        figure(h_thermal_fig); thermal_diff_obj.plotTemp; drawnow;
    end

    % CEM43 iso update — R uses T_max so that any voxel that has previously
    % exceeded 43 °C permanently accumulates at R=0.5, even after cooling.
    % Irreversible damage (T>=57 ever reached) is recorded as Inf. See doc_simulations-thermal.md.
    tmp_obj.cem43_iso = tmp_obj.cem43_iso + ...
        cem43_iso_increment(thermal_diff_obj.T, T_max, sched.dur(seg_i));

    % Update focal / max
    T_cur         = thermal_diff_obj.T;
    CEM43_cur     = thermal_diff_obj.cem43;
    CEM43_iso_cur = tmp_obj.cem43_iso;
    cur_timepoint = cur_timepoint + 1;
    if plane_keep(cur_timepoint)
        cur_plane = cur_plane + 1;
        if ndims(T_max) == 3
            T_focal(:,:,cur_plane)         = squeeze(T_cur(:,tr.trans_pos(2),:));
            CEM43_focal(:,:,cur_plane)     = squeeze(CEM43_cur(:,tr.trans_pos(2),:));
            CEM43_iso_focal(:,:,cur_plane) = squeeze(CEM43_iso_cur(:,tr.trans_pos(2),:));
        else
            T_focal(:,:,cur_plane)         = T_cur;
            CEM43_focal(:,:,cur_plane)     = CEM43_cur;
            CEM43_iso_focal(:,:,cur_plane) = CEM43_iso_cur;
        end
    end
    T_max         = max(T_max,         T_cur);
    CEM43_max     = max(CEM43_max,     CEM43_cur);
    CEM43_iso_max = max(CEM43_iso_max, CEM43_iso_cur);
    timeseries = thermal_update_timeseries(parameters, ts_masks, timeseries, T_cur, CEM43_cur, CEM43_iso_cur);
end

% Protocol on/off sequence for thermal_plot_sim: the initial state, then one
% entry per segment (time at its end, cumulative steps); 'plane' marks the
% entries that have a focal-plane frame
time_status_seq = struct( ...
    'status',   [{'off'}, seq_status], ...
    'time',     num2cell([0, (sched.t0 + sched.dur)']), ...
    'step',     num2cell([0, cumsum(sched.n_steps)']), ...
    'recorded', num2cell(ones(1, n_seg + 1)), ...
    'plane',    num2cell(double(plane_keep)));

% Trim unused timepoints (if too many were originally required)
T_focal         = T_focal(:,:,1:cur_plane);
CEM43_focal     = CEM43_focal(:,:,1:cur_plane);
CEM43_iso_focal = CEM43_iso_focal(:,:,1:cur_plane);

% Gather results into output struct (gathering from GPU if needed)
results_heating.maxT                 = gather(T_max);
results_heating.focal_planeT         = gather(T_focal);
results_heating.heating_endT         = gather(T_cur);
results_heating.CEM43                = gather(CEM43_max);
results_heating.focal_planeCEM43     = gather(CEM43_focal);
results_heating.CEM43_end            = gather(CEM43_cur);
results_heating.timeseries           = timeseries;
results_heating.CEM43_iso            = gather(CEM43_iso_max);
results_heating.focal_planeCEM43_iso = gather(CEM43_iso_focal);
results_heating.CEM43_iso_end        = gather(CEM43_iso_cur);

fprintf('Thermal simulation complete. Recorded %d timepoints.\n', cur_timepoint);

end