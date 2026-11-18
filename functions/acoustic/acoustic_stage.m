function [sensor_data, parameters, kgrid, kwave_medium, medium_masks, ...
          segmentation, source, source_labels, acoustic_provenance] = ...
    acoustic_stage(parameters, kgrid, kwave_medium, source, sensor, ...
                   medium_masks, segmentation, source_labels, ...
                   acoustic_provenance, planimg, trans_pos, focus_pos, pseudoCT)
% ACOUSTIC_STAGE  Run or load the acoustic pressure simulation.
%
% Handles a prior pressure map (path.acoustic_input), cache load/save, and backend
% dispatch (k-Wave default / jwave / BabelBrain).  Must be called after
% source setup and acoustic_water_baseline.
%
% When run_acoustic_sims is enabled and no valid cache exists, runs the
% simulation and caches the result.  When a cache exists and overwriting is
% not requested, loads the cached result.
%
% Use as:
%   [sensor_data, parameters, kgrid, kwave_medium, medium_masks, ...
%    segmentation, source, source_labels, acoustic_provenance] = ...
%       acoustic_stage(parameters, kgrid, kwave_medium, source, sensor, ...
%                      medium_masks, segmentation, source_labels, ...
%                      acoustic_provenance, planimg, trans_pos, focus_pos, pseudoCT)
%
% Inputs:
%   parameters        - pipeline parameters struct
%   kgrid             - kWaveGrid object ([] for BabelBrain backend)
%   kwave_medium      - medium properties struct
%   source            - k-Wave source struct ([] for BabelBrain backend)
%   sensor            - k-Wave sensor struct ([] for BabelBrain backend)
%   medium_masks      - tissue mask struct
%   segmentation      - segmentation volume
%   source_labels     - source label volume ([] for BabelBrain backend)
%   acoustic_provenance - provenance struct (FOV offsets, tx calibration)
%   planimg           - planning image struct with affine transforms
%   trans_pos         - transducer position in current sim-grid voxels
%   focus_pos         - focus position in current sim-grid voxels
%   pseudoCT          - pseudo-CT volume ([] if not available)
%
% Outputs:
%   sensor_data       - struct with p_max_all field (Pa) (also for a prior pressure
%                       map, resampled onto the grid), or [] when unavailable
%   parameters        - updated with state.acoustics_available and
%                       state.acoustic_dry_run flags
%   kgrid, kwave_medium, medium_masks, segmentation, source, source_labels -
%                       updated by k-Wave path (expanded for axisymmetric runs,
%                       also on cache reload); unchanged for other backends

arguments
    parameters        (1,1) struct
    kgrid
    kwave_medium      (1,1) struct
    source
    sensor
    medium_masks
    segmentation
    source_labels
    acoustic_provenance (1,1) struct
    planimg           (1,1) struct
    trans_pos
    focus_pos
    pseudoCT
end

if isfield(parameters.io, 'acoustic_cache_affix')
    acoustic_file_affix = parameters.io.acoustic_cache_affix;
else
    acoustic_file_affix = parameters.io.output_affix;
end
filename_sensor_data = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_results%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, acoustic_file_affix));

parameters.state.acoustics_available = 0;
parameters.state.acoustic_dry_run    = false;
parameters.state.recomputed.acoustic = false;   % set by RUN_ACOUSTIC_BACKEND (OUTPUTS_FOLLOW_RECOMPUTE)

% Settings the acoustic field depends on: written to the pressure map's
% sidecar (WRITE_ACOUSTIC_SIDECAR) and compared when a later run starts from it.
input_signature = cache_signature(parameters, 'acoustic');
acoustic_provenance.input_signature = input_signature;

% --- Prior pressure map (path.acoustic_input): no simulation -----------------
% The map is resampled onto this grid, so acoustic analysis and the thermal
% stage run as after a simulation (LOAD_ACOUSTIC_INPUT).
acoustic_input = acoustic_input_file(parameters);
if ~isempty(acoustic_input)
    [sensor_data, acoustic_provenance] = load_acoustic_input( ...
        acoustic_input, parameters, planimg, acoustic_provenance);
    parameters.state.acoustics_available = 1;
    return;
end

run_new_sim = isfield(parameters.modules, 'run_acoustic_sims') && ...
    parameters.modules.run_acoustic_sims && ...
    confirm_overwriting(filename_sensor_data, parameters) && ...
    (parameters.simulation.interactive == 0 || ...
     confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'));

if run_new_sim

    [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, ...
     source, source_labels, acoustic_provenance] = run_acoustic_backend( ...
        parameters, kgrid, kwave_medium, source, sensor, medium_masks, ...
        segmentation, source_labels, acoustic_provenance, planimg, ...
        trans_pos, focus_pos, pseudoCT, filename_sensor_data);
    if parameters.state.acoustic_dry_run
        % Inputs written, nothing simulated: there is no field to scale or analyse.
        parameters.modules.run_acoustic_analysis = 0;
        return;
    end
    sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);
    parameters.state.acoustics_available = 1;

elseif exist(filename_sensor_data, 'file')
    disp('Skipping acoustic simulation, loading existing output file.')
    live_provenance = acoustic_provenance;
    try
        [sensor_data, kgrid_cached, acoustic_provenance, ~, parameters] = ...
            cache_acoustic('load', filename_sensor_data, parameters, struct(), [], live_provenance);
        if user_driven_reuse(parameters)
            % Like a reused grid/medium cache: the cached field must come from
            % this run's settings (io.cache_signature_check).
            stored_sig = [];
            if isfield(acoustic_provenance, 'input_signature')
                stored_sig = acoustic_provenance.input_signature;
            end
            check_cache_signature(filename_sensor_data, stored_sig, input_signature, parameters);
        end
        acoustic_provenance.input_signature = input_signature;
        % The cached field of an axisymmetric run is already expanded: expand
        % the grid, positions and in-scope medium the same way.
        [parameters, kgrid, kwave_medium, medium_masks, segmentation, source, source_labels] = ...
            expand_axisymmetric_inputs(parameters, sensor_data, kgrid, kgrid_cached, kwave_medium, ...
                                       medium_masks, segmentation, source, source_labels, filename_sensor_data);
    catch ME
        if strcmp(ME.identifier, 'cache_acoustic:staleGrid') && cache_remap_enabled(parameters)
            % io.acoustic_cache_remap: map the cached field onto this run's grid
            % through world space instead of re-simulating or stopping.
            [sensor_data, kgrid, acoustic_provenance] = remap_stale_cache( ...
                filename_sensor_data, parameters, planimg, kgrid, live_provenance);
            acoustic_provenance.input_signature = input_signature;
        elseif strcmp(ME.identifier, 'cache_acoustic:staleGrid') && isempty(source)
            % No source to re-simulate with (e.g. a job that only reads a
            % combined field): the grid settings differ from the cache's run.
            differing = stale_cache_setting_diffs(filename_sensor_data, input_signature, parameters);
            error('acoustic_stage:staleCacheNoSource', ...
                ['%s\nThis run sets up no source (modules.run_source_setup = 0), so it cannot ' ...
                 're-simulate. %s'], ME.message, differing);
        elseif strcmp(ME.identifier, 'cache_acoustic:staleGrid')
            % Re-run with the configured backend (not unconditionally k-Wave),
            % scaled exactly like a fresh run.
            fprintf('[pipeline] Acoustic cache grid mismatch — re-running acoustic simulation (%s backend).\n', ...
                acoustic_backend_name(parameters));
            acoustic_provenance = live_provenance;   % not the stale cache's
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, ...
             source, source_labels, acoustic_provenance] = run_acoustic_backend( ...
                parameters, kgrid, kwave_medium, source, sensor, medium_masks, ...
                segmentation, source_labels, acoustic_provenance, planimg, ...
                trans_pos, focus_pos, pseudoCT, filename_sensor_data);
            sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);
        else
            rethrow(ME);
        end
    end
    parameters.state.acoustics_available = 1;
    % After axisymmetric expansion, the cache holds an expanded sensor_data but
    % source_labels has the pre-expansion shape.  Clear mismatched labels.
    if ~isempty(source_labels) && ~isequal(size(source_labels), size(sensor_data.p_max_all))
        warning('prestus_pipeline:sourceLabelsMismatch', ...
            ['source_labels size %s does not match cached sensor_data size %s. ' ...
             'Clearing source_labels — transducer bowl overlay absent. ' ...
             'Re-run with run_grid_setup=1 to regenerate.'], ...
            mat2str(size(source_labels)), mat2str(size(sensor_data.p_max_all)));
        source_labels = zeros(size(sensor_data.p_max_all), 'like', source_labels);
    end

else
    disp('No acoustic simulation available or requested ... skipping analysis')
    sensor_data = [];
    parameters.state.acoustics_available = 0;
    parameters.modules.run_acoustic_analysis = 0;
end

end

%% ------------------------------------------------------------------------
function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, ...
          source, source_labels, acoustic_provenance] = run_acoustic_backend( ...
        parameters, kgrid, kwave_medium, source, sensor, medium_masks, ...
        segmentation, source_labels, acoustic_provenance, planimg, ...
        trans_pos, focus_pos, pseudoCT, filename_sensor_data)
% RUN_ACOUSTIC_BACKEND  Dispatch one acoustic simulation to the configured
% backend (k-Wave default / jwave / BabelBrain) and cache the result.
% Shared by the fresh-run and stale-cache re-run paths so both honour
% parameters.acoustic_backend.

use_babelbrain = strcmpi(acoustic_backend_name(parameters), 'babelbrain');
use_jwave      = strcmpi(acoustic_backend_name(parameters), 'jwave');
% A dry run (jwave/BabelBrain) writes the solver inputs and returns no field;
% it is never cached, so a later real run cannot load the empty result.
parameters.state.acoustic_dry_run = is_acoustic_dry_run(parameters);
parameters.state.recomputed.acoustic = true;
save_cache = should_save_output(parameters.io, 'save_acoustic_matrices') && ...
    ~parameters.state.acoustic_dry_run;
% The grid the solver runs on. acoustic_wrapper returns a new, expanded grid
% without a time array after an axisymmetric run, so run_stats reads this one.
kgrid_sim = kgrid;
axisym_sim = isfield(parameters.grid, 'axisymmetric') && isequal(parameters.grid.axisymmetric, 1);
t_backend = tic;

if use_babelbrain
    sensor_data = babelbrain_wrapper( ...
        parameters, medium_masks, kwave_medium, ...
        focus_pos, trans_pos, planimg, pseudoCT);
    if isfield(sensor_data, 'babelbrain_version')
        acoustic_provenance.babelbrain_version = sensor_data.babelbrain_version;
        fprintf('[pipeline] BabelBrain %s  (%s)\n', ...
            sensor_data.babelbrain_version.version_txt, ...
            sensor_data.babelbrain_version.git_describe);
    end
    if isfield(sensor_data, 'freefield_isppa_wcm2') && isfinite(sensor_data.freefield_isppa_wcm2)
        acoustic_provenance.freefield_isppa_wcm2             = sensor_data.freefield_isppa_wcm2;
        acoustic_provenance.freefield_source                 = 'babelbrain_water_sim';
        acoustic_provenance.tx_calibration.isppa_wcm2        = sensor_data.freefield_isppa_wcm2;
        acoustic_provenance.tx_calibration.isppa_source      = 'babelbrain_water_sim';
        % A k-Wave baseline profile from acoustic_water_baseline is at the
        % unscaled PRESTUS drive level, not BabelBrain's; drop it.
        if isfield(acoustic_provenance, 'freefield_axial')
            acoustic_provenance = rmfield(acoustic_provenance, 'freefield_axial');
        end
    end
    sensor_data = backend_p_complex_files(parameters, sensor_data, filename_sensor_data, 'babelbrain', acoustic_provenance);
    if save_cache
        cache_acoustic('save', filename_sensor_data, parameters, sensor_data, [], ...
            acoustic_provenance, struct('parameters', parameters, 'backend', 'babelbrain'), kwave_medium);
        fprintf('[pipeline] babelbrain sensor_data cached: %s\n', filename_sensor_data);
    end

elseif use_jwave
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric
        error('[pipeline] jwave backend does not support axisymmetric grids.');
    end
    sensor_data = jwave_wrapper(parameters, kgrid, kwave_medium, source, source_labels, sensor);
    % Record the PML jwave actually used (its 'auto' differs from k-Wave's) so the
    % report states it, mirroring acoustic_wrapper for the k-Wave path.
    if isfield(sensor_data, 'jwave_pml_size')
        parameters.grid.pml_size_effective = sensor_data.jwave_pml_size;
    end
    sensor_data = backend_p_complex_files(parameters, sensor_data, filename_sensor_data, 'jwave', acoustic_provenance);
    if save_cache
        jw_solver = 'helmholtz';
        if isfield(parameters, 'jwave') && isfield(parameters.jwave, 'solver')
            jw_solver = parameters.jwave.solver;
        end
        cache_acoustic('save', filename_sensor_data, parameters, sensor_data, kgrid, ...
            acoustic_provenance, ...
            struct('parameters', parameters, 'backend', 'jwave', 'solver', jw_solver), ...
            kwave_medium);
        fprintf('[pipeline] jwave sensor_data cached: %s\n', filename_sensor_data);
    end

else
    % k-Wave (default)
    [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
        acoustic_wrapper(parameters, kgrid, kwave_medium, source, sensor, ...
                         medium_masks, filename_sensor_data, segmentation, ...
                         source_labels, acoustic_provenance);
end

% Size and duration for resource estimates (run_stats.jsonl); BabelBrain
% simulates on its own grid, so only its duration is recorded.
if ~parameters.state.acoustic_dry_run
    run_stats('add', 'acoustic_s', toc(t_backend));
    if ~use_babelbrain && isa(kgrid_sim, 'kWaveGrid')
        dims = [kgrid_sim.Nx, kgrid_sim.Ny, kgrid_sim.Nz];
        nt = kgrid_sim.Nt;
        if ~isnumeric(nt), nt = NaN; end   % no time array (frequency-domain jwave)
        pml = NaN;
        if isfield(parameters.grid, 'pml_size_effective')
            pml = max(double(parameters.grid.pml_size_effective(:)));
        end
        run_stats('set', 'acoustic_voxels', kgrid_sim.total_grid_points, 'acoustic_nt', nt, ...
            'acoustic_dims', dims(1:kgrid_sim.dim), 'acoustic_pml', pml, ...
            'axisymmetric', double(axisym_sim));
    end
end

end

%% ------------------------------------------------------------------------
function tf = cache_remap_enabled(parameters)
% io.acoustic_cache_remap on, for a 3-D (not axisymmetric) run
tf = isfield(parameters.io, 'acoustic_cache_remap') && ~isempty(parameters.io.acoustic_cache_remap) && ...
     logical(parameters.io.acoustic_cache_remap) && ...
     ~(isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric);
end

function [sensor_data, kgrid, prov] = remap_stale_cache(filename, parameters, planimg, kgrid_live, live_provenance)
% REMAP_STALE_CACHE  Load an acoustic cache written on another grid and map its
% field onto this run's grid (io.acoustic_cache_remap). Both grids are placed
% in world space by their planning images (GRID_WORLD_AFFINE: the cache run's
% planimg sidecar, this run's planimg); FOV boxes by their offsets. Pressures
% are resampled as p^2 (REMAP_ACOUSTIC_FIELD). The mapping, its scale and the
% fraction of this grid the cached field covers are recorded in
% acoustic_provenance.cache_remap and reported as a warning: the field was
% computed through the cache run's head model, which may differ from this
% run's near tissue boundaries.
p_load = parameters;
p_load.state.acoustic_cache_remap_load = true;
[sensor_data, kgrid_cached, prov, info_cached] = ...
    cache_acoustic('load', filename, p_load, struct(), [], live_provenance);

A_src = grid_world_affine(info_cached);
A_dst = grid_world_affine(planimg);
if isempty(A_src) || isempty(A_dst)
    error('acoustic_stage:remapNoAffine', ...
        ['io.acoustic_cache_remap: cannot place both grids in world space (the cache run''s planimg ' ...
         'sidecar or this run''s T1 header is missing):\n  %s'], filename);
end
[M, minfo] = grid_remap_affine(A_src, A_dst, fov_offset_of(prov), fov_offset_of(live_provenance));

src_dims = size(sensor_data.p_max_all);
src_dims = [src_dims ones(1, 3 - numel(src_dims))];
dst_dims = double(parameters.grid.dims(:)');
coverage = NaN;
for f = fieldnames(sensor_data)'
    v = sensor_data.(f{1});
    if ~(isnumeric(v) || islogical(v)) || ndims(v) < 3 || ~isequal(size(v, [1 2 3]), src_dims)
        continue
    end
    [sensor_data.(f{1}), rinfo] = remap_acoustic_field(v, M, dst_dims, field_quantity(f{1}, v));
    coverage = rinfo.coverage;
end

% k-Wave grid of this run: the cached one has the cache run's size (with the
% same PML padding, if it was stored padded); BabelBrain caches store none.
kgrid = kgrid_live;
if ~isempty(kgrid_cached) && isobject(kgrid_cached)
    pad = [kgrid_cached.Nx kgrid_cached.Ny kgrid_cached.Nz] - src_dims;
    N   = dst_dims + pad;
    dx  = parameters.grid.resolution_mm * 1e-3;
    kgrid = kWaveGrid(N(1), dx, N(2), dx, N(3), dx);
    try
        kgrid.setTime(kgrid_cached.Nt, kgrid_cached.dt);
    catch
    end
end

for fn = {'fov_offset_ac', 'full_ac_dims'}
    if isfield(live_provenance, fn{1})
        prov.(fn{1}) = live_provenance.(fn{1});
    elseif isfield(prov, fn{1})
        prov = rmfield(prov, fn{1});
    end
end
prov.cache_remap = struct('source_file', filename, 'M', M, 'scale', minfo.scale, ...
    'rigid', minfo.rigid, 'shift_vox', minfo.shift_vox, 'coverage', coverage, ...
    'source_dims', src_dims, 'dims', dst_dims);
warning('acoustic_stage:cacheRemapped', ...
    ['Acoustic cache from grid [%s] remapped onto this grid [%s] (io.acoustic_cache_remap): ' ...
     'shift [%s] voxels, scale [%s], %.1f %% of this grid covered. The field was computed ' ...
     'through the cache run''s head model:\n  %s'], num2str(src_dims), num2str(dst_dims), ...
    num2str(minfo.shift_vox, '%.2f '), num2str(minfo.scale, '%.4f '), 100 * coverage, filename);
end

function off = fov_offset_of(prov)
off = [1 1 1];
if isstruct(prov) && isfield(prov, 'fov_offset_ac') && ~isempty(prov.fov_offset_ac)
    off = double(prov.fov_offset_ac(:)');
end
end

function q = field_quantity(name, v)
% How REMAP_ACOUSTIC_FIELD resamples a sensor_data field
if islogical(v)
    q = 'mask';
elseif ~isreal(v)
    q = 'complex';
elseif startsWith(name, 'p_min')
    q = 'negative';
elseif startsWith(name, {'p_max', 'p_rms', 'p_amp', 'p_intensity_equiv', 'p_async'})
    q = 'pressure';                       % peaks and amplitudes (>= 0)
else
    q = 'linear';                         % e.g. p_final: a signed snapshot
end
end

function [parameters, kgrid, kwave_medium, medium_masks, segmentation, source, source_labels] = ...
        expand_axisymmetric_inputs(parameters, sensor_data, kgrid, kgrid_cached, kwave_medium, ...
                                   medium_masks, segmentation, source, source_labels, filename)
% EXPAND_AXISYMMETRIC_INPUTS  Bring the half-plane state of an axisymmetric
% run to the expanded state of the reloaded field, as ACOUSTIC_WRAPPER does
% after a simulation (ACOUSTIC_CONVERT_AXISYMMETRY). Variables that are not
% half-plane arrays (e.g. no source in a job that only reads the field) are
% returned unchanged. Non-axisymmetric runs only take the cached kgrid.

mode = axisymmetric_output_mode(parameters);
if isempty(mode) || ~isfield(parameters.grid, 'axisym_bilateral_dims') || ...
        ~isstruct(sensor_data) || ~isfield(sensor_data, 'p_max_all')
    kgrid = kgrid_cached;
    return
end

% The cached field decides the form; a 2-D plane cannot feed 3-D heating.
cached_mode = '2d';
if ndims(sensor_data.p_max_all) == 3; cached_mode = '3d'; end
if ~strcmp(cached_mode, mode)
    heating = isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims == 1;
    if heating && strcmp(cached_mode, '2d')
        error('acoustic_stage:axisymCacheForm', ...
            ['%s holds the 2-D axisymmetric plane, but heating needs the 3-D expansion. ' ...
             'Re-run the acoustic simulation with heating enabled or grid.axisymmetric_output: 3d.'], filename);
    end
    warning('acoustic_stage:axisymCacheForm', ...
        'Cached axisymmetric field is %s, this run would expand to %s; using the cached form.', ...
        upper(cached_mode), upper(mode));
    mode = cached_mode;
end

% Placeholders for variables that are not half-plane arrays, so the
% conversion updates parameters and kgrid in any case
half = parameters.grid.dims;
is_half = @(A) (isnumeric(A) || islogical(A)) && isequal(size(A), half);
real_seg   = is_half(segmentation);
real_masks = is_half(medium_masks);
real_src   = isstruct(source) && isfield(source, 'p_mask') && is_half(source.p_mask);
real_lbl   = is_half(source_labels);
seg_in = segmentation;   if ~real_seg,   seg_in = zeros(half, 'uint8'); end
msk_in = medium_masks;   if ~real_masks, msk_in = zeros(half, 'uint8'); end
src_in = source;         if ~real_src,   src_in = struct('p_mask', zeros(half)); end
lbl_in = source_labels;  if ~real_lbl,   lbl_in = zeros(half); end
med_in = kwave_medium;
for f = fieldnames(med_in)'
    if ~isscalar(med_in.(f{1})) && ~is_half(med_in.(f{1}))
        med_in = rmfield(med_in, f{1});   % not on the half-plane: leave as is
    end
end
sd_placeholder = struct('p_final', zeros(half, 'single'), 'p_max_all', zeros(half, 'single'));
kgrid_in = kgrid;
if isempty(kgrid_in); kgrid_in = struct(); end   % replaced by the conversion

[~, parameters, seg_out, msk_out, med_out, kgrid_out, src_out, lbl_out] = ...
    acoustic_convert_axisymmetry(parameters, sd_placeholder, seg_in, msk_in, med_in, ...
                                 kgrid_in, src_in, lbl_in, mode);

if real_seg,   segmentation  = seg_out; end
if real_masks, medium_masks  = msk_out; end
if real_src,   source        = src_out; end
if real_lbl,   source_labels = lbl_out; end
for f = fieldnames(med_out)'
    kwave_medium.(f{1}) = med_out.(f{1});
end
kgrid = kgrid_cached;
if ~isa(kgrid, 'kWaveGrid'); kgrid = kgrid_out; end

% Same response as a stale 3-D cache (CACHE_ACOUSTIC): the caller re-runs
% the simulation on its unchanged half-plane variables.
if ~isequal(size(sensor_data.p_max_all), parameters.grid.dims)
    error('cache_acoustic:staleGrid', ...
        ['Acoustic cache was built on a different grid [%s] than the expanded axisymmetric ' ...
         'grid [%s] of this run:\n  %s'], ...
        num2str(size(sensor_data.p_max_all)), num2str(parameters.grid.dims), filename);
end
end

%% ------------------------------------------------------------------------
function msg = stale_cache_setting_diffs(filename, current_sig, parameters)
% Names the settings in which this run differs from the run that wrote the
% acoustic cache (e.g. headmodel.csf_expansion changes the head crop).
msg = ['Check that its grid and head-model settings (e.g. grid.acoustic_fov_*, ' ...
       'headmodel.*) match the run that wrote the cache.'];
try
    s = load(filename, 'acoustic_provenance');
    stored = [];
    if isfield(s, 'acoustic_provenance') && isfield(s.acoustic_provenance, 'input_signature')
        stored = s.acoustic_provenance.input_signature;
    end
    if isempty(stored); return; end
    p = parameters;
    p.io.cache_signature_check = 'warn';
    ws = warning('off', 'prestus:cache:signatureMismatch');
    restore = onCleanup(@() warning(ws));
    [~, diffs] = check_cache_signature(filename, stored, current_sig, p);
    if ~isempty(diffs)
        msg = sprintf('Settings that differ from the run that wrote the cache: %s.', ...
            strjoin(diffs(1:min(end, 10)), ', '));
    end
catch
    % keep the generic advice: the cache's signature could not be read
end
end

%% ------------------------------------------------------------------------
function tf = user_driven_reuse(parameters)
% True when a run reloads its own acoustic cache (io.overwrite_files), not a
% field another job wrote for it: io.acoustic_cache_affix (multi-ISPPA,
% async combine) or the caches of a chained pipeline (io.reuse_upstream_caches).
tf = ~(isfield(parameters.io, 'acoustic_cache_affix') && ~isempty(parameters.io.acoustic_cache_affix)) && ...
     ~(isfield(parameters.io, 'reuse_upstream_caches') && ~isempty(parameters.io.reuse_upstream_caches));
end

function name = acoustic_backend_name(parameters)
% ACOUSTIC_BACKEND_NAME  Configured acoustic backend, defaulting to 'kwave'.
    name = 'kwave';
    if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
        name = char(parameters.acoustic_backend);
    end
end

function sensor_data = backend_p_complex_files(parameters, sensor_data, filename_sensor_data, backend, acoustic_provenance)
% io.save_p_complex for jwave and BabelBrain: the same p_complex file as the
% k-Wave path (WRITE_P_COMPLEX_FILES). jwave refers phases to t = 0 of its
% time-stepped simulation; BabelBrain to the start of its own recording window
% (one FFT bin over its last cycles), so SUPERPOSE_COMPLEX_FIELDS combines its
% fields only with each other.
    if ~(isfield(parameters.io, 'save_p_complex') && parameters.io.save_p_complex) || ...
            ~isfield(sensor_data, 'p_complex') || isempty(sensor_data.p_complex)
        return
    end
    meta = struct('phase_reference', 'simulation_t0', ...
        'freq_hz', parameters.transducer(1).freq_hz, 'field_names', {{'p_complex'}}, ...
        'record_t0_s', NaN, 'dt_s', NaN, 'method', '', 'n_samples', NaN, 'p_max_from_fit', false, ...
        'drive_level', 'simulated', 'grid_dims', [], 'fov_offset_ac', [], 'full_ac_dims', []);
    switch backend
        case 'jwave'
            solver = 'helmholtz';
            if isfield(parameters, 'jwave') && isfield(parameters.jwave, 'solver')
                solver = parameters.jwave.solver;
            end
            if strcmp(solver, 'helmholtz')
                meta.method = 'helmholtz';       % frequency-domain field; p_max_all is |P| as well
                meta.p_max_from_fit = true;
            else
                meta.method = 'steady_state_window';   % running least-squares sums over the window
            end
            if isfield(sensor_data, 'p_complex_window')
                w = sensor_data.p_complex_window;
                meta.record_t0_s = (double(w.record_start_index) - 1) * w.dt;
                meta.dt_s        = w.dt;
                meta.n_samples   = w.n_samples;
                sensor_data = rmfield(sensor_data, 'p_complex_window');
            end
        case 'babelbrain'
            meta.phase_reference = 'babelbrain_window';
            meta.method = 'babelbrain_fft';
            meta.p_max_from_fit = true;          % p_max_all is BabelBrain's p_amp
    end
    if isstruct(acoustic_provenance) && isfield(acoustic_provenance, 'fov_offset_ac')
        meta.fov_offset_ac = acoustic_provenance.fov_offset_ac;
        meta.full_ac_dims  = acoustic_provenance.full_ac_dims;
    end
    sensor_data.p_complex = cast(sensor_data.p_complex, 'single');
    sensor_data = write_p_complex_files(parameters, sensor_data, filename_sensor_data, meta);
end
