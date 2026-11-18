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
        if strcmp(ME.identifier, 'cache_acoustic:staleGrid') && isempty(source)
            % No source to re-simulate with (e.g. a job that only reads a
            % combined field): the grid settings differ from the cache's run.
            error('acoustic_stage:staleCacheNoSource', ...
                ['%s\nThis run sets up no source (modules.run_source_setup = 0), so it cannot ' ...
                 're-simulate. Check that its grid settings (e.g. grid.acoustic_fov_*) match ' ...
                 'the run that wrote the cache.'], ME.message);
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
    if ~use_babelbrain && isa(kgrid, 'kWaveGrid')
        run_stats('set', 'acoustic_voxels', kgrid.total_grid_points, 'acoustic_nt', kgrid.Nt);
    end
end

end

%% ------------------------------------------------------------------------
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
