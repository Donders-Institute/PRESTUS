function [sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
        cache_acoustic(mode, filename, parameters, sensor_data, kgrid, ...
                       acoustic_provenance, acoustic_info, kwave_medium)
% CACHE_ACOUSTIC  Save or load the slim acoustic results cache.
%
% Centralises all cache I/O for the acoustic stage so that
% acoustic_wrapper and prestus_pipeline stay free of cache bookkeeping.
%
% SAVE mode
%   Writes a slim cache containing only the variables that cannot be
%   reconstructed from the grid/medium caches or from parameters:
%     sensor_data, kgrid, acoustic_provenance, acoustic_info
%
%   Async exception: combine_async_intensity reads kwave_medium directly
%   from each per-transducer cache file, so it is included in async runs.
%
% LOAD mode
%   Reads the slim variables from the cache file.  For legacy full caches
%   the extra fields are silently ignored; the caller's in-scope variables
%   (from the grid/medium/source stages) are always preferred.
%   Calls apply_isppa_scaling so the caller receives a correctly scaled
%   field. Axisymmetric runs: the cached field is already expanded; the
%   caller (ACOUSTIC_STAGE) expands parameters and the grid/medium variables.
%   When the cached free-field baseline is missing or unreadable, the
%   baseline of the current run (acoustic_provenance passed on load) is
%   used for that scaling.
%
% Use as:
%   % save
%   cache_acoustic('save', filename, parameters, sensor_data, kgrid, ...
%                  acoustic_provenance, acoustic_info, kwave_medium)
%
%   % load
%   [sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
%       cache_acoustic('load', filename, parameters)
%   % load, falling back to this run's free-water baseline
%   [...] = cache_acoustic('load', filename, parameters, struct(), [], live_provenance)
%
% Input:
%   mode              - 'save' or 'load'
%   filename          - full path to the cache .mat file
%   parameters        - PRESTUS config struct
%   sensor_data       - (save only) struct with p_max_all [Pa]
%   kgrid             - (save only) kWaveGrid object
%   acoustic_provenance - save: struct with freefield_isppa_wcm2;
%                       load: (optional) provenance of the current run, used
%                       when the cached baseline is missing or unreadable
%   acoustic_info     - (save only) struct with parameters snapshot
%   kwave_medium      - (save only) medium properties; included in async mode
%
% Output (load mode only):
%   sensor_data       - loaded and ISPPA-scaled pressure struct
%   kgrid             - kWaveGrid object
%   acoustic_provenance - provenance struct
%   acoustic_info     - info struct
%   parameters        - parameters (unchanged on load)

arguments
    mode                (1,:) char   {mustBeMember(mode, {'save','load'})}
    filename            (1,:) char
    parameters          (1,1) struct
    sensor_data               struct = struct()
    kgrid                            = []
    acoustic_provenance       struct = struct()
    acoustic_info             struct = struct()
    kwave_medium              struct = struct()
end

if strcmp(mode, 'save')
    acoustic_cache_save(filename, parameters, sensor_data, kgrid, ...
                        acoustic_provenance, acoustic_info, kwave_medium);
    return
end

[sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
    acoustic_cache_load(filename, parameters, acoustic_provenance);

end

% =========================================================================

function acoustic_cache_save(filename, parameters, sensor_data, kgrid, ...
                              acoustic_provenance, acoustic_info, kwave_medium)

% include kwave_medium when explicitly requested.
% async_transducer_pipeline sets save_kwave_medium=1 for per-transducer runs
% so combine_async_intensity has access to per-voxel density/sound_speed.
% Record the grid dimensions at save time so the load path can detect a
% stale cache (e.g. grid was regenerated after the acoustic run).
acoustic_grid_dims = parameters.grid.dims; %#ok<NASGU>

if should_save_output(parameters.io, 'save_kwave_medium')
    save(filename, ...
        'sensor_data', ...
        'kgrid', ...
        'kwave_medium', ...
        'acoustic_info', ...
        'acoustic_provenance', ...
        'acoustic_grid_dims', '-v7.3')
else
    save(filename, ...
        'sensor_data', ...
        'kgrid', ...
        'acoustic_info', ...
        'acoustic_provenance', ...
        'acoustic_grid_dims', '-v7.3')
end

end

% =========================================================================

function [sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
        acoustic_cache_load(filename, parameters, live_provenance)

% Load only the slim variables the cache is guaranteed to contain.
% Legacy full caches may contain additional fields; they are ignored here
% so that in-scope variables from the grid/medium/source stages are used.
cached = load(filename, 'sensor_data', 'kgrid', 'acoustic_provenance', 'acoustic_info', ...
              'acoustic_grid_dims');
required = {'sensor_data', 'kgrid', 'acoustic_provenance', 'acoustic_info'};
missing = required(~isfield(cached, required));
if ~isempty(missing)
    error('cache_acoustic:corruptCache', ...
        ['Acoustic cache is missing variable(s): %s\n' ...
         'Delete the cache file and re-run the acoustic stage:\n  %s'], ...
        strjoin(missing, ', '), filename);
end
sensor_data         = cached.sensor_data;
kgrid               = cached.kgrid;
acoustic_provenance = cached.acoustic_provenance;
acoustic_info       = cached.acoustic_info;

% Legacy caches may hold acoustic_provenance.freefield_isppa_wcm2 (and its
% copy in tx_calibration.isppa_wcm2) as a serialized scalar gpuArray
% (water_baseline previously cast with double() alone, which does not move
% data off the device). Such a scalar can fail to deserialize, and the first
% touch (isa, gather, isfinite) throws an opaque "Unrecognized field name
% gpuArrayReal". Coerce to host doubles; unreadable values become NaN so the
% caller can substitute the baseline computed in the current run.
if isfield(acoustic_provenance, 'freefield_isppa_wcm2')
    acoustic_provenance.freefield_isppa_wcm2 = ...
        host_scalar(acoustic_provenance.freefield_isppa_wcm2, 'freefield_isppa_wcm2');
end
if isfield(acoustic_provenance, 'tx_calibration') && ...
        isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2')
    acoustic_provenance.tx_calibration.isppa_wcm2 = ...
        host_scalar(acoustic_provenance.tx_calibration.isppa_wcm2, 'tx_calibration.isppa_wcm2');
end

% Grid compatibility check.
% If acoustic_grid_dims was recorded at save time, compare it against the
% current grid.  A mismatch means the cache was produced on a different grid —
% e.g. from a different run, a different subject segmentation, or an
% independent per-transducer run with a different grid extent.
%
% The cache file is NOT deleted: combine_async_intensity can remap fields
% from different grids onto a common reference.  The mismatch is signalled
% via a named error so callers can decide whether to re-run or proceed.
if isfield(cached, 'acoustic_grid_dims') && numel(parameters.grid.dims) == 3
    if ~isequal(cached.acoustic_grid_dims(:)', parameters.grid.dims(:)')
        fprintf(['[cache_acoustic] Acoustic cache grid [%s] does not match current ' ...
                 'grid [%s].\n  Cache retained for potential grid remapping: %s\n'], ...
            num2str(cached.acoustic_grid_dims(:)'), ...
            num2str(parameters.grid.dims(:)'), filename);
        error('cache_acoustic:staleGrid', ...
            ['Acoustic cache was built on a different grid [%s] than the current ' ...
             'grid [%s].\nCache file retained — the caller can re-run the acoustic ' ...
             'simulation or use combine_async_intensity to remap the field.'], ...
            num2str(cached.acoustic_grid_dims(:)'), num2str(parameters.grid.dims(:)'));
    end
end

% Axisymmetric runs: the cached sensor_data is already expanded to 2-D / 3-D,
% while parameters and the grid/medium variables still describe the
% half-plane. ACOUSTIC_STAGE expands those (ACOUSTIC_CONVERT_AXISYMMETRY).

% Fill a missing cached baseline from the current run before scaling, so the
% target can still be met (apply_isppa_scaling errors without a baseline).
acoustic_provenance = fill_baseline_from_live(acoustic_provenance, live_provenance);
sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);

end

%% ------------------------------------------------------------------------
function v = host_scalar(v, name)
try
    v = double(gather(v));
catch ME
    warning('cache_acoustic:unreadableProvenance', ...
        'Cached acoustic_provenance.%s is unreadable (%s); treating it as missing.', ...
        name, ME.message);
    v = NaN;
end
end

%% ------------------------------------------------------------------------
function provenance = fill_baseline_from_live(provenance, live)
% The cached free-field ISPPA describes the drive level of the cached field and
% takes precedence. When it is missing or unreadable, fall back to the baseline
% computed in the current run (acoustic_water_baseline, same source and grid).
cached_ok = isfield(provenance, 'freefield_isppa_wcm2') && ...
    isscalar(provenance.freefield_isppa_wcm2) && isfinite(provenance.freefield_isppa_wcm2);
live_ok = isfield(live, 'freefield_isppa_wcm2') && ...
    isscalar(live.freefield_isppa_wcm2) && isfinite(live.freefield_isppa_wcm2);
if ~cached_ok && live_ok
    fprintf('[acoustic_stage] Using free-water baseline from this run (%.2f W/cm²); cached value unavailable.\n', ...
        live.freefield_isppa_wcm2);
    provenance.freefield_isppa_wcm2 = live.freefield_isppa_wcm2;
    if isfield(live, 'freefield_source')
        provenance.freefield_source = live.freefield_source;
    end
    if isfield(live, 'tx_calibration')
        provenance.tx_calibration = live.tx_calibration;
    end
    if isfield(live, 'freefield_axial')
        provenance.freefield_axial = live.freefield_axial;
    end
end
end
