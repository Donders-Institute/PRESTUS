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
%   Applies the axisymmetric coordinate fix-up when needed, and calls
%   apply_isppa_scaling so the caller receives a correctly scaled field.
%
% Use as:
%   % save
%   cache_acoustic('save', filename, parameters, sensor_data, kgrid, ...
%                  acoustic_provenance, acoustic_info, kwave_medium)
%
%   % load
%   [sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
%       cache_acoustic('load', filename, parameters)
%
% Input:
%   mode              - 'save' or 'load'
%   filename          - full path to the cache .mat file
%   parameters        - PRESTUS config struct (updated on load for axisym)
%   sensor_data       - (save only) struct with p_max_all [Pa]
%   kgrid             - (save only) kWaveGrid object
%   acoustic_provenance - (save only) struct with freefield_isppa_wcm2
%   acoustic_info     - (save only) struct with parameters snapshot
%   kwave_medium      - (save only) medium properties; included in async mode
%
% Output (load mode only):
%   sensor_data       - loaded and ISPPA-scaled pressure struct
%   kgrid             - kWaveGrid object
%   acoustic_provenance - provenance struct
%   acoustic_info     - info struct
%   parameters        - updated parameters (axisym coord fix applied)

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
    acoustic_cache_load(filename, parameters);

end

% =========================================================================

function acoustic_cache_save(filename, parameters, sensor_data, kgrid, ...
                              acoustic_provenance, acoustic_info, kwave_medium)

if is_async_mode(parameters)
    save(filename, ...
        'sensor_data', ...
        'kgrid', ...
        'kwave_medium', ...
        'acoustic_info', ...
        'acoustic_provenance', '-v7.3')
else
    save(filename, ...
        'sensor_data', ...
        'kgrid', ...
        'acoustic_info', ...
        'acoustic_provenance', '-v7.3')
end

end

% =========================================================================

function [sensor_data, kgrid, acoustic_provenance, acoustic_info, parameters] = ...
        acoustic_cache_load(filename, parameters)

% Load only the slim variables the cache is guaranteed to contain.
% Legacy full caches may contain additional fields; they are ignored here
% so that in-scope variables from the grid/medium/source stages are used.
cached = load(filename, 'sensor_data', 'kgrid', 'acoustic_provenance', 'acoustic_info');
sensor_data         = cached.sensor_data;
kgrid               = cached.kgrid;
acoustic_provenance = cached.acoustic_provenance;
acoustic_info       = cached.acoustic_info;

% Axisymmetric coordinate fix-up.
% acoustic_wrapper saves sensor_data after expanding the axisymmetric
% result to 2-D / 3-D, but parameters still holds the pre-expansion dims
% and transducer positions when we bypass acoustic_wrapper via cache.
if numel(parameters.grid.dims) == 2 && ...
        isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric
    expand_to_3d = strcmp(parameters.simulation.medium, 'phantom') || ...
                   (isfield(parameters.modules, 'run_heating_sims') && ...
                    parameters.modules.run_heating_sims == 1);
    if ~expand_to_3d
        % The cached sensor_data was already expanded to 2D bilateral by
        % convert_axisymmetric_to_2d, which restores the original dims and
        % positions.  Read those from the bilateral snapshot stored by
        % grid_axisymmetry so the loaded parameters match.
        if isfield(parameters.grid, 'axisym_bilateral_dims')
            bilateral_dims = parameters.grid.axisym_bilateral_dims;
            for ti = 1:numel(parameters.transducer)
                parameters.transducer(ti).trans_pos = ...
                    fliplr(parameters.transducer(ti).trans_pos_bilateral);
                parameters.transducer(ti).focus_pos = ...
                    fliplr(parameters.transducer(ti).focus_pos_bilateral);
            end
            parameters.grid.dims = fliplr(bilateral_dims);
        end
    end
end

sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);

end
