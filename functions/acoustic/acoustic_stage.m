function [sensor_data, parameters, kgrid, kwave_medium, medium_masks, ...
          segmentation, source, source_labels, acoustic_provenance] = ...
    acoustic_stage(parameters, kgrid, kwave_medium, source, sensor, ...
                   medium_masks, segmentation, source_labels, ...
                   acoustic_provenance, planimg, trans_pos, focus_pos, pseudoCT)
% ACOUSTIC_STAGE  Run or load the acoustic pressure simulation (Stage 6).
%
% Handles the external acoustic NIfTI shortcut, cache load/save, and backend
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
%   sensor_data       - struct with p_max_all field (Pa), or char NIfTI path
%                       for external acoustic input, or [] when unavailable
%   parameters        - updated with state.acoustics_available flag
%   kgrid, kwave_medium, medium_masks, segmentation, source, source_labels -
%                       updated by k-Wave path; unchanged for other backends

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

% --- External acoustic NIfTI (e.g. from BabelBrain or a prior PRESTUS run) ---
% When io.external_acoustic_nifti is set, all acoustic sim stages are skipped
% and the NIfTI path is passed directly to thermal_simulation.
has_external_acoustic = isfield(parameters.io, 'external_acoustic_nifti') && ...
    ~isempty(parameters.io.external_acoustic_nifti);

if has_external_acoustic
    ext_nii = parameters.io.external_acoustic_nifti;
    if ~isfile(ext_nii)
        error('prestus_pipeline:externalAcousticNotFound', ...
            'io.external_acoustic_nifti not found: %s', ext_nii);
    end
    fprintf('[pipeline] Using external acoustic NIfTI: %s\n', ext_nii);
    sensor_data = ext_nii;
    parameters.state.acoustics_available = 1;
    parameters.modules.run_acoustic_analysis = 0;
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
    sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);
    parameters.state.acoustics_available = 1;

elseif exist(filename_sensor_data, 'file')
    disp('Skipping acoustic simulation, loading existing output file.')
    try
        [sensor_data, kgrid, acoustic_provenance, ~, parameters] = ...
            cache_acoustic('load', filename_sensor_data, parameters);
    catch ME
        if strcmp(ME.identifier, 'cache_acoustic:staleGrid')
            % Re-run with the configured backend (not unconditionally k-Wave),
            % scaled exactly like a fresh run.
            fprintf('[pipeline] Acoustic cache grid mismatch — re-running acoustic simulation (%s backend).\n', ...
                acoustic_backend_name(parameters));
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

if use_babelbrain
    sensor_data = babelbrain_wrapper( ...
        parameters, medium_masks, kwave_medium, ...
        focus_pos, trans_pos, planimg, struct(), pseudoCT);
    if isfield(sensor_data, 'babelbrain_version')
        acoustic_provenance.babelbrain_version = sensor_data.babelbrain_version;
        fprintf('[pipeline] BabelBrain %s  (%s)\n', ...
            sensor_data.babelbrain_version.version_txt, ...
            sensor_data.babelbrain_version.git_describe);
    end
    if isfield(sensor_data, 'freefield_isppa_wcm2') && isfinite(sensor_data.freefield_isppa_wcm2)
        acoustic_provenance.freefield_isppa_wcm2             = sensor_data.freefield_isppa_wcm2;
        acoustic_provenance.tx_calibration.isppa_wcm2        = sensor_data.freefield_isppa_wcm2;
        acoustic_provenance.tx_calibration.isppa_source      = 'babelbrain_water_sim';
    end
    if should_save_output(parameters.io, 'save_acoustic_matrices')
        cache_acoustic('save', filename_sensor_data, parameters, sensor_data, [], ...
            acoustic_provenance, struct('parameters', parameters, 'backend', 'babelbrain'), kwave_medium);
        fprintf('[pipeline] babelbrain sensor_data cached: %s\n', filename_sensor_data);
    end

elseif use_jwave
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric
        error('[pipeline] jwave backend does not support axisymmetric grids.');
    end
    if isfield(parameters.grid, 'use_kWaveArray') && parameters.grid.use_kWaveArray
        error('[pipeline] jwave backend requires direct source geometry (grid.use_kWaveArray: 0).');
    end
    sensor_data = jwave_wrapper(parameters, kgrid, kwave_medium, source, source_labels, sensor);
    if should_save_output(parameters.io, 'save_acoustic_matrices')
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

end

function name = acoustic_backend_name(parameters)
% ACOUSTIC_BACKEND_NAME  Configured acoustic backend, defaulting to 'kwave'.
    name = 'kwave';
    if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
        name = char(parameters.acoustic_backend);
    end
end
