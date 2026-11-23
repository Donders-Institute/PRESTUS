function [results_heating, time_status_seq, kwave_medium, medium_masks, ...
          segmentation, planimg, highlighted_pos, parameters] = ...
    thermal_stage(parameters, sensor_data, kgrid, kwave_medium, medium_plus, ...
                  medium_masks, segmentation, planimg, source, sensor, ...
                  trans_pos, focus_pos, acoustic_provenance, fov_state, highlighted_pos)
% THERMAL_STAGE  Run or load the bio-heat thermal simulation.
%
% Handles acoustic-FOV restore, thermal-grid setup, backend dispatch
% (k-Wave default / BabelBrain), and cache save/load.  Returns results
% together with the (possibly restored) medium_masks, segmentation, planimg,
% and highlighted_pos so that downstream thermal_analysis and nifti_thermal
% operate in a consistent coordinate space.
%
% FOV restore logic (controlled by fov_state):
%   fov_state.active == false              — no FOV was applied; pass through
%   fov_state.active && ~thermal_use_fov   — restore full-grid state before sim
%   fov_state.active &&  thermal_use_fov   — stay in FOV coords; adjust planimg
%
% Use as:
%   [results_heating, time_status_seq, kwave_medium, medium_masks, ...
%    segmentation, planimg, highlighted_pos, parameters] = ...
%       thermal_stage(parameters, sensor_data, kgrid, kwave_medium, medium_plus, ...
%                     medium_masks, segmentation, planimg, source, sensor, ...
%                     trans_pos, focus_pos, acoustic_provenance, fov_state)
%   [...] = thermal_stage(..., highlighted_pos)
%
% Inputs:
%   fov_state       - struct from acoustic_fov_crop (fields: active, offset,
%                     full_dims, kwave_medium_full, medium_masks_full,
%                     segmentation_full, trans_pos_full, focus_pos_full)
%   highlighted_pos - focus voxel in current grid coords ([] if not available)
%   All other inputs mirror the pipeline workspace after acoustic stage.
%
% Outputs:
%   results_heating  - struct with temperature and CEM43 timeseries, on the
%                      acoustic grid (a separate thermal grid's results are
%                      mapped back: grid.thermal_resolution_mm / thermal_fov_mm)
%   time_status_seq  - protocol on/off event sequence for thermal_plot_sim
%   kwave_medium     - input kwave_medium updated with temp_0 from thermal grid
%   medium_masks     - tissue masks in the coordinate space used by thermal sim
%   segmentation     - segmentation in the coordinate space used by thermal sim
%   planimg          - possibly updated (adjusted affine for thermal_use_acoustic_fov)
%   highlighted_pos  - translated to sim-grid coords if FOV was restored
%   parameters       - updated with state.heating_available

arguments
    parameters        (1,1) struct
    sensor_data
    kgrid
    kwave_medium      (1,1) struct
    medium_plus
    medium_masks
    segmentation
    planimg           (1,1) struct
    source
    sensor
    trans_pos
    focus_pos
    acoustic_provenance (1,1) struct
    fov_state         (1,1) struct
    highlighted_pos   = []
end

parameters.state.heating_available = 0;
results_heating  = struct();
time_status_seq  = [];

if ~(isfield(parameters.modules, 'run_heating_sims') && ...
        parameters.modules.run_heating_sims && ...
        parameters.state.acoustics_available == 1)
    return;
end

if isfield(parameters.io, 'thermal_cache_affix')
    thermal_file_affix = parameters.io.thermal_cache_affix;
else
    thermal_file_affix = parameters.io.output_affix;
end
filename_heating_data = fullfile(parameters.io.dir_cache, ...
    sprintf('sub-%03d_%s_heating_res%s.mat', ...
    parameters.subject_id, parameters.simulation.medium, thermal_file_affix));

run_new_sim = confirm_overwriting(filename_heating_data, parameters) && ...
    (parameters.simulation.interactive == 0 || ...
     confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No'));

parameters.state.recomputed.thermal = logical(run_new_sim);   % OUTPUTS_FOLLOW_RECOMPUTE

if run_new_sim

    thermal_use_acoustic_fov = fov_state.active && ...
        isfield(parameters.grid, 'thermal_use_acoustic_fov') && ...
        parameters.grid.thermal_use_acoustic_fov;

    if fov_state.active && ~thermal_use_acoustic_fov
        % Restore full-grid state — acoustic FOV scope ends here.
        if ~isempty(sensor_data) && isstruct(sensor_data)
            if isfield(sensor_data, 'p_max_all')
                sensor_data.p_max_all = fov_backproject_volume( ...
                    single(sensor_data.p_max_all), fov_state.offset, fov_state.full_dims);
            end
            for fn = {'p_max_async', 'p_min_async', 'p_min_all'}
                if isfield(sensor_data, fn{1})
                    sensor_data.(fn{1}) = fov_backproject_volume( ...
                        single(sensor_data.(fn{1})), fov_state.offset, fov_state.full_dims);
                end
            end
        end
        kwave_medium  = fov_state.kwave_medium_full;
        medium_masks  = fov_state.medium_masks_full;
        segmentation  = fov_state.segmentation_full;
        trans_pos     = fov_state.trans_pos_full;
        focus_pos     = fov_state.focus_pos_full;
        parameters.grid.dims = fov_state.full_dims;
        parameters = fov_restore_transducers(parameters, fov_state);
        if ~isempty(highlighted_pos)
            highlighted_pos = fov_to_grid_pos(highlighted_pos, fov_state.offset);
        end

    elseif thermal_use_acoustic_fov
        % Stay in acoustic FOV coords for thermal sim; adjust planimg affines.
        M_ac2fov        = eye(4);
        M_ac2fov(4,1:3) = -(fov_state.offset - 1);
        planimg.transf  = planimg.transf * M_ac2fov;
        if isfield(planimg, 'inv_transf')
            M_fov2ac            = eye(4);
            M_fov2ac(4,1:3)     = (fov_state.offset - 1);
            planimg.inv_transf  = M_fov2ac * planimg.inv_transf;
        end
        fprintf('[pipeline] Thermal FOV = acoustic FOV (offset [%s], dims [%s]).\n', ...
            num2str(fov_state.offset), num2str(parameters.grid.dims));
    end

    % Axisymmetric runs: the acoustic stage expanded the positions in
    % parameters; the pipeline's trans_pos/focus_pos are the half-plane ones.
    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1 && ...
            isfield(parameters.transducer(1), 'trans_pos_bilateral') && ...
            numel(parameters.transducer(1).trans_pos) == numel(parameters.grid.dims)
        trans_pos = parameters.transducer(1).trans_pos;
        focus_pos = parameters.transducer(1).focus_pos;
    end

    kwave_medium.temp_0              = medium_plus.temp_0;
    % Heat source from absorption_np [Np/m]: fraction x attenuation formed per
    % tissue and smoothed as one map (medium_setup). Multiplying the separately
    % smoothed attenuation and fraction maps overshoots both tissues where one
    % rises and the other falls (skin -> bone: up to 8x in the last skin voxel
    % at 1 mm FWHM). The fraction path is kept for medium caches without
    % absorption_np.
    if isfield(medium_plus, 'absorption_np')
        kwave_medium.absorption_np = medium_plus.absorption_np;
    else
        kwave_medium.absorption_fraction = medium_plus.absorption_fraction;
    end
    % async transducers at different frequencies: absorption at each frequency
    if isfield(medium_plus, 'absorption_power') && numel(unique([parameters.transducer.freq_hz])) > 1
        kwave_medium.absorption_power = medium_plus.absorption_power;
    end
    % medium_plus holds full-grid maps; with thermal_use_acoustic_fov the medium
    % stays on the acoustic FOV, so crop them to it
    if thermal_use_acoustic_fov
        for fn = {'temp_0', 'absorption_np', 'absorption_fraction', 'absorption_power'}
            if isfield(kwave_medium, fn{1})
                kwave_medium.(fn{1}) = crop_to_fov(kwave_medium.(fn{1}), fov_state.offset, ...
                    parameters.grid.dims, fov_state.full_dims);
            end
        end
    end

    use_babelbrain_acoustic = isfield(parameters, 'acoustic_backend') && ...
        strcmpi(parameters.acoustic_backend, 'babelbrain');
    has_external_acoustic   = ischar(sensor_data);  % external NIfTI path
    use_separate_grid = has_external_acoustic || fov_state.active || use_babelbrain_acoustic || ...
        (isfield(parameters.grid, 'thermal_resolution_mm') && ...
         ~isempty(parameters.grid.thermal_resolution_mm) && ...
         parameters.grid.thermal_resolution_mm ~= parameters.grid.resolution_mm) || ...
        (isfield(parameters.grid, 'thermal_fov_mm') && ~isempty(parameters.grid.thermal_fov_mm));

    use_bb_thermal = isfield(parameters, 'thermal_backend') && ...
        strcmpi(parameters.thermal_backend, 'babelbrain');

    if use_separate_grid
        % k-Wave heating needs the pressure on the thermal grid: it gets the
        % medium's crop and resampling, as p^2 so that each thermal voxel
        % holds the mean intensity (heat source) over the voxels it covers.
        % BabelBrain thermal reads the acoustic-grid pressure itself.
        medium_in  = kwave_medium;
        resample_p = ~use_bb_thermal && isstruct(sensor_data) && isfield(sensor_data, 'p_max_all');
        if resample_p
            medium_in.p_sq_for_thermal = single(gather(sensor_data.p_max_all)).^2;
        end
        [kgrid_th, kwave_medium_th, medium_masks_th, parameters_th, transf_th] = ...
            thermal_grid_setup(parameters, planimg, medium_in, medium_masks, trans_pos, focus_pos);
        clear medium_in
        if resample_p
            if isfield(parameters_th.grid, 'thermal_fov_start_ac')   % grids differ
                sensor_data.p_max_all = sqrt(max(kwave_medium_th.p_sq_for_thermal, 0));
            end
            kwave_medium_th = rmfield(kwave_medium_th, 'p_sq_for_thermal');
        end
        sensor_th      = struct('mask', ones(parameters_th.grid.dims));
        source_th      = struct('Q',    zeros(parameters_th.grid.dims));
    else
        kgrid_th        = kgrid;
        kwave_medium_th = kwave_medium;
        medium_masks_th = medium_masks;
        parameters_th   = parameters;
        is_layered      = contains(parameters.simulation.medium, 'layered');
        assert(~is_layered || (isfield(planimg, 'transf') && ~isempty(planimg.transf)), ...
            ['thermal_stage:missingPlanimg — planimg.transf required for thermal NIfTI export. ' ...
             'Re-run with run_grid_setup=1 to regenerate grid cache and planimg sidecar.']);
        transf_th = planimg.transf;
        if ~isempty(source) && isstruct(source) && ~isempty(sensor) && isstruct(sensor)
            sensor_th = sensor;
            source_th = source;
        else
            sensor_th = struct('mask', ones(parameters_th.grid.dims));
            source_th = struct('Q',    zeros(parameters_th.grid.dims));
        end
    end

    if use_bb_thermal
        % Propagate free-water ISPPA for the wrapper's transmission log line.
        % The Isppa actually passed to BabelBrain is the in-brain peak of
        % sensor_data (see babelbrain_thermal_wrapper), so pass the medium and
        % masks on the ACOUSTIC grid (same grid as sensor_data), not *_th.
        if ~isfield(parameters_th.transducer(1), 'tx_calibration') || ...
                ~isfield(parameters_th.transducer(1).tx_calibration, 'isppa_wcm2') || ...
                ~isfinite(parameters_th.transducer(1).tx_calibration.isppa_wcm2)
            if isfield(acoustic_provenance, 'freefield_isppa_wcm2') && ...
                    isfinite(acoustic_provenance.freefield_isppa_wcm2)
                parameters_th.transducer(1).tx_calibration.isppa_wcm2   = acoustic_provenance.freefield_isppa_wcm2;
                parameters_th.transducer(1).tx_calibration.isppa_source = 'acoustic_provenance';
            elseif isfield(acoustic_provenance, 'tx_calibration') && ...
                    isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2') && ...
                    isfinite(acoustic_provenance.tx_calibration.isppa_wcm2)
                parameters_th.transducer(1).tx_calibration = acoustic_provenance.tx_calibration;
            end
        end
        results_heating = babelbrain_thermal_wrapper(parameters_th, sensor_data, planimg, ...
            kwave_medium, medium_masks);
        kwaveDiffusion  = [];
        % Build k-Wave-compatible time_status_seq from protocol timing so
        % thermal_plot_sim can shade on/off periods on the temperature curve.
        if isfield(results_heating, 'timeseries') && ...
                isfield(results_heating.timeseries, 'TimeProfileTarget') && ...
                ~isempty(results_heating.timeseries.TimeProfileTarget)
            % One entry per BabelBrain time sample (the time axis of its
            % per-tissue series), 'on' within a pulse train of the protocol
            % BabelBrain ran (BABELBRAIN_THERMAL_TIMING)
            t_ax_bb = results_heating.timeseries.TimeProfileTarget(:)';
            bt      = babelbrain_thermal_timing(parameters_th.timing);
            period  = bt.duration_s + bt.duration_off_s;
            k_rep   = floor(t_ax_bb / max(period, eps));
            on      = k_rep < bt.repetitions & (t_ax_bb - k_rep * period) < bt.duration_s;
            status  = repmat({'off'}, size(t_ax_bb));
            status(on) = {'on'};
            time_status_seq = struct('time', num2cell(t_ax_bb), 'recorded', 1, 'status', status);
        else
            time_status_seq = struct([]);
        end
    else
        t_thermal = tic;
        [kwaveDiffusion, time_status_seq, results_heating] = ...
            thermal_simulation(parameters_th, sensor_data, kgrid_th, kwave_medium_th, ...
                               sensor_th, source_th, transf_th, medium_masks_th);
        % Size and duration for resource estimates (run_stats.jsonl)
        run_stats('add', 'thermal_s', toc(t_thermal), ...
            'thermal_steps', time_status_seq(end).step);
        run_stats('set', 'thermal_voxels', kgrid_th.total_grid_points);
    end

    % A separate thermal grid: the results go back onto the acoustic grid, on
    % which the cache, thermal_analysis and nifti_thermal (masks, planimg) work
    if isfield(parameters_th.grid, 'thermal_fov_start_ac')
        results_heating = thermal_results_to_acoustic(results_heating, parameters, ...
            parameters_th, kwave_medium.temp_0);
    end

    if should_save_output(parameters.io, 'save_thermal_matrices')
        save(filename_heating_data, ...
            'kwaveDiffusion', 'time_status_seq', 'sensor_th', ...
            'results_heating', 'kwave_medium_th', '-v7.3');
    else
        disp('Not saving thermal simulation output matrices ...')
    end

    if isfield(kwave_medium_th, 'temp_0') && ~isfield(parameters_th.grid, 'thermal_fov_start_ac')
        kwave_medium.temp_0 = kwave_medium_th.temp_0;
    end

    parameters.state.heating_available = 1;

elseif exist(filename_heating_data, 'file')
    disp('Skipping thermal simulation, loading existing output file.')
    time_status_seq = [];   % for caches without time_status_seq
    loaded = load(filename_heating_data);
    results_heating = loaded.results_heating;
    if isfield(loaded, 'time_status_seq')
        time_status_seq = loaded.time_status_seq;
    end
    if isfield(loaded, 'kwave_medium_th') && isfield(loaded.kwave_medium_th, 'temp_0')
        kwave_medium.temp_0 = loaded.kwave_medium_th.temp_0;
        % a separate thermal grid's temp_0 is not on the acoustic grid
        if ~isequal(size(kwave_medium.temp_0), double(parameters.grid.dims)) && ...
                isstruct(medium_plus) && isfield(medium_plus, 'temp_0') && ...
                isequal(size(medium_plus.temp_0), double(parameters.grid.dims))
            kwave_medium.temp_0 = medium_plus.temp_0;
        end
    end
    parameters.state.heating_available = 1;

else
    % Reached only when an interactive user declines the thermal simulation
    warning('prestus_pipeline:thermalNoAcoustics', ...
        'Thermal simulation not run and no heating cache found (%s); skipping thermal analysis.', ...
        filename_heating_data);
    parameters.state.heating_available  = 0;
    parameters.modules.run_thermal_analysis = 0;
end

end

% =========================================================================
function res = thermal_results_to_acoustic(res, parameters, parameters_th, temp_0_ac)
% Map thermal-grid result volumes (3-D maps and [x z t] focal planes) onto the
% acoustic grid. Nearest-neighbour, so every thermal voxel value (and so each
% peak) is kept. Outside the thermal FOV: the baseline temperature temp_0_ac
% for temperatures, 0 for everything else (CEM43, pressure). The maps take the
% class of temp_0_ac: the thermal grid works in single (thermal_grid_setup),
% and results should not change class with the thermal grid option.
lo      = parameters_th.grid.thermal_fov_start_ac;
hi      = parameters_th.grid.thermal_fov_end_ac;
box     = hi - lo + 1;
ac_dims = double(parameters.grid.dims);
th_dims = double(parameters_th.grid.dims);
y_ac    = parameters.transducer(1).trans_pos(2);   % focal plane, as thermal_analysis reads it
temp_0_ac = gather(temp_0_ac);
if isscalar(temp_0_ac)
    temp_0_ac = repmat(temp_0_ac, ac_dims);
end
temp_fields = {'maxT', 'heating_endT', 'focal_planeT'};
fns = fieldnames(res);
for i = 1:numel(fns)
    f = fns{i};
    v = res.(f);
    if ~isnumeric(v) || isempty(v)
        continue
    end
    v = cast(gather(v), 'like', temp_0_ac);
    is_temp = ismember(f, temp_fields);
    if startsWith(f, 'focal_plane')
        if ~isequal(size(v, [1 2]), th_dims([1 3]))
            continue
        end
        nt = size(v, 3);
        if is_temp
            out = repmat(cast(squeeze(temp_0_ac(:, y_ac, :)), 'like', v), [1 1 nt]);
        else
            out = zeros([ac_dims([1 3]) nt], 'like', v);
        end
        out(lo(1):hi(1), lo(3):hi(3), :) = imresize3(v, [box([1 3]) nt], 'nearest');
    elseif isequal(size(v), th_dims)
        if is_temp
            out = cast(temp_0_ac, 'like', v);
        else
            out = zeros(ac_dims, 'like', v);
        end
        out(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = imresize3(v, box, 'nearest');
    else
        continue
    end
    res.(f) = out;
end
end

function v = crop_to_fov(v, offset, dims, full_dims)
% Crop a full-grid map (first three dimensions) to the acoustic FOV; scalars and
% maps already on the FOV pass through.
sz = size(v);
if numel(sz) < 3 || ~isequal(sz(1:3), double(full_dims(:)'))
    return
end
e = offset + dims - 1;
v = v(offset(1):e(1), offset(2):e(2), offset(3):e(3), :);
end
