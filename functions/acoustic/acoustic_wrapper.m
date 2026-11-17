function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    acoustic_wrapper(...
    parameters, ...
    kgrid, ...
    kwave_medium, ...
    source, ...
    sensor, ...
    medium_masks, ...
    filename_sensor_data, ...
    segmentation, ...
    source_labels, ...
    acoustic_provenance)
% ACOUSTIC_WRAPPER  Run the acoustic simulation and convert axisymmetric output to 2D or 3D
%
% Calls acoustic_simulation with PMLInside=false and an optional skull-edge
% display mask. After the k-Wave run, detects axisymmetric grid output and
% expands it: phantom simulations and runs with heating always expand to 3D
% (radial expansion via radialExpand2DTo3D); acoustic-only non-phantom runs
% expand to a bilateral 2D cross-section instead. Complex fields
% (io.save_p_complex) go to their own files (WRITE_P_COMPLEX_FILES); the
% simulation matrices are cached with CACHE_ACOUSTIC (io.save_acoustic_matrices).
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, kwave_medium, ...
%    kgrid, source, source_labels] = ...
%       acoustic_wrapper(parameters, kgrid, kwave_medium, source, sensor, ...
%                        medium_masks, filename_sensor_data, segmentation, source_labels)
%       acoustic_wrapper(..., acoustic_provenance)
%
% Input:
%   parameters          - (1,1) struct, PRESTUS config
%   kgrid               - (1,1) kWaveGrid, simulation grid
%   kwave_medium        - (1,1) struct, spatial medium property maps
%   source              - (1,1) struct, k-Wave source definition
%   sensor              - (1,1) struct, k-Wave sensor definition
%   medium_masks        - (:,:) or (:,:,:) numeric, tissue layer label map
%   filename_sensor_data- (1,:) char, path for optional .mat output file
%   segmentation        - (:,:) or (:,:,:) numeric, CHARM tissue label map
%   source_labels       - (:,:) or (:,:,:) numeric, transducer element label map
%   acoustic_provenance - (1,1) struct, optional; FOV crop offsets (fov_offset_ac,
%                         full_ac_dims) for the p_complex files, saved in the cache
%
% Output:
%   sensor_data  - (1,1) struct, k-Wave output expanded to 2D or 3D [Pa]
%   parameters   - (1,1) struct, updated grid.dims and transducer positions
%   segmentation - (:,:,:) or (:,:) numeric, expanded segmentation volume
%   medium_masks - (:,:,:) or (:,:) numeric, expanded medium mask volume
%   kwave_medium - (1,1) struct, expanded medium property maps
%   kgrid        - (1,1) kWaveGrid, updated 2D or 3D grid
%   source       - (1,1) struct, expanded source with p_mask
%   source_labels- (:,:,:) or (:,:) numeric, expanded element label map
%
% See also: ACOUSTIC_SIMULATION, ACOUSTIC_ANALYSIS, CONVERT_AXISYMMETRIC_TO_3D,
%           CONVERT_AXISYMMETRIC_TO_2D

arguments
    parameters           (1,1) struct
    kgrid                (1,1)
    kwave_medium         (1,1) struct
    source               (1,1) struct
    sensor               (1,1) struct
    medium_masks         {mustBeNumericOrLogical}
    filename_sensor_data (1,:) char
    segmentation         {mustBeNumericOrLogical}
    source_labels        {mustBeNumericOrLogical}
    acoustic_provenance  (1,1) struct = struct()
end

    % PML is always applied outside the user-defined grid (PMLInside=false).
    % The returned sensor data (p_max_all, p_final) is automatically cropped
    % by k-Wave to the original grid dimensions; callers never see the PML voxels.
    % grid.pml_size may be a positive integer or the string 'auto' (k-Wave selects
    % the optimal size for FFT efficiency; only valid with PMLInside=false).
    %
    % Resolve 'auto' to the actual size via getOptimalPMLSize so the effective
    % value is available for reporting and reproducibility.
    [effective_pml, pml_alpha] = kwave_pml_settings(parameters, kgrid);
    parameters.grid.pml_size_effective = effective_pml;

    % --- Simulation summary (mirrors jwave_wrapper format) ---
    grid_dims   = parameters.grid.dims;
    n_voxels    = prod(grid_dims);
    n_src       = nnz(source.p_mask);
    freq_kHz    = parameters.transducer(1).freq_hz / 1e3;
    precision   = char(parameters.simulation.precision);
    code_type   = parameters.simulation.code_type;

    pml_scalar  = max(effective_pml(:));
    if isscalar(effective_pml) || numel(unique(effective_pml(:))) == 1
        pml_size_str = num2str(pml_scalar);
    else
        pml_size_str = sprintf('[%s] → max=%d', num2str(effective_pml(:)'), pml_scalar);
    end
    % k-Wave PML is outside the grid: actual allocated dims are grid + 2*PML per axis
    pml_total_dims = grid_dims(:)' + 2 * pml_scalar;

    c  = kwave_medium.sound_speed;
    rh = kwave_medium.density;
    is_elastic = isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'kwave_elastic');
    fprintf('[acoustic_wrapper] ---- k-Wave simulation summary ----\n');
    if is_elastic
        fprintf('[acoustic_wrapper]   backend   : %s  (elastic: pstdElastic, shear waves in skull)\n', code_type);
    else
        fprintf('[acoustic_wrapper]   backend   : %s\n', code_type);
    end
    fprintf('[acoustic_wrapper]   precision : %s\n', precision);
    fprintf('[acoustic_wrapper]   grid      : [%s] = %dk voxels\n', num2str(grid_dims), round(n_voxels/1e3));
    fprintf('[acoustic_wrapper]   grid res  : %.2f mm\n', parameters.grid.resolution_mm);
    fprintf('[acoustic_wrapper]   time      : Nt=%d  dt=%.3e s\n', kgrid.Nt, kgrid.dt);
    fprintf('[acoustic_wrapper]   source    : %d active voxels  freq=%.1f kHz\n', n_src, freq_kHz);
    fprintf('[acoustic_wrapper]   PML size  : %s voxels/side  location=outside-grid  total=[%s]\n', ...
        pml_size_str, num2str(pml_total_dims));
    if isscalar(c)
        fprintf('[acoustic_wrapper]   medium    : homogeneous  c=%.0f m/s  rho=%.0f kg/m3\n', c, rh);
    else
        fprintf('[acoustic_wrapper]   medium    : heterogeneous  c=[%.0f–%.0f] m/s  rho=[%.0f–%.0f] kg/m3\n', ...
            min(c(:)), max(c(:)), min(rh(:)), max(rh(:)));
    end
    fprintf('[acoustic_wrapper] ------------------------------------\n');

    % --- GPU / RAM memory estimate ---
    % Estimate the k-Wave computational grid (grid + 2×PML per axis, before
    % FFT-optimal rounding). For a 3D heterogeneous absorbing layered medium,
    % k-Wave allocates roughly 40 single-precision arrays of this size.
    if ismember(code_type, {'matlab_gpu', 'cpp_gpu'})
        bytes_per_elem = 4 * (1 + strcmp(precision, 'double'));
        if isscalar(effective_pml)
            pml_per_ax = repmat(double(effective_pml), 1, numel(grid_dims));
        else
            pml_per_ax = double(effective_pml(:)');
        end
        comp_dims_est = grid_dims(:)' + 2 * pml_per_ax;
        comp_voxels   = prod(comp_dims_est);
        n_kwave_arrays = 40;
        if is_elastic
            % 9 split stress + 9 split velocity fields, Lamé/viscosity maps,
            % staggered-grid copies and PML operators (rough upper bound)
            n_kwave_arrays = 70;
        end
        gpu_mem_gb  = n_kwave_arrays * comp_voxels * bytes_per_elem / 1e9;
        fprintf('[acoustic_wrapper]   est. comp grid : [%s] voxels (before FFT rounding)\n', num2str(comp_dims_est));
        fprintf('[acoustic_wrapper]   est. GPU mem   : ~%.1f GB  (%d arrays × %.0fM vox × %d B)\n', ...
            gpu_mem_gb, n_kwave_arrays, comp_voxels/1e6, bytes_per_elem);
        if strcmp(code_type, 'matlab_gpu')
            % matlab_gpu: k-Wave casts every medium/source array from CPU to
            % gpuArray inside kspaceFirstOrder_dataCast. During each cast,
            % the CPU original and the GPU copy coexist, so MATLAB RAM spikes
            % by approximately the GPU footprint on top of what is already
            % allocated (source weights, medium maps, etc.).
            fprintf('[acoustic_wrapper]   ⚠  matlab_gpu DataCast: CPU+GPU copies coexist (~%.1f GB extra RAM)\n', gpu_mem_gb);
            fprintf('[acoustic_wrapper]      Recommended job memory ≥ current_workspace + %.0f GB.\n', ceil(gpu_mem_gb));
            fprintf('[acoustic_wrapper]   ℹ  cpp_gpu avoids this spike: medium/source arrays are\n');
            fprintf('[acoustic_wrapper]      serialised to HDF5 on disk before the C++ binary\n');
            fprintf('[acoustic_wrapper]      allocates GPU memory in its own process, so MATLAB''s\n');
            fprintf('[acoustic_wrapper]      heap does not hold CPU+GPU copies simultaneously.\n');
        end
        fprintf('[acoustic_wrapper] ------------------------------------\n');
    end

    % Pass the resolved integer vector so k-Wave uses exactly the same PML
    % that is recorded in pml_size_effective (passing 'auto' would let k-Wave
    % resolve it independently, risking a mismatch with our stored value).
    kwave_input_args = struct('PMLInside', false, ...
        'PMLSize', effective_pml, ...
        'PlotPML', true);
    if ~isempty(pml_alpha)
        kwave_input_args.PMLAlpha = pml_alpha;   % else k-Wave's 2
    end

    if contains(parameters.simulation.medium, {'layered'}) && ...
            any(ismember(fieldnames(parameters.layers), {'skull'}))
        % Extract the skull edge ...
        mask = tissuemask_binary(parameters, medium_masks);
        skull_log  = logical(mask.skull);
        skull_edge = skull_log & (convn(single(~skull_log), ones(3,3,3,'single'), 'same') > 0);
        % ... to set as display mask
        kwave_input_args.DisplayMask = skull_edge;
    end

    if parameters.modules.run_source_setup==0
        error('Source setup not requested. Not able to proceed with acoustic simulation.')
    end

    % Complex pressure field only when explicitly requested: the solver records
    % the steady-state time series at every grid point (ACOUSTIC_SIMULATION).
    save_p_complex = isfield(parameters.io, 'save_p_complex') && parameters.io.save_p_complex;
    p_complex_freqs = [];
    if save_p_complex
        p_complex_freqs = unique([parameters.transducer.freq_hz]);
    end

    % perform the acoustic simulation
    % kwave_elastic: add shear speed/absorption in the skull. The elastic medium
    % is a separate struct so the returned kwave_medium stays the fluid medium
    % consumed by thermal simulation and reporting.
    use_elastic = isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'kwave_elastic');
    if use_elastic
        if contains(parameters.simulation.medium, {'layered'})
            sim_medium = medium_elastic_properties(parameters, kwave_medium, medium_masks);
        else
            sim_medium = medium_elastic_properties(parameters, kwave_medium, []);
        end
        save_elastic_debug_maps(parameters, sim_medium);
    else
        sim_medium = kwave_medium;
    end
    sensor_data = acoustic_simulation(kgrid, sim_medium, source, sensor, kwave_input_args, parameters, p_complex_freqs);
    clear sim_medium

    % Complex pressure amplitude over the recorded steady-state window
    % (STEADY_STATE_COMPLEX_AMPLITUDE, via ACOUSTIC_SIMULATION).
    % p_complex = amp * exp(1i * phase) preserves interference information for
    % coherent multi-transducer superposition. A single-frequency setup produces
    % sensor_data.p_complex; a multi-frequency setup produces one field per unique
    % frequency, keyed as p_complex_<freq_Hz>Hz. The phase is referred to t = 0
    % of the simulation, where every run's source signal starts, so fields from
    % separate runs (e.g. one per transducer, possibly on FOV sub-grids with a
    % different dt) can be added coherently (SUPERPOSE_COMPLEX_FIELDS).
    if isfield(sensor_data, 'p_complex_fit')
        unique_freqs = p_complex_freqs;
        p_window = sensor_data.p_complex_window;   % set by ACOUSTIC_SIMULATION
        sensor_data = rmfield(sensor_data, 'p_complex_window');
        t0_s = (double(p_window.record_start_index) - 1) * kgrid.dt;
        p_complex_keys = cell(1, numel(unique_freqs));
        if numel(unique_freqs) > 1
            fprintf('Multi-frequency setup: extracting p_complex at %d frequencies (%s Hz).\n', ...
                numel(unique_freqs), num2str(unique_freqs));
        end
        for fi = 1:numel(unique_freqs)
            % one value per sensor point; the sensor mask is the whole grid, so
            % the points are in linear-index order
            pc = cast(sensor_data.p_complex_fit(:, fi), char(parameters.simulation.precision));
            if numel(pc) == prod(parameters.grid.dims)
                pc = reshape(pc, parameters.grid.dims);
            end
            if numel(unique_freqs) == 1
                p_complex_keys{fi} = 'p_complex';
            else
                p_complex_keys{fi} = sprintf('p_complex_%gHz', unique_freqs(fi));
            end
            sensor_data.(p_complex_keys{fi}) = pc;
        end
        clear pc
        sensor_data = rmfield(sensor_data, 'p_complex_fit');
        % Saved with each p_complex file (one entry per frequency below).
        p_complex_meta = struct( ...
            'phase_reference', 'simulation_t0', ...
            'freq_hz',         unique_freqs, ...
            'field_names',     {p_complex_keys}, ...
            'record_t0_s',     t0_s, ...
            'dt_s',            kgrid.dt, ...
            'method',          p_window.method, ...           % 'quarter_period' | 'steady_state_window'
            'n_samples',       p_window.n_samples, ...        % recorded time steps used by the fit
            'p_max_from_fit',  p_window.p_max_from_fit, ...   % p_max_all / p_min_all are +/- |p_complex|
            'drive_level',     'simulated', ...   % before any ISPPA rescaling
            'grid_dims',       [], ...
            'fov_offset_ac',   [], ...
            'full_ac_dims',    []);
        if nargin >= 10 && isstruct(acoustic_provenance) && isfield(acoustic_provenance, 'fov_offset_ac')
            p_complex_meta.fov_offset_ac = acoustic_provenance.fov_offset_ac;
            p_complex_meta.full_ac_dims  = acoustic_provenance.full_ac_dims;
        end
    end

    % convert media and results to 2D/3D (if axisymmetry was used); the form
    % (AXISYMMETRIC_OUTPUT_MODE) is shared with the acoustic cache reload
    axisym_mode = axisymmetric_output_mode(parameters);
    if ~isempty(axisym_mode)
        [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            acoustic_convert_axisymmetry(parameters, sensor_data, segmentation, medium_masks, ...
                                         kwave_medium, kgrid, source, source_labels, axisym_mode);
        if strcmp(axisym_mode, '2d')
            warn('Axisymmetry setup was requested for acoustic simulation without follow-up heating simulation. Follow-up thermal simulations loading these 2D results may not work as expected (e.g., 3D thermal diffusion).')
        end
    end

    % keep 'parameters' in info so not to confuse future runs when saving
    acoustic_info.parameters = parameters;
    acoustic_info.kwave_input_args = kwave_input_args;

    % Save p_complex field(s) to separate cache file(s) when explicitly requested
    % (one per unique frequency; WRITE_P_COMPLEX_FILES).
    if exist('p_complex_meta', 'var')
        sensor_data = write_p_complex_files(parameters, sensor_data, filename_sensor_data, p_complex_meta);
    end

    if ~should_save_output(parameters.io, 'save_acoustic_matrices')
        disp('Not saving acoustic output matrices ...')
    else
        acoustic_info.parameters = parameters;
        acoustic_info.kwave_input_args = kwave_input_args;
        cache_acoustic('save', filename_sensor_data, parameters, sensor_data, ...
                       kgrid, acoustic_provenance, acoustic_info, kwave_medium);
    end
end

function save_elastic_debug_maps(parameters, elastic_medium)
% [DEBUG] Save the simulation-grid shear speed / shear absorption maps.
    if ~debug_on(parameters, 'medium', 'volumes') || ~isfield(parameters.io, 'dir_debug_medium')
        return
    end
    try
        f_speed = fullfile(parameters.io.dir_debug_medium, ['matrix_sound_speed_shear' char(parameters.io.output_affix) '.nii']);
        niftiwrite(single(elastic_medium.sound_speed_shear), f_speed, 'Compressed', true);
        debug_register(parameters, 'medium', f_speed, 'Elastic skull: shear sound speed [m/s], simulation grid', 'volume');
        f_alpha = fullfile(parameters.io.dir_debug_medium, ['matrix_alpha_coeff_shear' char(parameters.io.output_affix) '.nii']);
        niftiwrite(single(elastic_medium.alpha_coeff_shear), f_alpha, 'Compressed', true);
        debug_register(parameters, 'medium', f_alpha, 'Elastic skull: shear attenuation, simulation grid', 'volume');
    catch ME
        warn('acoustic_wrapper:elasticDebugNifti', 'Could not save elastic debug maps: %s', ME.message);
    end
end
