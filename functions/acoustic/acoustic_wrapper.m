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
% expand to a bilateral 2D cross-section instead. The full set of simulation
% matrices is optionally saved to disk as a -v7.3 .mat file.
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, kwave_medium, ...
%    kgrid, source, source_labels] = ...
%       acoustic_wrapper(parameters, kgrid, kwave_medium, source, sensor, ...
%                        medium_masks, filename_sensor_data, segmentation, source_labels)
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
    pml_size_val = parameters.grid.pml_size;
    if (ischar(pml_size_val) || isstring(pml_size_val)) && strcmp(pml_size_val, 'auto')
        axisym = isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1;
        if axisym
            effective_pml = getOptimalPMLSize(kgrid, [], 'WSWA');
        else
            effective_pml = getOptimalPMLSize(kgrid);
        end
        parameters.grid.pml_size_effective = effective_pml;
    else
        effective_pml = double(pml_size_val);
        parameters.grid.pml_size_effective = effective_pml;
    end

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
            fprintf('[acoustic_wrapper]      Recommended PBS mem ≥ current_workspace + %.0f GB.\n', ceil(gpu_mem_gb));
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

    % Record the steady-state time series only when the complex pressure field
    % is explicitly requested — storing p for an all-grid sensor is large.
    save_p_complex = isfield(parameters.io, 'save_p_complex') && parameters.io.save_p_complex;
    if save_p_complex && ~ismember('p', sensor.record)
        sensor.record{end+1} = 'p';
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
    sensor_data = acoustic_simulation(kgrid, sim_medium, source, sensor, kwave_input_args, parameters);
    clear sim_medium

    % Extract complex pressure amplitude from the recorded steady-state time series.
    % p_complex = amp * exp(1i * phase) preserves interference information for
    % coherent multi-transducer superposition. For multi-frequency setups, the DFT
    % isolates each tone's component independently, so extractAmpPhase is run once
    % per unique frequency. A single-frequency setup produces the original
    % sensor_data.p_complex field; a multi-frequency setup produces one field per
    % unique frequency, keyed as p_complex_<freq_Hz>Hz.
    % extractAmpPhase returns the phase relative to the first recorded sample;
    % subtracting 2*pi*f*t0 refers it to t = 0 of the simulation, where every
    % run's source signal starts. Fields from separate runs (e.g. one per
    % transducer, possibly on FOV sub-grids with a different dt) can then be
    % added coherently (SUPERPOSE_COMPLEX_FIELDS).
    if isfield(sensor_data, 'p')
        freqs_all    = [parameters.transducer.freq_hz];
        unique_freqs = unique(freqs_all);
        t0_s = 0;
        if isfield(sensor, 'record_start_index') && ~isempty(sensor.record_start_index)
            t0_s = (double(sensor.record_start_index) - 1) * kgrid.dt;
        end
        p_complex_keys = cell(1, numel(unique_freqs));
        if numel(unique_freqs) > 1
            fprintf('Multi-frequency setup: extracting p_complex at %d frequencies (%s Hz).\n', ...
                numel(unique_freqs), num2str(unique_freqs));
        end
        for fi = 1:numel(unique_freqs)
            [amp, phase] = extractAmpPhase(sensor_data.p, 1/kgrid.dt, unique_freqs(fi), 'Dim', ndims(sensor_data.p));
            % k-Wave returns p as [sensor points x time]; the sensor mask is the
            % whole grid, so the points are in linear-index order.
            if numel(amp) == prod(parameters.grid.dims) && ~isequal(size(amp), parameters.grid.dims)
                amp   = reshape(amp,   parameters.grid.dims);
                phase = reshape(phase, parameters.grid.dims);
            end
            if numel(unique_freqs) == 1
                p_complex_keys{fi} = 'p_complex';
            else
                p_complex_keys{fi} = sprintf('p_complex_%gHz', unique_freqs(fi));
            end
            sensor_data.(p_complex_keys{fi}) = amp .* exp(1i .* (phase - 2 * pi * unique_freqs(fi) * t0_s));
        end
        sensor_data = rmfield(sensor_data, 'p');
        % Saved with each p_complex file (one entry per frequency below).
        p_complex_meta = struct( ...
            'phase_reference', 'simulation_t0', ...
            'freq_hz',         unique_freqs, ...
            'field_names',     {p_complex_keys}, ...
            'record_t0_s',     t0_s, ...
            'dt_s',            kgrid.dt, ...
            'drive_level',     'simulated', ...   % before any ISPPA rescaling
            'grid_dims',       [], ...
            'fov_offset_ac',   [], ...
            'full_ac_dims',    []);
        if nargin >= 10 && isstruct(acoustic_provenance) && isfield(acoustic_provenance, 'fov_offset_ac')
            p_complex_meta.fov_offset_ac = acoustic_provenance.fov_offset_ac;
            p_complex_meta.full_ac_dims  = acoustic_provenance.full_ac_dims;
        end
    end

    % convert media and results to 2D/3D (if axisymmetry was used)
    if numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        % grid.axisymmetric_output 'auto' (default): phantom simulations always
        % expand to 3D so that output NIfTIs have consistent dimensions regardless
        % of whether heating is requested; all other mediums expand to 3D only when
        % heating is requested (3D thermal diffusion requires it). '2d' keeps the
        % mirrored plane (no heating), '3d' always expands.
        heating = isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims == 1;
        out_mode = 'auto';
        if isfield(parameters.grid, 'axisymmetric_output') && ~isempty(parameters.grid.axisymmetric_output)
            out_mode = lower(char(parameters.grid.axisymmetric_output));
        end
        switch out_mode
            case 'auto'; expand_to_3d = strcmp(parameters.simulation.medium, 'phantom') || heating;
            case '3d';   expand_to_3d = true;
            case '2d'
                if heating
                    error('acoustic_wrapper:axisymOutput2dHeating', ...
                        'grid.axisymmetric_output = ''2d'' cannot be combined with heating simulations (3D required).');
                end
                expand_to_3d = false;
            otherwise
                error('acoustic_wrapper:axisymOutput', ...
                    'grid.axisymmetric_output must be ''auto'', ''2d'' or ''3d'' (got ''%s'').', out_mode);
        end
        if expand_to_3d
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels);
        else
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels);
            warn('Axisymmetry setup was requested for acoustic simulation without follow-up heating simulation. Follow-up thermal simulations loading these 2D results may not work as expected (e.g., 3D thermal diffusion).')
        end
    end

    % keep 'parameters' in info so not to confuse future runs when saving
    acoustic_info.parameters = parameters;
    acoustic_info.kwave_input_args = kwave_input_args;

    % Save p_complex field(s) to separate cache file(s) when explicitly requested.
    % Single-frequency: saves p_complex as before.
    % Multi-frequency: saves one file per unique frequency with a _<freq>Hz suffix.
    p_complex_fields = fieldnames(sensor_data);
    p_complex_fields = p_complex_fields(startsWith(p_complex_fields, 'p_complex'));
    if ~isempty(p_complex_fields)
        [fdir, fname, fext] = fileparts(filename_sensor_data);
        for fi = 1:numel(p_complex_fields)
            fld = p_complex_fields{fi};
            if isfield(parameters.io, 'save_p_complex') && parameters.io.save_p_complex
                % For backward compatibility, the single-frequency file keeps the original name.
                if strcmp(fld, 'p_complex')
                    filename_p_complex = fullfile(fdir, [strrep(fname, '_results', '_p_complex'), fext]);
                else
                    freq_suffix = fld(length('p_complex')+1:end); % e.g. '_500Hz'
                    filename_p_complex = fullfile(fdir, [strrep(fname, '_results', ['_p_complex' freq_suffix]), fext]);
                end
                p_complex = sensor_data.(fld);
                p_complex_info = rmfield(p_complex_meta, 'field_names');
                p_complex_info.freq_hz   = p_complex_meta.freq_hz(strcmp(p_complex_meta.field_names, fld));
                p_complex_info.grid_dims = size(p_complex);
                save(filename_p_complex, 'p_complex', 'p_complex_info', '-v7.3');
                fprintf('Saved complex pressure field (%s) to %s\n', fld, filename_p_complex);
            end
            sensor_data = rmfield(sensor_data, fld);
        end
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
    if parameters.simulation.debug ~= 1 || ~isfield(parameters.io, 'dir_debug_medium')
        return
    end
    try
        niftiwrite(single(elastic_medium.sound_speed_shear), ...
            fullfile(parameters.io.dir_debug_medium, ['matrix_sound_speed_shear' char(parameters.io.output_affix)]), ...
            'Compressed', true);
        niftiwrite(single(elastic_medium.alpha_coeff_shear), ...
            fullfile(parameters.io.dir_debug_medium, ['matrix_alpha_coeff_shear' char(parameters.io.output_affix)]), ...
            'Compressed', true);
    catch ME
        warn('acoustic_wrapper:elasticDebugNifti', 'Could not save elastic debug maps: %s', ME.message);
    end
end
