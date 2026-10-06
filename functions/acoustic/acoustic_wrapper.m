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
    fprintf('[acoustic_wrapper] ---- k-Wave simulation summary ----\n');
    fprintf('[acoustic_wrapper]   backend   : %s\n', code_type);
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
    sensor_data = acoustic_simulation(kgrid, kwave_medium, source, sensor, kwave_input_args, parameters);

    % Extract complex pressure amplitude from the recorded steady-state time series.
    % p_complex = amp * exp(1i * phase) preserves interference information for
    % coherent multi-transducer superposition. For multi-frequency setups, the DFT
    % isolates each tone's component independently, so extractAmpPhase is run once
    % per unique frequency. A single-frequency setup produces the original
    % sensor_data.p_complex field; a multi-frequency setup produces one field per
    % unique frequency, keyed as p_complex_<freq_Hz>Hz.
    if isfield(sensor_data, 'p')
        freqs_all    = [parameters.transducer.freq_hz];
        unique_freqs = unique(freqs_all);
        if numel(unique_freqs) == 1
            [amp, phase] = extractAmpPhase(sensor_data.p, 1/kgrid.dt, unique_freqs(1), 'Dim', ndims(sensor_data.p));
            sensor_data.p_complex = amp .* exp(1i .* phase);
        else
            fprintf('Multi-frequency setup: extracting p_complex at %d frequencies (%s Hz).\n', ...
                numel(unique_freqs), num2str(unique_freqs));
            for fi = 1:numel(unique_freqs)
                [amp, phase] = extractAmpPhase(sensor_data.p, 1/kgrid.dt, unique_freqs(fi), 'Dim', ndims(sensor_data.p));
                fkey = sprintf('p_complex_%gHz', unique_freqs(fi));
                sensor_data.(fkey) = amp .* exp(1i .* phase);
            end
        end
        sensor_data = rmfield(sensor_data, 'p');
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
                p_complex = sensor_data.(fld); %#ok<NASGU>
                save(filename_p_complex, 'p_complex', '-v7.3');
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