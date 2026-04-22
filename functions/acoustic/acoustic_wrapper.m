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
    source_labels)
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
    parameters          (1,1) struct
    kgrid               (1,1)
    kwave_medium        (1,1) struct
    source              (1,1) struct
    sensor              (1,1) struct
    medium_masks        (:,:) {mustBeNumeric}
    filename_sensor_data(1,:) char
    segmentation        (:,:) {mustBeNumeric}
    source_labels       (:,:) {mustBeNumeric}
end

    disp('Specifying and Starting acoustic simulations...')

    % PML is always applied outside the user-defined grid (PMLInside=false).
    % The returned sensor data (p_max_all, p_final) is automatically cropped
    % by k-Wave to the original grid dimensions; callers never see the PML voxels.
    % grid.pml_size may be a positive integer or the string 'auto' (k-Wave selects
    % the optimal size for FFT efficiency; only valid with PMLInside=false).
    %
    % Resolve 'auto' to the actual size via getOptimalPMLSize so the effective
    % value is available for reporting and reproducibility.
    if ischar(parameters.grid.pml_size) && strcmp(parameters.grid.pml_size, 'auto')
        axisym = isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1;
        if axisym
            effective_pml = getOptimalPMLSize(kgrid, [], 'WSWA');
        else
            effective_pml = getOptimalPMLSize(kgrid);
        end
        parameters.grid.pml_size_effective = effective_pml;
        fprintf('PML size (auto-selected): [%s]\n', num2str(effective_pml));
    else
        effective_pml = parameters.grid.pml_size;
        parameters.grid.pml_size_effective = effective_pml;
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
        skull_edge = edge3(mask.skull, 'approxcanny', 0.1);
        % ... to set as display mask
        kwave_input_args.DisplayMask = skull_edge;
    end

    if parameters.modules.run_source_setup==0
        error('Source setup not requested. Not able to proceed with acoustic simulation.')
    end

    % perform the acoustic simulation
    sensor_data = acoustic_simulation(kgrid, kwave_medium, source, sensor, kwave_input_args, parameters);

    % convert media and results to 2D/3D (if axisymmetry was used)
    if numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        % Phantom simulations always expand to 3D so that output NIfTIs have
        % consistent dimensions regardless of whether heating is requested.
        % All other mediums expand to 3D only when heating is requested
        % (3D thermal diffusion requires it); otherwise retain 2D output.
        expand_to_3d = strcmp(parameters.simulation.medium, 'phantom') || ...
                       (isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims == 1);
        if expand_to_3d
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels);
        else
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels);
            warning('Axisymmetry setup was requested for acoustic simulation without follow-up heating simulation. Follow-up thermal simulations loading these 2D results may not work as expected (e.g., 3D thermal diffusion).')
        end
    end

    % keep 'parameters' in info so not to confuse future runs when saving
    acoustic_info.parameters = parameters;
    acoustic_info.kwave_input_args = kwave_input_args;

    if ~should_save_output(parameters.io, 'save_acoustic_matrices')
        disp('Not saving acoustic output matrices ...')
    else
        % keep 'parameters' as a copy so not to confuse future runs
        acoustic_info.parameters = parameters;
        acoustic_info.kwave_input_args = kwave_input_args;
        save(filename_sensor_data, ...
            'sensor_data', ...
            'kgrid', ...
            'kwave_medium', ...
            'source', ...
            'sensor', ...
            'segmentation', ...
            'source_labels', ...
            'medium_masks', ...
            'acoustic_info' ,'-v7.3')
    end