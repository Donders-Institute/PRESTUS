function [sensor_data] = acoustic_simulation_elastic(kgrid, medium, source, sensor, input_args, parameters)
% ACOUSTIC_SIMULATION_ELASTIC  Run a k-Wave elastic (shear-wave) simulation via pstdElastic2D/3D
%
% Backend for acoustic_backend: 'kwave_elastic'. Models the skull as a
% solid supporting compressional and shear waves (Gao et al., 2026, IEEE
% TBME 73(3):1298), instead of the fluid skull of kspaceFirstOrder*.
%
% The PRESTUS fluid inputs are translated as follows:
%   medium  - if it does not already carry elastic fields, it is converted
%             by MEDIUM_ELASTIC_PROPERTIES without a skull mask (cS = 0;
%             e.g. the free-water baseline). Only the fields accepted by
%             pstdElastic are passed on.
%   source  - the pressure source (p_mask, p) becomes an equivalent normal
%             stress source: s_mask = p_mask, sxx = syy (= szz) = -p.
%             k-Wave scales both by 2*dt*c0/dx and pstdElastic records
%             p = -(sxx + syy + szz)/3, so the injected pressure matches the
%             fluid solver in a fluid region.
%   sensor  - unchanged (p_max_all, p_min_all, p_final, p are all supported).
%
% Only the MATLAB solvers exist for elastic waves: cpp_* code types fall
% back to their matlab_* equivalent. Axisymmetric grids are not supported.
% pstdElastic is a PSTD scheme without k-space correction; keep
% grid.source_cfl <= 0.3 (Gao et al. used 0.3; the PRESTUS default 0.15 is fine).
%
% Use as:
%   sensor_data = acoustic_simulation_elastic(kgrid, medium, source, sensor, input_args, parameters)
%
% Input:
%   kgrid      - kWaveGrid
%   medium     - elastic medium (from MEDIUM_ELASTIC_PROPERTIES) or fluid medium
%   source     - k-Wave pressure source (p_mask, p[, p_mode])
%   sensor     - k-Wave sensor (mask, record, record_start_index)
%   input_args - struct of k-Wave options (PMLInside, PMLSize, PlotSim, ...)
%   parameters - PRESTUS config (simulation.code_type, simulation.precision, grid)
%
% Output:
%   sensor_data - k-Wave sensor output (p_max_all, p_final, ...) [Pa]
%
% See also: ACOUSTIC_SIMULATION, MEDIUM_ELASTIC_PROPERTIES, PSTDELASTIC3D, PSTDELASTIC2D

arguments
    kgrid      (1,1)
    medium     (1,1) struct
    source     (1,1) struct
    sensor     (1,1) struct
    input_args (1,1) struct
    parameters (1,1) struct
end

    if isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        error('acoustic_simulation_elastic:axisymmetric', ...
            'The kwave_elastic backend does not support axisymmetric grids (set grid.axisymmetric: 0).');
    end

    %% Medium
    if ~isfield(medium, 'sound_speed_shear')
        medium = medium_elastic_properties(parameters, medium, []);
    end
    keep = {'sound_speed_compression', 'sound_speed_shear', 'density', ...
            'alpha_coeff_compression', 'alpha_coeff_shear'};
    medium = rmfield(medium, setdiff(fieldnames(medium), keep));

    %% Source: pressure -> normal stress
    if ~isfield(source, 'p_mask') || ~isfield(source, 'p')
        error('acoustic_simulation_elastic:source', ...
            'Expected a pressure source (source.p_mask, source.p).');
    end
    neg_p = source.p .* -1;   % CWSourceMatrix supports scalar times, not uminus
    elastic_source = struct('s_mask', source.p_mask, 'sxx', neg_p, 'syy', neg_p);
    if kgrid.dim == 3
        elastic_source.szz = neg_p;
    end
    if isfield(source, 'p_mode')
        mode = char(source.p_mode);
        % 'additive-no-correction' has no elastic equivalent
        if startsWith(mode, 'additive'), mode = 'additive'; end
        elastic_source.s_mode = mode;
    end
    clear neg_p source

    %% Solver options
    % C++/CUDA binaries do not implement the elastic equations.
    code_type = char(parameters.simulation.code_type);
    switch code_type
        case 'cpp_gpu'
            warning('acoustic_simulation_elastic:noCpp', ...
                'No C++/CUDA elastic solver in k-Wave; running pstdElastic with matlab_gpu instead.');
            code_type = 'matlab_gpu';
        case 'cpp_cpu'
            warning('acoustic_simulation_elastic:noCpp', ...
                'No C++ elastic solver in k-Wave; running pstdElastic with matlab_cpu instead.');
            code_type = 'matlab_cpu';
    end
    precision = char(parameters.simulation.precision);
    switch code_type
        case 'matlab_gpu', input_args.DataCast = ['gpuArray-', precision];
        case 'matlab_cpu', input_args.DataCast = precision;
        otherwise
            error('Unsupported code_type: %s. Supported options: matlab_cpu, matlab_gpu, cpp_cpu, cpp_gpu', code_type);
    end
    % Options of the fluid solvers that pstdElastic does not accept
    input_args = rmfield(input_args, intersect(fieldnames(input_args), ...
        {'SaveToDisk', 'DataName', 'DataPath', 'DeleteData', 'FunctionName', ...
         'DeviceNum', 'NumThreads', 'Axisymmetric'}));
    input_args_cell = zip_fields(input_args);

    fprintf('[acoustic_simulation_elastic] pstdElastic%dD  (%s, %s)\n', kgrid.dim, code_type, precision);
    switch kgrid.dim
        case 3
            sensor_data = pstdElastic3D(kgrid, medium, elastic_source, sensor, input_args_cell{:});
        case 2
            sensor_data = pstdElastic2D(kgrid, medium, elastic_source, sensor, input_args_cell{:});
        otherwise
            error('acoustic_simulation_elastic:dim', 'Elastic simulations require a 2D or 3D grid.');
    end
end
