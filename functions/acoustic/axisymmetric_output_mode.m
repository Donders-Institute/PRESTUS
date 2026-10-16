function mode = axisymmetric_output_mode(parameters)
% AXISYMMETRIC_OUTPUT_MODE  Form into which axisymmetric results are expanded
%
% Decides once, for both a fresh simulation (ACOUSTIC_WRAPPER) and a reloaded
% acoustic cache (ACOUSTIC_STAGE), whether the [Naxial x Nr] half-plane is
% expanded to 3D or mirrored into the bilateral 2D plane.
% grid.axisymmetric_output 'auto' (default): phantom simulations always expand
% to 3D so that output NIfTIs have consistent dimensions regardless of whether
% heating is requested; all other mediums expand to 3D only when heating is
% requested (3D thermal diffusion requires it). '2d' keeps the mirrored plane
% (no heating), '3d' always expands.
%
% Use as:
%   mode = axisymmetric_output_mode(parameters)
%
% Input:
%   parameters - PRESTUS config (grid.dims, grid.axisymmetric,
%                grid.axisymmetric_output, simulation.medium,
%                modules.run_heating_sims)
%
% Output:
%   mode - '3d' or '2d'; '' when the grid is not an axisymmetric half-plane
%          (grid.axisymmetric off or grid.dims not 2-D)
%
% See also: ACOUSTIC_CONVERT_AXISYMMETRY, CONVERT_AXISYMMETRIC_TO_3D,
%           CONVERT_AXISYMMETRIC_TO_2D

arguments
    parameters (1,1) struct
end

    mode = '';
    if ~(numel(parameters.grid.dims) == 2 && isfield(parameters.grid, 'axisymmetric') && ...
            parameters.grid.axisymmetric == 1)
        return
    end

    heating = isfield(parameters, 'modules') && isfield(parameters.modules, 'run_heating_sims') && ...
        parameters.modules.run_heating_sims == 1;
    out_mode = 'auto';
    if isfield(parameters.grid, 'axisymmetric_output') && ~isempty(parameters.grid.axisymmetric_output)
        out_mode = lower(char(parameters.grid.axisymmetric_output));
    end
    switch out_mode
        case 'auto'
            if strcmp(parameters.simulation.medium, 'phantom') || heating
                mode = '3d';
            else
                mode = '2d';
            end
        case '3d'
            mode = '3d';
        case '2d'
            if heating
                error('acoustic_wrapper:axisymOutput2dHeating', ...
                    'grid.axisymmetric_output = ''2d'' cannot be combined with heating simulations (3D required).');
            end
            mode = '2d';
        otherwise
            error('acoustic_wrapper:axisymOutput', ...
                'grid.axisymmetric_output must be ''auto'', ''2d'' or ''3d'' (got ''%s'').', out_mode);
    end
end
