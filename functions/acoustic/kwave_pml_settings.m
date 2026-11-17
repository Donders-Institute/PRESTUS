function [pml, pml_alpha] = kwave_pml_settings(parameters, kgrid)
% KWAVE_PML_SETTINGS  PML size and absorption for a k-Wave run (PML outside the grid)
%
% Resolves grid.pml_size: a number is used as given; 'auto' is resolved with
% getOptimalPMLSize over grid.pml_size_range (k-Wave's [10, 40] when unset),
% with the 'WSWA' option on axisymmetric grids (grid.axisymmetric = 1).
% pml_alpha is grid.pml_alpha, or [] for k-Wave's default (2).
%
% Use as:
%   [pml, pml_alpha] = kwave_pml_settings(parameters, kgrid)
%
% Input:
%   parameters - PRESTUS config (grid.pml_size, grid.pml_size_range,
%                grid.pml_alpha, grid.axisymmetric)
%   kgrid      - kWaveGrid of the run
%
% Output:
%   pml       - PML size [grid points], scalar or one value per dimension
%   pml_alpha - PML absorption [Nepers per grid point], or []
%
% See also: ACOUSTIC_WRAPPER, REVERSE_POINT_SOURCE_FIELD, RESOLVE_JWAVE_PML_SIZE

pml = parameters.grid.pml_size;
if (ischar(pml) || isstring(pml)) && strcmp(pml, 'auto')
    axisym = isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1;
    pml_range = [];   % k-Wave default [10, 40]
    if isfield(parameters.grid, 'pml_size_range'); pml_range = parameters.grid.pml_size_range; end
    if axisym
        pml = getOptimalPMLSize(kgrid, pml_range, 'WSWA');
    else
        pml = getOptimalPMLSize(kgrid, pml_range);
    end
else
    pml = double(pml);
end

pml_alpha = [];
if isfield(parameters.grid, 'pml_alpha') && ~isempty(parameters.grid.pml_alpha)
    pml_alpha = parameters.grid.pml_alpha;
end
end
