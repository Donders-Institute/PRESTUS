function pml = resolve_pml_size(pml_size, kgrid)
% RESOLVE_PML_SIZE  Resolve grid.pml_size to a numeric array for k-Wave
%
% yaml.loadFile returns MATLAB string objects (not char), so the 'auto'
% value arrives as string("auto"). This function handles char, string, and
% numeric inputs uniformly and resolves 'auto' via getOptimalPMLSize.
%
% Use as:
%   pml = resolve_pml_size(pml_size, kgrid)
%
% Input:
%   pml_size - parameters.grid.pml_size (numeric, char, or string)
%   kgrid    - kWaveGrid object (needed only when pml_size == 'auto')
%
% Output:
%   pml - double scalar or vector suitable for the k-Wave PMLSize argument

if (ischar(pml_size) || isstring(pml_size)) && strcmp(pml_size, 'auto')
    axisym = nargin >= 2 && ~isempty(kgrid) && ...
             isprop(kgrid, 'dim') && kgrid.dim == 2;
    if axisym
        pml = getOptimalPMLSize(kgrid, [], 'WSWA');
    else
        pml = getOptimalPMLSize(kgrid);
    end
else
    pml = double(pml_size);
end

end
