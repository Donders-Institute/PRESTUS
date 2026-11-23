function hu_max = medium_pct_hu_max(parameters, hu_skull)
% MEDIUM_PCT_HU_MAX  Densest-bone HU for the max-normalised pCT mappings
%
% Returns the whole-head reference that PCT_MAP_PROPERTIES resolved with
% PCT_HU_REFERENCE (held in parameters.pct.hu_reference_value for the
% duration of the mapping). Without it, e.g. for mapping curves over a
% synthetic HU range, the largest of the given skull HU.
%
% Use as:
%   hu_max = medium_pct_hu_max(parameters, hu_skull)
%
% Input:
%   parameters - PRESTUS config
%   hu_skull   - HU of the skull voxels being mapped
%
% Output:
%   hu_max - reference HU
%
% See also: PCT_HU_REFERENCE, PCT_MAP_PROPERTIES, MEDIUM_PCT_DENSITY

    arguments
        parameters (1,1) struct
        hu_skull   {mustBeNumeric}
    end
    hu_max = get_field(parameters, {'pct', 'hu_reference_value'}, []);
    if isempty(hu_max)
        hu_max = double(max(hu_skull(:)));
        return
    end
    hu_max = double(hu_max);
    data_max = double(max(hu_skull(:)));
    if data_max > hu_max
        % Expected with pct.hu_reference_percentile below 100; the mappings
        % limit these voxels to the densest-bone value (MEDIUM_PCT_POROSITY)
        fprintf('[medium_pct_hu_max] %d skull voxels above the HU reference (%.0f HU, up to %.0f HU) map as the densest bone.\n', ...
            nnz(hu_skull > hu_max), hu_max, data_max);
    end
end
