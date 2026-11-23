function phi = medium_pct_porosity(parameters, hu)
% MEDIUM_PCT_POROSITY  Skull porosity from HU for the porosity mixture mappings
%
% phi = 1 - HU / HU_ref (Aubry et al., 2003), with HU_ref the densest-bone
% reference (MEDIUM_PCT_HU_MAX, PCT_HU_REFERENCE), limited to [0, 1]: HU
% above the reference map to the densest bone (phi = 0) and HU below 0 to
% water (phi = 1). With pct.hu_reference_percentile below 100, some skull
% voxels lie above the reference; the limit keeps their density, sound speed,
% and attenuation within the water-to-bone range instead of extrapolating.
%
% Use as:
%   phi = medium_pct_porosity(parameters, hu)
%
% Input:
%   parameters - PRESTUS config
%   hu         - HU of the skull voxels being mapped
%
% Output:
%   phi - porosity in [0, 1], same size as hu
%
% See also: MEDIUM_PCT_HU_MAX, MEDIUM_PCT_DENSITY, MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION

    arguments
        parameters (1,1) struct
        hu         {mustBeNumeric}
    end
    phi = min(max(1 - double(hu) / medium_pct_hu_max(parameters, hu), 0), 1);
end
