function [sig, skull_mapping] = pct_signature(parameters)
% PCT_SIGNATURE  Settings a pseudo-CT depends on, and the resolved skull mapping
%
% Recorded in the pseudo-CT's sidecar (pseudoCT.json) when it is generated and
% compared when an existing pseudo-CT is reused (PREPROC_SEGMENTATION). Only
% the UTE→HU skull mapping changes pseudoCT.nii.gz; the pct.mapping_* keys act
% later, in the medium stage, and pct.debug only keeps intermediate files.
%
% Use as:
%   [sig, skull_mapping] = pct_signature(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (pct.skull_mapping, pct.custom_*)
%
% Output:
%   sig           - struct: signature_version, skull_mapping
%   skull_mapping - mapping passed to pct_create_pseudoCT.sh ('kosciessa', ...,
%                   or 'custom:<slope>:<intercept>')
%
% See also: PCT_CREATE_PCT_RUN, CHECK_SIDECAR_SIGNATURE

arguments
    parameters (1,1) struct
end

    skull_mapping = 'kosciessa';
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'skull_mapping') ...
            && ~isempty(parameters.pct.skull_mapping)
        skull_mapping = char(parameters.pct.skull_mapping);
    end
    if strcmp(skull_mapping, 'custom')
        % user-defined linear map, passed to the shell script as custom:<slope>:<intercept>
        assert(isfield(parameters.pct, 'custom_slope') && isfield(parameters.pct, 'custom_intercept') && ...
            isscalar(parameters.pct.custom_slope) && isscalar(parameters.pct.custom_intercept), ...
            'pct.skull_mapping = ''custom'' needs scalar pct.custom_slope and pct.custom_intercept');
        skull_mapping = sprintf('custom:%.10g:%.10g', parameters.pct.custom_slope, parameters.pct.custom_intercept);
    end

    sig = struct('signature_version', 1, 'skull_mapping', skull_mapping);
end
