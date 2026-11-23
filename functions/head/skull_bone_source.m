function source = skull_bone_source(parameters)
% SKULL_BONE_SOURCE  Where this run's skull layer comes from
%
% pct.bone_source: 'segmentation' (the segmentation's skull labels, holes
% filled by SKULL_FILL_HOLES), 'hu' (an HU threshold on the simulation grid,
% SKULL_FROM_HU; needs pct.enabled, else the segmentation is used), or 'mask'
% (a bone mask file, PCT_LOAD_BONE_MASK). 'hu' and 'mask' are reconciled with
% the segmentation by SKULL_RECONCILE.
%
% Use as:
%   source = skull_bone_source(parameters)
%
% Input:
%   parameters - PRESTUS parameters (pct.bone_source, pct.enabled)
%
% Output:
%   source - 'segmentation', 'hu' or 'mask'
%
% See also: SKULL_RECONCILE, SKULL_FROM_HU, PCT_LOAD_BONE_MASK

    source = char(get_field(parameters, {'pct', 'bone_source'}, 'segmentation'));
    switch source
        case 'segmentation'
        case 'hu'
            if ~isequal(get_field(parameters, {'pct', 'enabled'}, 0), 1)
                source = 'segmentation';
            end
        case 'mask'
        otherwise
            error('skull_bone_source:badSource', ...
                'pct.bone_source must be ''segmentation'', ''hu'' or ''mask'' (got ''%s'').', source);
    end
end
