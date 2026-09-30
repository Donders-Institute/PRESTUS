function cs = manual_coordinate_space(parameters)
% MANUAL_COORDINATE_SPACE  Coordinate space of transducer.trans_pos / focus_pos
%
% Returns the space in which transducer.trans_pos and transducer.focus_pos
% are currently expressed. Only placement.mode = 'manual' (the default)
% honours placement.manual.coordinate_space; every other placement mode
% writes voxel indices itself, so 'voxel' is returned for them.
%
% Use as:
%   cs = manual_coordinate_space(parameters)
%
% Input:
%   parameters - (1,1) simulation configuration struct
%
% Output:
%   cs - 'voxel' (1-based T1 voxel indices) | 'ras' (subject RAS+ mm) |
%        'mni' (MNI152 mm)
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, FOCAL_DISTANCE_CALCULATION

arguments
    parameters (1,1) struct
end

    cs = 'voxel';
    if ~isfield(parameters, 'placement')
        return;
    end
    pl = parameters.placement;

    mode = 'manual';
    if isfield(pl, 'mode') && ~isempty(pl.mode)
        mode = lower(char(pl.mode));
    end
    if ~strcmp(mode, 'manual') || ~isfield(pl, 'manual') || ...
            ~isfield(pl.manual, 'coordinate_space') || isempty(pl.manual.coordinate_space)
        return;
    end

    cs = lower(char(pl.manual.coordinate_space));
    if ~any(strcmp(cs, {'voxel', 'ras', 'mni'}))
        error('PRESTUS:placement:badCoordinateSpace', ...
            ['Unknown placement.manual.coordinate_space ''%s''. ' ...
             'Use ''voxel'', ''ras'', or ''mni''.'], cs);
    end
end
