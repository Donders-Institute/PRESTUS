function validate_any_required_field(s, candidates, id, label)
% VALIDATE_ANY_REQUIRED_FIELD  Error clearly unless at least one candidate field is set.
%
% Use for config values with more than one valid source (e.g. an explicit
% override vs. a value derivable from a sibling config block), where at
% least one, but not necessarily all, must be provided.
%
% Use as:
%   validate_any_required_field(parameters.transducer(1), ...
%       {'focal_distance_ep', 'focal_distance_bowl'}, ...
%       'PRESTUS:localite:missingConfig', 'transducer(1)')
%
% Input:
%   s          - struct to check
%   candidates - cell array of field name strings; at least one must be
%                present and non-empty
%   id         - MException identifier
%   label      - dotted config path prefix used in the error message
%
% See also: VALIDATE_REQUIRED_FIELDS

    arguments
        s (1,1) struct
        candidates (1,:) cell
        id (1,:) char
        label (1,:) char
    end

    has_any = any(cellfun(@(f) isfield(s, f) && ~isempty(s.(f)), candidates));

    if ~has_any
        error(id, 'At least one of the following field(s) must be set: %s', strjoin( ...
            cellfun(@(f) sprintf('%s.%s', label, f), candidates, 'UniformOutput', false), ', '));
    end
end
