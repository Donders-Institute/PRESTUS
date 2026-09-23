function validate_required_fields(s, required, id, label)
% VALIDATE_REQUIRED_FIELDS  Error clearly unless every required field is set.
%
% Checks each field name in `required` against struct `s` and raises a
% single error listing everything that is missing or empty, instead of
% failing on the first miss.
%
% Use as:
%   validate_required_fields(parameters.placement.plantus, ...
%       {'script_path', 'env_path', 'target_name'}, ...
%       'PRESTUS:plantus:missingConfig', 'placement.plantus')
%
% Input:
%   s        - struct to check
%   required - cell array of field name strings
%   id       - MException identifier, e.g. 'PRESTUS:plantus:missingConfig'
%   label    - dotted config path prefix used in the error message,
%              e.g. 'placement.plantus'
%
% See also: VALIDATE_ANY_REQUIRED_FIELD

    arguments
        s (1,1) struct
        required (1,:) cell
        id (1,:) char
        label (1,:) char
    end

    missing = required(~cellfun(@(f) isfield(s, f) && ~isempty(s.(f)), required));

    if ~isempty(missing)
        error(id, 'Missing required field(s): %s', strjoin( ...
            cellfun(@(f) sprintf('%s.%s', label, f), missing, 'UniformOutput', false), ', '));
    end
end
