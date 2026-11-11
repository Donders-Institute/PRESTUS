function v = get_nonempty_field(s, path, default)
% GET_NONEMPTY_FIELD  Value of a (nested) struct field, or a default when it is missing or empty
%
% As GET_FIELD, but an empty value ('' or []) also returns the default, as
% for config keys left empty to mean "use the default".
%
% Use as:
%   v = get_nonempty_field(parameters, {'transducer', 'condition'}, 'default')
%
% Input:
%   s       - struct
%   path    - field name (char) or field path (cell of char)
%   default - value returned when a level is missing or the value is empty
%             (default [])
%
% Output:
%   v - the field value, or default
%
% See also: GET_FIELD

    if nargin < 3
        default = [];
    end
    v = get_field(s, path, default);
    if isempty(v)
        v = default;
    end
end
