function v = get_field(s, path, default)
% GET_FIELD  Value of a (nested) struct field, or a default when it is missing
%
% Reads s.a.b.c without erroring when a level is missing or not a struct.
% An empty value is returned as is; use GET_NONEMPTY_FIELD to treat empty
% as missing. A struct array is read from its first element.
%
% Use as:
%   v = get_field(parameters, 'platform', 'matlab')
%   v = get_field(parameters, {'grid', 'resolution_mm'}, [])
%
% Input:
%   s       - struct
%   path    - field name (char) or field path (cell of char)
%   default - value returned when a level is missing (default [])
%
% Output:
%   v - the field value, or default
%
% See also: GET_NONEMPTY_FIELD

    if nargin < 3
        default = [];
    end
    v = default;
    if ischar(path) || isstring(path)
        path = cellstr(path);
    end
    for k = 1:numel(path)
        if ~isstruct(s) || isempty(s) || ~isfield(s, path{k})
            return
        end
        s = s(1).(path{k});
    end
    v = s;
end
