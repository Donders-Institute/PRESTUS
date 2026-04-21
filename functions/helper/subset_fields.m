function res_struct = subset_fields(orig_struct, fieldnames)
% SUBSET_FIELDS  Copy a named subset of fields from one struct to a new struct
%
% Use as:
%   res_struct = subset_fields(orig_struct, fieldnames)
%
% Input:
%   orig_struct - source struct
%   fieldnames  - cell array of field name strings to copy
%
% Output:
%   res_struct  - new struct containing only the requested fields
%
% See also: RMFIELD, FIELDNAMES

arguments
    orig_struct (1,1) struct
    fieldnames  (1,:) cell
end
    res_struct = struct();
    for x = fieldnames
        res_struct.(char(x)) = orig_struct.(char(x));
    end
end