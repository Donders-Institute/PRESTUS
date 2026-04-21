function cell_arr = zip_fields(input_struct)

% ZIP_FIELDS  Convert a struct into an alternating {field, value, field, value, ...} cell array
%
% Useful for passing struct contents as name-value pairs to functions that
% accept varargin-style arguments.
%
% Use as:
%   cell_arr = zip_fields(input_struct)
%
% Input:
%   input_struct - scalar struct with any field types
%
% Output:
%   cell_arr - [1×2N] cell array: {field1, value1, field2, value2, ...}
%
% See also: FIELDNAMES, STRUCT2CELL

    % Initialize the output cell array
    cell_arr = {};

    % Get all field names from the input structure
    fields = fieldnames(input_struct);

    % Iterate through each field and append its name and value to the cell array
    for i = 1:length(fields)
        fname = fields{i};
        cell_arr = [cell_arr, fname]; % Append field name
        cell_arr = [cell_arr, input_struct.(fname)]; % Append field value
    end
end