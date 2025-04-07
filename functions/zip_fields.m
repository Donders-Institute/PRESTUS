function cell_arr = zip_fields(input_struct)

% ZIP_FIELDS Converts a structure's fields and values into a cell array.
%
% This function takes a structure (`input_struct`) and creates a cell array where 
% each pair of elements corresponds to a field name and its associated value. 
%
% Input:
%   input_struct - Struct whose fields and values will be converted into a cell array.
%
% Output:
%   cell_arr     - Cell array containing field names and their corresponding values in alternating order.

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