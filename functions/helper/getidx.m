function idx = getidx(parameters_fields, string)

% GETIDX Retrieves indices for requested tissues from a parameter structure.
%
% This function searches for tissue labels in a parameter structure (`parameters_fields`) 
% that match a specified string (`string`). It returns the indices associated with 
% the matching labels.
%
% Input:
%   parameters_fields - Struct containing tissue labels and their associated indices.
%                       Example: `parameters.layer_labels`.
%   string            - String specifying the tissue name to search for.
%                       Example: `'skin'`.
%
% Output:
%   idx               - Array of indices corresponding to the matching tissue labels.
%
% Notes:
%   - If no matching labels are found, the function returns an empty array.

    % Initialize output array
    idx = [];

    % Get all field names from the parameter structure
    labels = fieldnames(parameters_fields);

    % Find indices of field names that contain the specified string
    idx_label = find(contains(labels, string));

    % If matching labels are found, retrieve their indices
    if ~isempty(idx_label)
        for numidx = 1:numel(idx_label)
            idx = [idx, parameters_fields.(labels{idx_label(numidx)})];
        end
    end
end