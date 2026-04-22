function idx = getidx(parameters_fields, string)
% GETIDX  Return tissue label indices matching a name string from a layers struct
%
% Searches the fieldnames of parameters_fields for entries that contain
% string and concatenates their values. Returns empty when nothing matches.
%
% Use as:
%   idx = getidx(parameters_fields, string)
%
% Input:
%   parameters_fields - struct whose fields map tissue names to label arrays
%                       (e.g., parameters.layers)
%   string            - tissue name to search for (e.g., 'skin')
%
% Output:
%   idx - concatenated label indices for all matching fields
%
% See also: CHARM_SEG_LABELS, CHECK_LAYERS

arguments
    parameters_fields (1,1) struct
    string            % char or cell array of char
end

    if ischar(string)
        string = {string};
    end

    % Initialize output array
    idx = [];

    % Get all field names from the parameter structure
    labels = fieldnames(parameters_fields);

    for si = 1:numel(string)
        idx_label = find(contains(labels, string{si}));
        for numidx = 1:numel(idx_label)
            idx = [idx, parameters_fields.(labels{idx_label(numidx)})];
        end
    end
end