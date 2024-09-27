function idx = getidx(parameters_fields, string)
% This function retrieves the indices for requested tissues
% e.g., parameters_fields = parameters.layer_labels;
% e.g., string = 'skin';
    idx = [];
    labels = fieldnames(parameters_fields);
    idx_label = find(contains(labels, string));
    if ~isempty(idx_label)
        for numidx = 1:numel(idx_label)
            idx = [idx, parameters_fields.(labels{idx_label(numidx)})];
        end
    end
end