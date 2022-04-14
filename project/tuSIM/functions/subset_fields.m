function res_struct = subset_fields(orig_struct, fieldnames)
% fieldnames is a cell array of fields to copy
    res_struct = struct();
    for x = fieldnames
        res_struct.(char(x)) = orig_struct.(char(x));
    end
end