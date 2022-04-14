function cell_arr = zip_fields(input_struct)
    cell_arr = {};
    fields = fieldnames(input_struct);
    for i = 1:length(fields)
        fname = fields{i};
        cell_arr = [cell_arr, fname];
        cell_arr = [cell_arr, input_struct.(fname)];
    end
end
