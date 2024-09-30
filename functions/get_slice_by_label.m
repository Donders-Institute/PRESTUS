function out_img = get_slice_by_label(img, slice_label, slice_n)
    slice_x = 1:size(img,1);
    slice_y = 1:size(img,2);
    slice_z = 1:size(img,3);

   
    if strcmp(slice_label, 'x')
        slice_x = slice_n;
    elseif strcmp(slice_label, 'y')
        slice_y = slice_n;
    elseif strcmp(slice_label, 'z')
        slice_z = slice_n;
    else
        error("slice must be a cell array with the first element of 'x','y', or 'z' and a second element a slice number")
    end
    out_img = squeeze(img(slice_x, slice_y, slice_z));

end