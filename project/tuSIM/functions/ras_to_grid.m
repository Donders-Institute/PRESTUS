function [grid_pos] = ras_to_grid(ras_pos, nii_header)
    grid_pos = round(nii_header.Transform.T' \ [ras_pos; 1]);
    grid_pos = grid_pos(1:3);
end



