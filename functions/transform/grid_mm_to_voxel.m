function vox = grid_mm_to_voxel(pos_mm, resolution_mm)
% GRID_MM_TO_VOXEL  Grid position in mm to a 1-based voxel index
%
% transducer.trans_pos_mm / focus_pos_mm are measured from the centre of grid
% voxel [1,1,1] (the voxel grid.origin_ras_mm refers to), so 0 mm is voxel 1
% and voxel i lies at (i - 1) * resolution_mm.
%
% Use as:
%   vox = grid_mm_to_voxel(pos_mm, resolution_mm)
%
% Input:
%   pos_mm        - [1xN] position [mm] from the centre of voxel 1
%   resolution_mm - grid resolution [mm]
%
% Output:
%   vox           - [1xN] voxel indices (row vector)
%
% See also: GRID_TRANSDUCER_LOCATION, GRID_TISSUE_SETUP

    vox = round(double(pos_mm(:)') / resolution_mm) + 1;
end
