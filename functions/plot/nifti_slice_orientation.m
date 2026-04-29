function [do_transpose, flip_rows, flip_cols] = nifti_slice_orientation(hdr, slice_axis)
% NIFTI_SLICE_ORIENTATION  Compute display transforms for a NIfTI orthogonal slice.
%
% Given a NIfTI header (from niftiinfo) and the axis that is held fixed,
% returns the three boolean flags needed to orient the resulting 2D slice
% for standard neurological display: superior at top, patient left on
% viewer's left.
%
% Use as:
%   [do_transpose, flip_rows, flip_cols] = nifti_slice_orientation(hdr, slice_axis)
%
% Input:
%   hdr        - NIfTI info struct from niftiinfo (must contain Transform.T)
%   slice_axis - char, 'x', 'y', or 'z' — axis that is fixed (sliced through)
%
% Output:
%   do_transpose - logical; transpose the 2D slice (swap row/col dims)
%   flip_rows    - logical; flipud (puts superior at top, row 1)
%   flip_cols    - logical; fliplr (puts patient left on viewer left)
%
% Background:
%   MATLAB's squeeze on an orthogonal slice places the first free voxel
%   dimension in rows and the second in columns.  For a y-slice of a
%   standard RAS volume the rows encode R-L and the columns encode I-S,
%   which is 90° off from a standard coronal view.  This function reads
%   the affine to decide whether a transpose and/or flips are needed so
%   the result is always orientation-independent.
%
%   The Transform.T convention follows MATLAB's niftiinfo: T' (i.e.
%   hdr.Transform.T transposed) is the 4x4 voxel-to-world matrix in the
%   standard column-vector sense, so row d of T gives the world-space
%   direction for increasing voxel index along dim d.
%
% See also: PLOT_OVERLAY, NIFTIINFO

    % Row d of T → world-space direction for voxel dim d
    T = hdr.Transform.T;

    % The two voxel dims that remain after fixing slice_axis
    switch lower(slice_axis)
        case 'x',  free = [2, 3];
        case 'y',  free = [1, 3];
        case 'z',  free = [1, 2];
        otherwise, error('slice_axis must be ''x'', ''y'', or ''z''.');
    end

    d1 = free(1);   % → rows after squeeze
    d2 = free(2);   % → cols after squeeze

    % World-space direction vectors (x=R, y=A, z=S in RAS)
    v1 = T(d1, 1:3);   % direction for row axis
    v2 = T(d2, 1:3);   % direction for col axis

    % Put the dimension with the largest z-component (I-S) in the row axis
    % so that a simple flipud can put superior at the top.
    if abs(v2(3)) > abs(v1(3))
        do_transpose = true;
        row_vec = v2;
        col_vec = v1;
    else
        do_transpose = false;
        row_vec = v1;
        col_vec = v2;
    end

    % flip_rows: superior (positive z) should be at row 1 (top of imagesc).
    % If row_vec(3) > 0 then row 1 is the most inferior voxel → flipud needed.
    flip_rows = row_vec(3) > 0;

    % flip_cols: patient's left (negative x in RAS) should be on viewer's left.
    % If col_vec(1) < 0 then col 1 is the rightmost voxel → fliplr needed.
    flip_cols = col_vec(1) < 0;

end
