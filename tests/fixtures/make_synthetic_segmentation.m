function seg = make_synthetic_segmentation(dims)
% MAKE_SYNTHETIC_SEGMENTATION  Create a synthetic SimNIBS-style segmentation volume.
%
%   seg = make_synthetic_segmentation(dims)
%
%   Produces a [dims(1) x dims(2) x dims(3)] uint8 label volume with the
%   following concentric shell structure (charm label values):
%
%     0  — background (water)
%     5  — skin        (outer shell, 3 voxels thick)
%     4  — skull       (2 voxels thick)
%     3  — CSF         (2 voxels thick)
%     2  — grey matter (2 voxels thick)
%     1  — white matter (core)
%
%   Useful for testing head preprocessing functions without real MRI data.

    if nargin < 1
        dims = [64, 64, 80];
    end

    seg = zeros(dims, 'uint8');
    c   = round(dims / 2);   % centre voxel

    % Radii (voxels) for each shell boundary
    r_skin_outer = min(dims)/2 - 2;
    r_skull_outer = r_skin_outer  - 3;
    r_csf_outer   = r_skull_outer - 2;
    r_gm_outer    = r_csf_outer   - 2;
    r_wm_outer    = r_gm_outer    - 2;

    [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    r = sqrt((X - c(1)).^2 + (Y - c(2)).^2 + (Z - c(3)).^2);

    seg(r <= r_skin_outer)  = 5;   % skin
    seg(r <= r_skull_outer) = 4;   % skull
    seg(r <= r_csf_outer)   = 3;   % CSF
    seg(r <= r_gm_outer)    = 2;   % grey matter
    seg(r <= r_wm_outer)    = 1;   % white matter

end
