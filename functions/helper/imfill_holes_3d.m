function out = imfill_holes_3d(BW)
% IMFILL_HOLES_3D  Fill interior holes in a 3-D binary volume
%
% Drop-in replacement for imfill(BW, 'holes') from the Image Processing
% Toolbox. Background voxels (BW==0) that are not 6-connected to the image
% boundary are considered holes and set to 1.
%
% Use as:
%   out = imfill_holes_3d(BW)
%
% Input:
%   BW  - logical 3-D array
%
% Output:
%   out - logical array, same size as BW, with holes filled
%
% See also: SKULL_RUBBER_WRAP, TP_CANDIDATE_MESH

    BW  = logical(BW);
    sz  = size(BW);
    if numel(sz) < 3; sz(end+1:3) = 1; end

    % Pad with a one-voxel background shell: every background voxel on the
    % image boundary is then face-connected to the shell, so seeding the
    % shell is equivalent to seeding the boundary background voxels, and
    % interior neighbour offsets can never leave the array or wrap.
    psz = sz + 2;
    bg  = true(psz);
    bg(2:end-1, 2:end-1, 2:end-1) = ~BW;

    visited = true(psz);
    visited(2:end-1, 2:end-1, 2:end-1) = false;
    frontier = find(visited);

    % Vectorised BFS with 6-connectivity (faces only) — holes use
    % 6-conn background. Each wavefront is expanded as one array operation,
    % so the loop runs once per BFS depth rather than once per voxel.
    offsets = [1, -1, psz(1), -psz(1), psz(1)*psz(2), -psz(1)*psz(2)];
    N = prod(psz);

    while ~isempty(frontier)
        nb = frontier(:) + offsets;          % [n x 6] linear neighbours
        nb = nb(:);
        nb = nb(nb >= 1 & nb <= N);          % shell seeds only
        nb = nb(bg(nb) & ~visited(nb));
        visited(nb) = true;
        frontier = unique(nb);
    end

    % Holes = unreachable background
    out = BW | ~visited(2:end-1, 2:end-1, 2:end-1);
end
