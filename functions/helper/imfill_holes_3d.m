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
    N   = prod(sz);

    bg  = ~BW;

    % Seed: all background voxels on the 6-face boundary
    face_mask = false(sz);
    face_mask(1,:,:)   = true; face_mask(end,:,:) = true;
    face_mask(:,1,:)   = true; face_mask(:,end,:) = true;
    face_mask(:,:,1)   = true; face_mask(:,:,end) = true;

    seeds = uint32(find(bg & face_mask));
    if isempty(seeds)
        out = true(sz);
        return
    end

    visited       = false(sz);
    visited(seeds) = true;

    % BFS with 6-connectivity (faces only) — holes use 6-conn background
    queue  = zeros(N, 1, 'uint32');
    queue(1:numel(seeds)) = seeds;
    head   = uint32(1);
    tail   = uint32(numel(seeds));

    offsets_d  = [1 -1 0  0  0  0];
    offsets_ri = [0  0 1 -1  0  0];
    offsets_c  = [0  0 0  0  1 -1];

    while head <= tail
        curr = queue(head); head = head + uint32(1);
        [cx,cy,cz] = ind2sub(sz, double(curr));

        for d = 1:6
            nx = cx + offsets_d(d);
            ny = cy + offsets_ri(d);
            nz = cz + offsets_c(d);
            if nx<1||nx>sz(1)||ny<1||ny>sz(2)||nz<1||nz>sz(3); continue; end
            ni = uint32(sub2ind(sz, nx, ny, nz));
            if bg(ni) && ~visited(ni)
                visited(ni)     = true;
                tail            = tail + uint32(1);
                queue(tail)     = ni;
            end
        end
    end

    % Holes = unreachable background
    out = BW | ~visited;
end
