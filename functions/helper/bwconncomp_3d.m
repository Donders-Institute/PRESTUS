function CC = bwconncomp_3d(BW, conn)
% BWCONNCOMP_3D  Find connected components in a 3-D binary volume
%
% Drop-in replacement for bwconncomp(BW, conn) from the Image Processing
% Toolbox. Returns a struct with the same fields used by PRESTUS callers:
%   CC.NumObjects   - number of connected components
%   CC.PixelIdxList - cell array of linear index vectors (one per component)
%   CC.ImageSize    - size of BW
%   CC.Connectivity - conn (6 or 26)
%
% Use as:
%   CC = bwconncomp_3d(BW)          % 26-connectivity (default)
%   CC = bwconncomp_3d(BW, 26)
%   CC = bwconncomp_3d(BW, 6)
%
% Input:
%   BW   - logical 3-D array
%   conn - connectivity: 6 (face) or 26 (face+edge+vertex, default)
%
% See also: SKULL_RUBBER_WRAP

    if nargin < 2; conn = 26; end
    BW = logical(BW);
    sz = size(BW);
    if numel(sz) < 3; sz(end+1:3) = 1; end

    % Neighbour offsets for chosen connectivity
    [ox,oy,oz] = ndgrid(-1:1, -1:1, -1:1);
    ox = ox(:); oy = oy(:); oz = oz(:);
    self = ox==0 & oy==0 & oz==0;
    ox(self) = []; oy(self) = []; oz(self) = [];
    if conn == 6
        keep = abs(ox)+abs(oy)+abs(oz) == 1;
        ox = ox(keep); oy = oy(keep); oz = oz(keep);
    end

    labeled       = zeros(sz, 'uint32');
    n_obj         = uint32(0);
    PixelIdxList  = {};
    seeds         = uint32(find(BW));

    % BFS queue pre-allocated to worst-case volume size
    max_q   = uint32(numel(BW));
    queue   = zeros(max_q, 1, 'uint32');

    for si = 1:numel(seeds)
        seed = seeds(si);
        if labeled(seed); continue; end

        n_obj       = n_obj + uint32(1);
        head        = uint32(1);
        tail        = uint32(1);
        queue(1)    = seed;
        labeled(seed) = n_obj;

        while head <= tail
            curr = queue(head); head = head + uint32(1);
            [cx,cy,cz] = ind2sub(sz, double(curr));

            nx = cx + ox;
            ny = cy + oy;
            nz = cz + oz;
            valid = nx>=1 & nx<=sz(1) & ny>=1 & ny<=sz(2) & nz>=1 & nz<=sz(3);
            nx = nx(valid); ny = ny(valid); nz = nz(valid);
            ni = uint32(sub2ind(sz, nx, ny, nz));

            for k = 1:numel(ni)
                nb = ni(k);
                if BW(nb) && labeled(nb) == 0
                    labeled(nb)  = n_obj;
                    tail         = tail + uint32(1);
                    queue(tail)  = nb;
                end
            end
        end

        PixelIdxList{n_obj} = queue(1:tail); %#ok<AGROW>
    end

    CC.Connectivity  = conn;
    CC.ImageSize     = sz;
    CC.NumObjects    = double(n_obj);
    CC.PixelIdxList  = PixelIdxList;
end
