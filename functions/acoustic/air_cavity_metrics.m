function m = air_cavity_metrics(pressure, aircav_mask, p_ref, beams, resolution_mm, warn_db, warn_beam_fraction)
% AIR_CAVITY_METRICS  How much of the beam enters internal air cavities
%
% Read-out for the air-cavity beam-path check (headmodel.air_cavities). The
% field is simulated with the cavities filled as their layer (water by
% default), so the pressure inside the cavity mask is the beam that would
% meet the air; real air reflects nearly all of it.
%
% Two measures:
%   rel_db    - peak pressure inside the cavities relative to p_ref (the
%               brain peak) [dB]: 20*log10(Psptp_aircav / p_ref)
%   beam_frac - per beam, along the axis from the exit plane to the focus in
%               1-voxel slabs, restricted to a cylinder of the aperture
%               radius: the -6 dB beam cross-section of a slab is its voxels
%               with p >= half the slab's peak; beam_frac is the share of it
%               inside the cavities. Slabs with fewer than 5 beam voxels are
%               ignored. The largest share over all slabs (and beams) is
%               reported, with its depth from the exit plane, and the share
%               over the whole prefocal -6 dB volume.
% The run is flagged when rel_db > warn_db or beam_frac_max > warn_beam_fraction.
%
% Use as:
%   m = air_cavity_metrics(pressure, aircav_mask, p_ref, beams, resolution_mm, ...
%                          warn_db, warn_beam_fraction)
%
% Input:
%   pressure           - temporal peak pressure map [Pa] (2D or 3D)
%   aircav_mask        - logical cavity mask, same size
%   p_ref              - reference pressure [Pa] (brain peak)
%   beams              - struct array, one per transducer, with fields
%                        trans_pos [1xD voxels], focus_pos [1xD voxels],
%                        exit_offset_vox (apex to exit plane along the axis [voxels]),
%                        radius_vox (aperture radius [voxels]; Inf = no lateral limit)
%   resolution_mm      - grid spacing [mm]
%   warn_db            - flag threshold on rel_db [dB]
%   warn_beam_fraction - flag threshold on beam_frac_max [0-1]
%
% Output:
%   m - struct: n_vox, Psptp_aircav [Pa], rel_db [dB], beam_frac_max [0-1],
%       beam_frac_depth_mm (from the exit plane), beam_frac_vol [0-1], flag (0/1)
%
% See also: AIR_CAVITY_MASK, ACOUSTIC_ANALYSIS

arguments
    pressure      {mustBeNumeric}
    aircav_mask   {mustBeNumericOrLogical}
    p_ref         (1,1) double
    beams         struct
    resolution_mm (1,1) double {mustBePositive}
    warn_db       (1,1) double
    warn_beam_fraction (1,1) double
end

    min_slab_vox = 5;
    aircav_mask = logical(aircav_mask);
    if ~isequal(size(aircav_mask), size(pressure))
        error('air_cavity_metrics:size', 'Cavity mask [%s] and pressure map [%s] differ in size.', ...
            num2str(size(aircav_mask)), num2str(size(pressure)));
    end

    m = struct('n_vox', nnz(aircav_mask), 'Psptp_aircav', NaN, 'rel_db', NaN, ...
               'beam_frac_max', 0, 'beam_frac_depth_mm', NaN, 'beam_frac_vol', 0, 'flag', 0);
    if m.n_vox == 0
        return
    end

    pressure = double(gather(pressure));
    m.Psptp_aircav = max(pressure(aircav_mask));
    if p_ref > 0
        m.rel_db = 20 * log10(m.Psptp_aircav / p_ref);
    end

    sz = size(pressure);
    nd = numel(sz);
    nb_all = 0; nc_all = 0;
    for b = 1:numel(beams)
        tp = double(beams(b).trans_pos(1:nd));
        fp = double(beams(b).focus_pos(1:nd));
        L = norm(fp - tp);
        if L == 0, continue; end
        u = (fp - tp) / L;
        e0 = double(beams(b).exit_offset_vox);
        R  = double(beams(b).radius_vox);
        % Work in the bounding box of the beam cylinder (exit plane to focus,
        % radius R), not the whole grid
        p0 = tp + u * e0;
        lo = max(1, floor(min(p0, fp) - min(R, max(sz))));
        hi = min(sz, ceil(max(p0, fp) + min(R, max(sz))));
        if any(hi < lo), continue; end
        bsz = hi - lo + 1;
        % Axial coordinate t and squared lateral distance r2 of every box voxel
        % [voxels], by implicit expansion over the grid axes (single precision)
        t = zeros(bsz, 'single');
        sq = zeros(bsz, 'single');
        for d = 1:nd
            shape = ones(1, max(nd, 2)); shape(d) = bsz(d);
            g = reshape(single((lo(d):hi(d)) - tp(d)), shape);
            t = t + u(d) * g;
            sq = sq + g.^2;
        end
        r2 = sq - t.^2;
        clear sq
        bidx = find(t >= e0 & t <= L & r2 <= R^2);
        if isempty(bidx), clear t r2; continue; end
        k  = floor(double(t(bidx)) - e0) + 1;         % 1-voxel slabs from the exit plane
        clear t r2
        % Box indices -> grid indices
        sub = cell(1, nd);
        [sub{:}] = ind2sub(bsz, bidx);
        for d = 1:nd, sub{d} = sub{d} + lo(d) - 1; end
        idx = sub2ind(sz, sub{:});
        clear sub bidx
        pv = pressure(idx);
        K  = max(k);
        slab_max = accumarray(k, pv, [K 1], @max, 0);
        in_beam  = pv >= 0.5 * slab_max(k) & pv > 0;
        in_cav   = aircav_mask(idx);
        nb = accumarray(k(in_beam), 1, [K 1]);
        nc = accumarray(k(in_beam & in_cav), 1, [K 1]);
        frac = zeros(K, 1);
        ok = nb >= min_slab_vox;
        frac(ok) = nc(ok) ./ nb(ok);
        [fmax, kmax] = max(frac);
        if fmax > m.beam_frac_max
            m.beam_frac_max = fmax;
            m.beam_frac_depth_mm = (kmax - 0.5) * resolution_mm;
        end
        nb_all = nb_all + sum(nb);
        nc_all = nc_all + sum(nc);
    end
    if nb_all > 0
        m.beam_frac_vol = nc_all / nb_all;
    end
    m.flag = double((isfinite(m.rel_db) && m.rel_db > warn_db) || m.beam_frac_max > warn_beam_fraction);
end
