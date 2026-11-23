function [rp, info] = fit_ring_profile(tran, prof, options)
% FIT_RING_PROFILE  Fit a within-ring velocity profile shared by all focal depths
%
% The uniform-ring model assumes every ring vibrates evenly across its width.
% For some transducers (notably several Sonic Concepts CTX) it cannot
% reproduce the measured axial profiles, in particular the ratio of the
% near-field peak to the focal peak, whatever phases and gains are chosen.
% This function fits how each ring vibrates across its width: every ring is
% split into NSub sub-rings with a complex weight that is the SAME at every
% focal depth (a hardware property), while each depth keeps only what the
% driver controls, one phase and gain per channel:
%
%   I_d(z) = | sum_e H_e(z) * u_e * c_d(channel(e)) |²
%
% u is fitted jointly over all depths with a smoothness penalty between
% neighbouring sub-rings (Smooth) and normalised so that each ring's
% area-weighted mean weight is 1 (its level and phase stay with the channel
% drive). The fit alternates between per-depth channel drives (u fixed) and
% u (drives fixed). Peak-normalised intensity, near field included from
% NearFieldSkipMm (default 4 mm).
%
% The profile is a model of the axial data: profiles measured on axis do not
% determine the surface distribution uniquely, so the off-axis field it
% implies is not validated by this fit (validate against lateral / 2D scans).
%
% Use as:
%   rp = fit_ring_profile(tran, prof)
%   [rp, info] = fit_ring_profile(tran, prof, options)
%
% Input:
%   tran    - annular transducer struct from load_equipment_config (equipment
%             geometry, without a profile)
%   prof    - measured profiles (READ_AXIAL_PROFILES): .depths_ep_mm, .z_ep_mm, .I
%   options - struct, all optional:
%     .NSub            - sub-rings per ring (default 4)
%     .Smooth          - weight of the smoothness penalty (default 0.3; 0 = none)
%     .NearFieldSkipMm - ignore z < this [mm from EP] (default 4)
%     .Iterations      - alternations (default 6)
%     .NStarts         - starts per depth for the channel drives (default 12)
%     .AmpBounds       - bounds on |u| and on channel gains (default [0.1 10])
%     .HeldOut         - true: also score each depth with u fitted on the
%                        other depths (leave-one-depth-out; default false)
%     .ForwardModel    - 'rayleigh' (default) | 'oneil'
%     .Profile         - a ring profile struct to score instead of fitting one
%                        (only the per-depth channel drives are fitted)
%     .OffAxis         - off-axis measurements fitted with the axial profiles
%                        (FIELD_TRANSFER_MATRIX), struct array with
%                        .depth_ep_mm (focal setting; must match a profile
%                        depth within 1 mm), .r_mm, .z_ep_mm (points: radial
%                        distance from the axis, distance from the exit
%                        plane), .I (intensity in the units of prof.I: it is
%                        normalised by that setting's axial peak, so the
%                        level relative to the axial scan counts), .weight
%                        (default 1). E.g. azimuthally averaged radial
%                        profiles of transverse hydrophone planes. Or the
%                        path of a plane-scan file, read with
%                        READ_OFFAXIS_SCANS (rescaled to prof's axial peaks)
%     .OffAxisOptions  - READ_OFFAXIS_SCANS options for a path (Kind, Weight,
%                        RMaxMm, BinMm, ZBinMm, Centre)
%     .Seed            - RNG seed (default 1)
%
% Output:
%   rp   - ring profile for global_model.ring_profile: .n_sub, .weight_abs,
%          .weight_phase_deg ([1 x rings*n_sub], inner to outer within each
%          ring), .smooth, .near_field_skip_mm, .fit_depths_ep_mm,
%          .r2_uniform_median, .r2_profile_median, .r2_heldout_median (NaN
%          unless HeldOut), .forward_model; with off-axis data also
%          .offaxis_source (file, '' for a struct), .offaxis_planes (labels)
%          and .offaxis_r2 (R² per plane)
%   info - .r2_uniform, .r2_profile, .r2_heldout [1 x N_d], .z_ep_mm,
%          .I_meas, .I_uniform, .I_profile [N_z x N_d] (peak-normalised, z >=
%          NearFieldSkipMm), .depths_ep_mm, .drives / .drives_uniform
%          [N_channels x N_d] complex channel drives (phase of channel 1 = 0),
%          .offaxis (per OffAxis entry: .r_mm, .z_ep_mm, .I_meas, .I_uniform,
%          .I_profile normalised by the axial peak, .r2_uniform, .r2_profile)
%
% See also: RING_PROFILE_TRANSDUCER, LIBRARY_TRANSDUCER, PLOT_RING_PROFILE,
%           CALIBRATE_LIBRARY_ENTRY, FIT_TRANSDUCER_MODEL

arguments
    tran (1,1) struct
    prof (1,1) struct
    options (1,1) struct = struct()
end

    if nargin < 3; options = struct(); end
    fixed   = opt(options, 'Profile', []);
    M       = opt(options, 'NSub', 4);
    if ~isempty(fixed); M = double(fixed.n_sub); end
    smooth  = opt(options, 'Smooth', 0.3);
    skip    = opt(options, 'NearFieldSkipMm', 4);
    n_it    = opt(options, 'Iterations', 6);
    n_st    = opt(options, 'NStarts', 12);
    ab      = opt(options, 'AmpBounds', [0.1 10]);
    fwd     = opt(options, 'ForwardModel', 'rayleigh');
    rng(opt(options, 'Seed', 1));

    assert(~(isfield(tran, 'ring_profile_applied') && tran.ring_profile_applied), ...
        'fit_ring_profile:applied', 'Pass the equipment geometry (without a ring profile).');
    a   = tran.transducer.annular;
    nr  = double(a.elem_n);
    n_ch = nr;
    if isfield(tran, 'n_elem') && ~isempty(tran.n_elem) && tran.n_elem < nr && mod(nr, tran.n_elem) == 0
        n_ch = double(tran.n_elem);
    end
    ring_of = kron(1:nr, ones(1, M));                    % sub-ring -> ring
    chan_of = kron(1:n_ch, ones(1, nr * M / n_ch));      % sub-ring -> channel

    % Uniform sub-ring transfer matrix on the fitted z range
    tr_u = ring_profile_transducer(tran, struct('n_sub', M));
    off  = annular_bowl_depth_mm(a);
    z    = prof.z_ep_mm(:); sel = z >= skip;
    H    = axial_transfer_matrix(rmfield_ann(tr_u, 'elem_weight'), z(sel) + off, 1500, 994, fwd);
    area = (tr_u.transducer.annular.elem_od_mm.^2 - tr_u.transducer.annular.elem_id_mm.^2);
    Y    = prof.I(sel, :);
    keep = any(isfinite(Y), 1);
    Y    = Y(:, keep); depths = prof.depths_ep_mm(keep);
    ymax = max(Y, [], 1);
    Y    = Y ./ ymax;
    nd   = size(Y, 2);
    S    = struct('H', H, 'chan', chan_of, 'ring', ring_of, 'n_ch', n_ch, 'nr', nr, 'M', M, ...
                  'area', area, 'n_st', n_st, 'ab', ab, 'smooth', smooth);
    % Off-axis measurements, attached to the depth of their focal setting
    S.X = cell(1, nd);
    OA  = opt(options, 'OffAxis', []);
    oa_src = '';
    if ischar(OA) || isstring(OA)                        % plane-scan file: READ_OFFAXIS_SCANS
        oa_src = char(OA);
        ro = opt(options, 'OffAxisOptions', struct());
        if ~isfield(ro, 'Axial'); ro.Axial = prof; end
        if ~isfield(ro, 'NearFieldSkipMm'); ro.NearFieldSkipMm = skip; end
        OA = read_offaxis_scans(oa_src, ro);
    end
    tr_h = rmfield_ann(tr_u, 'elem_weight');
    for j = 1:numel(OA)
        [dmin, k] = min(abs(depths - OA(j).depth_ep_mm));
        assert(dmin <= 1, 'fit_ring_profile:offAxisDepth', ...
            'OffAxis(%d) is at the %.1f mm setting; no axial profile within 1 mm.', j, OA(j).depth_ep_mm);
        Hx = field_transfer_matrix(tr_h, OA(j).r_mm(:), OA(j).z_ep_mm(:) + off, 1500, 994);
        x  = struct('H', Hx, 'y', double(OA(j).I(:)) / ymax(k), 'weight', field_or(OA(j), 'weight', 1), ...
                    'r_mm', OA(j).r_mm(:), 'z_ep_mm', OA(j).z_ep_mm(:), 'entry', j, ...
                    'label', char(field_or(OA(j), 'label', sprintf('off-axis %d', j))));
        S.X{k} = [S.X{k}, x];
    end

    % Uniform rings (reference) and the joint fit
    u1 = ones(1, nr * M);
    [C0, r2u, Iu] = fit_drives(S, Y, u1);
    if isempty(fixed)
        [u, C, r2p, Ip] = fit_joint(S, Y, C0, n_it);
    else
        u = double(fixed.weight_abs(:)') .* exp(1i * double(fixed.weight_phase_deg(:)') * pi / 180);
        [C, r2p, Ip] = fit_drives(S, Y, u);
        smooth = field_or(fixed, 'smooth', NaN);
    end
    r2h = nan(1, nd);
    if opt(options, 'HeldOut', false) && nd >= 3
        for k = 1:nd
            o = setdiff(1:nd, k);
            uk = fit_joint(sub(S, o), Y(:, o), C0(:, o), n_it);
            [~, r2h(k)] = fit_drives(sub(S, k), Y(:, k), uk);
        end
    end

    rp = struct('n_sub', M, 'weight_abs', round(abs(u), 4), ...
        'weight_phase_deg', round(angle(u) * 180 / pi, 2), 'smooth', smooth, ...
        'near_field_skip_mm', skip, 'fit_depths_ep_mm', depths(:)', ...
        'r2_uniform_median', median(r2u), 'r2_profile_median', median(r2p), ...
        'r2_heldout_median', median(r2h), 'forward_model', fwd);
    info = struct('r2_uniform', r2u, 'r2_profile', r2p, 'r2_heldout', r2h, 'z_ep_mm', z(sel), ...
        'I_meas', Y, 'I_uniform', Iu, 'I_profile', Ip, 'depths_ep_mm', depths(:)', ...
        'drives', C, 'drives_uniform', C0, 'offaxis', offaxis_report(S, u, u1, C, C0));
    if ~isempty(info.offaxis)                            % provenance of the off-axis data
        rp.offaxis_source = oa_src;
        rp.offaxis_planes = strjoin({info.offaxis.label}, '; ');
        rp.offaxis_r2 = round([info.offaxis.r2_profile], 4);
    end
end

function S = sub(S, cols)
% S restricted to some depths (columns of Y)
    S.X = S.X(cols);
end

function out = offaxis_report(S, u, u1, C, C0)
    out = struct('r_mm', {}, 'z_ep_mm', {}, 'I_meas', {}, 'I_uniform', {}, 'I_profile', {}, ...
        'r2_uniform', {}, 'r2_profile', {}, 'depth_index', {}, 'label', {});
    for k = 1:numel(S.X)
        if isempty(S.X{k}); continue; end
        [~, ~, Iu] = depth_res(S, k, u1, C0(:, k), []);
        [~, ~, Ip] = depth_res(S, k, u, C(:, k), []);
        for b = 1:numel(S.X{k})
            x = S.X{k}(b);
            ss = sum((x.y - mean(x.y)).^2);
            out(end + 1) = struct('r_mm', x.r_mm, 'z_ep_mm', x.z_ep_mm, 'I_meas', x.y, 'I_uniform', Iu{b}, ...
                'I_profile', Ip{b}, 'r2_uniform', 1 - sum((x.y - Iu{b}).^2) / ss, ...
                'r2_profile', 1 - sum((x.y - Ip{b}).^2) / ss, 'depth_index', k, 'label', x.label); %#ok<AGROW>
        end
    end
end

%% ---------------------------------------------------------------------------
function I = model_I(S, u, c)
% Peak-normalised axial intensity for sub-ring weights u and channel drive c
    w = u(:) .* c(S.chan(:));
    I = abs(S.H * w).^2;
    I = I / max(I);
end

function [r, Ia, Ix] = depth_res(S, k, u, c, y)
% Residual of depth k: axial (peak-normalised) and off-axis blocks (normalised
% by the same axial model peak, so their level relative to the axis counts)
    w  = u(:) .* c(S.chan(:));
    pa = abs(S.H * w).^2; m = max(pa);
    Ia = pa / m;
    if isempty(y); r = []; else; r = Ia - y; end
    Ix = {};
    if ~isempty(S.X) && numel(S.X) >= k && ~isempty(S.X{k})
        for x = S.X{k}
            ix = abs(x.H * w).^2 / m;
            Ix{end + 1} = ix; %#ok<AGROW>
            if ~isempty(y); r = [r; sqrt(x.weight / numel(x.y) * numel(y)) * (ix - x.y)]; end %#ok<AGROW>
        end
    end
end

function [C, r2, I] = fit_drives(S, Y, u)
% Per-depth channel drive (phase of channel 1 fixed) for fixed u
    nd = size(Y, 2); n = S.n_ch;
    C = complex(zeros(n, nd)); r2 = nan(1, nd); I = nan(size(Y));
    lb = [-inf(1, n - 1), log(S.ab(1)) * ones(1, n)]; ub = [inf(1, n - 1), log(S.ab(2)) * ones(1, n)];
    opts = optimoptions('lsqnonlin', 'Display', 'off');
    unpack = @(x) exp(x(n:end))' .* exp(1i * [0, x(1:n - 1)])';
    for k = 1:nd
        y = Y(:, k); best = inf; xb = [];
        f = @(x) depth_res(S, k, u, unpack(x), y);
        for s = 1:S.n_st
            x0 = [2 * pi * rand(1, n - 1), 0.2 * randn(1, n)];
            [x, rn] = lsqnonlin(f, x0, lb, ub, opts);
            if rn < best; best = rn; xb = x; end
        end
        C(:, k) = unpack(xb);
        I(:, k) = model_I(S, u, C(:, k));
        r2(k) = 1 - sum((y - I(:, k)).^2) / sum((y - mean(y)).^2);   % axial only
    end
end

function [u, C, r2, I] = fit_joint(S, Y, C, n_it)
% Alternate: shared u given the drives, drives given u
    u = ones(1, S.nr * S.M);
    for it = 1:n_it
        u = fit_u(S, Y, C, u);
        [C, r2, I] = fit_drives(S, Y, u);
    end
end

function u = fit_u(S, Y, C, u0)
% Shared sub-ring weights; smoothness between neighbouring sub-rings of a ring
    ne = numel(u0); nd = size(Y, 2);
    unpack = @(x) exp(x(ne + 1:end)) .* exp(1i * x(1:ne));
    same_ring = S.ring(1:end - 1) == S.ring(2:end);
    rough = @(uu) sqrt(S.smooth) * [real(diff(uu)); imag(diff(uu))];
    function r = res(x)
        uu = unpack(x);
        r = cell(nd, 1);
        for k = 1:nd
            r{k} = depth_res(S, k, uu, C(:, k), Y(:, k));
        end
        r = vertcat(r{:});
        if S.smooth > 0
            d = rough(uu);
            r = [r; reshape(d(:, same_ring), [], 1)];
        end
    end
    x0 = [angle(u0), log(abs(u0))];
    lb = [-inf(1, ne), log(S.ab(1)) * ones(1, ne)]; ub = [inf(1, ne), log(S.ab(2)) * ones(1, ne)];
    x  = lsqnonlin(@res, x0, lb, ub, optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 200));
    u  = unpack(x);
    % Normalise each ring: area-weighted mean weight = 1 (level and phase stay with the drive)
    for i = 1:S.nr
        e = S.ring == i;
        u(e) = u(e) / (sum(S.area(e) .* u(e)) / sum(S.area(e)));
    end
end

function tran = rmfield_ann(tran, f)
    if isfield(tran.transducer.annular, f); tran.transducer.annular = rmfield(tran.transducer.annular, f); end
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function v = field_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
