function [oa, planes] = read_offaxis_scans(src, options)
% READ_OFFAXIS_SCANS  Hydrophone plane and volume scans as off-axis data for FIT_RING_PROFILE
%
% Reads transverse (XY) and axial-plane (XZ / YZ) hydrophone scans and 3D
% volume scans, and reduces them to the radially symmetric samples an
% annular array model predicts: a transverse plane becomes its azimuthally
% averaged radial profile, an axial plane its (|x|, z) map with both sides
% folded together, and a volume its (r, z) map, azimuthally averaged per
% z bin around one axis parallel to z.
% The result is the .OffAxis input of FIT_RING_PROFILE (via the
% RingProfileOptions of CALIBRATION_FIT / CALIBRATE_LIBRARY_ENTRY).
%
% Each plane keeps its focal setting and absolute intensity: FIT_RING_PROFILE
% compares it with the axial peak at that setting, so the plane must be in
% the units and at the drive of the axial profiles. When the scans come with
% their own axial profiles (Sonic Concepts export) and options.Axial gives
% the profiles being fitted, each plane is rescaled by the ratio of the two
% axial peaks at its setting, so planes and axial profiles from different
% sessions combine.
%
% Formats:
%   Sonic Concepts "FUS Metrology" export (FUS_Metrology_Export_BableBrain.mat,
%     cell array .data with .type, .focus [mm from EP], .x/.y/.z [mm], .ISPPA.data
%     [W/cm²]): the .mat of the planes, its folder, the export folder holding
%     Planes/ and Axial/, or a .zip of it.
%   CSV with columns focus_ep_mm, x_mm, y_mm, z_ep_mm, intensity_wcm2 (one row
%     per point; a scan is the points sharing focus_ep_mm and, if present, a
%     label column). Points need not lie on a grid.
%   A scan's kind follows from its points: constant z = transverse plane,
%   constant x or y = axial plane, all three varying = volume.
%
% Use as:
%   oa = read_offaxis_scans(src)
%   [oa, planes] = read_offaxis_scans(src, options)
%
% Input:
%   src     - file or folder (see Formats)
%   options - struct, all optional:
%     .Kind            - 'all' (default) | 'transverse' (XY planes) | 'axial'
%                        (XZ / YZ planes) | 'volume' (3D scans); a cellstr
%                        selects several
%     .RMaxMm          - radial range used [mm] (default 30)
%     .BinMm           - radial bin of transverse planes and of axial planes [mm]
%                        (default 0.5)
%     .ZBinMm          - axial bin of axial planes and volumes [mm] (default 1;
%                        the field matrix grows with the number of points)
%     .Weight          - weight of each plane relative to one axial profile
%                        (scalar or one per plane; default 1)
%     .Centre          - 'origin' (default: the scan's axis is x = y = 0) |
%                        'centroid' (transverse planes and volumes: intensity
%                        centroid of the top 50 % of the scan; a volume's axis
%                        is taken parallel to z, tilt is not corrected)
%     .Axial           - axial profiles being fitted (READ_AXIAL_PROFILES) to
%                        rescale the planes to (default: no rescaling)
%     .NearFieldSkipMm - peak search beyond this [mm] for the rescaling (default 4)
%
% Output:
%   oa     - struct array for FIT_RING_PROFILE .OffAxis: .depth_ep_mm, .r_mm,
%            .z_ep_mm, .I [W/cm²], .weight, .label, .kind ('transverse' | 'axial' | 'volume'),
%            .scale (rescaling applied, 1 = none)
%   planes - per scan, for plotting: .label, .kind, .depth_ep_mm, .u, .v (grid
%            axes: x and y, x and z, or for a volume r and z), .I (map
%            [numel(v) x numel(u)], W/cm², rescaled; a volume's (r, z) map),
%            .z_ep_mm (transverse), .centre_mm
%
% See also: FIT_RING_PROFILE, FIELD_TRANSFER_MATRIX, READ_AXIAL_PROFILES, CALIBRATION_FIT

    if nargin < 2; options = struct(); end
    kind_sel = cellstr(opt(options, 'Kind', 'all'));
    if any(strcmp(kind_sel, 'all')); kind_sel = {'transverse', 'axial', 'volume'}; end
    assert(all(ismember(kind_sel, {'transverse', 'axial', 'volume'})), 'read_offaxis_scans:kind', ...
        'Kind must be ''all'', ''transverse'', ''axial'' or ''volume'' (or a list of them).');
    rmax  = opt(options, 'RMaxMm', 30);
    bin   = opt(options, 'BinMm', 0.5);
    zbin  = opt(options, 'ZBinMm', 1);
    centre_mode = char(opt(options, 'Centre', 'origin'));
    skip  = opt(options, 'NearFieldSkipMm', 4);

    [raw, ref] = read_points(char(src));
    keep = true(1, numel(raw));
    for p = 1:numel(raw)
        keep(p) = ismember(raw(p).kind, kind_sel);
    end
    raw = raw(keep);
    assert(~isempty(raw), 'read_offaxis_scans:none', 'No %s scans in %s.', strjoin(kind_sel, '/'), char(src));
    w = opt(options, 'Weight', 1);
    if isscalar(w); w = repmat(w, 1, numel(raw)); end
    assert(numel(w) == numel(raw), 'read_offaxis_scans:weight', ...
        'Weight needs one value or one per plane (%d planes).', numel(raw));

    meas = opt(options, 'Axial', []);
    oa = struct('depth_ep_mm', {}, 'r_mm', {}, 'z_ep_mm', {}, 'I', {}, 'weight', {}, ...
                'label', {}, 'kind', {}, 'scale', {});
    planes = struct('label', {}, 'kind', {}, 'depth_ep_mm', {}, 'u', {}, 'v', {}, 'I', {}, ...
                    'z_ep_mm', {}, 'centre_mm', {});
    for p = 1:numel(raw)
        P = raw(p);
        scale = rescale_factor(P.focus, ref, meas, skip);
        I = P.I * scale;
        c = [0 0];
        if strcmp(centre_mode, 'centroid') && ~strcmp(P.kind, 'axial')
            top = I >= 0.5 * max(I);
            c = [sum(P.x(top) .* I(top)), sum(P.y(top) .* I(top))] / sum(I(top));
        end
        switch P.kind
            case 'transverse'
                r  = hypot(P.x - c(1), P.y - c(2));
                rb = (0:bin:rmax)';
                Ir = arrayfun(@(q) mean(I(abs(r - q) <= bin / 2)), rb);
                ok = isfinite(Ir);
                z0 = median(P.z);
                oa(end + 1) = struct('depth_ep_mm', P.focus, 'r_mm', rb(ok), ...
                    'z_ep_mm', repmat(z0, nnz(ok), 1), 'I', Ir(ok), 'weight', w(p), ...
                    'label', P.label, 'kind', 'transverse', 'scale', scale); %#ok<AGROW>
                [u, v, M] = to_map(P.x, P.y, I);
                planes(end + 1) = struct('label', P.label, 'kind', 'transverse', 'depth_ep_mm', P.focus, ...
                    'u', u, 'v', v, 'I', M, 'z_ep_mm', z0, 'centre_mm', c); %#ok<AGROW>
            case 'axial'
                lat = P.x; if range(P.x) < range(P.y); lat = P.y; end
                r  = abs(lat);
                ir = round(r / bin); iz = round(P.z / zbin);
                ok = r <= rmax + bin / 2;
                [g, ~, idx] = unique([ir(ok), iz(ok)], 'rows');
                Ig = accumarray(idx, I(ok), [], @mean);
                oa(end + 1) = struct('depth_ep_mm', P.focus, 'r_mm', g(:, 1) * bin, ...
                    'z_ep_mm', g(:, 2) * zbin, 'I', Ig, 'weight', w(p), ...
                    'label', P.label, 'kind', 'axial', 'scale', scale); %#ok<AGROW>
                [u, v, M] = to_map(lat, P.z, I);
                planes(end + 1) = struct('label', P.label, 'kind', 'axial', 'depth_ep_mm', P.focus, ...
                    'u', u, 'v', v, 'I', M, 'z_ep_mm', NaN, 'centre_mm', [0 0]); %#ok<AGROW>
            case 'volume'                % azimuthal average per (r, z) bin around the axis
                r  = hypot(P.x - c(1), P.y - c(2));
                ir = round(r / bin); iz = round(P.z / zbin);
                ok = r <= rmax + bin / 2;
                [g, ~, idx] = unique([ir(ok), iz(ok)], 'rows');
                Ig = accumarray(idx, I(ok), [], @mean);
                oa(end + 1) = struct('depth_ep_mm', P.focus, 'r_mm', g(:, 1) * bin, ...
                    'z_ep_mm', g(:, 2) * zbin, 'I', Ig, 'weight', w(p), ...
                    'label', P.label, 'kind', 'volume', 'scale', scale); %#ok<AGROW>
                [u, v, M] = to_map(g(:, 1) * bin, g(:, 2) * zbin, Ig);
                planes(end + 1) = struct('label', P.label, 'kind', 'volume', 'depth_ep_mm', P.focus, ...
                    'u', u, 'v', v, 'I', M, 'z_ep_mm', NaN, 'centre_mm', c); %#ok<AGROW>
        end
    end
end

%% Readers -------------------------------------------------------------------------
function [raw, ref] = read_points(src)
% Points per plane (.label, .kind, .focus, .x, .y, .z, .I as columns) and the
% reference axial profiles that came with them (prof struct or [])
    ref = [];
    assert(isfile(src) || isfolder(src), 'read_offaxis_scans:missing', 'Not found: %s', src);
    [~, ~, ext] = fileparts(src);
    if strcmpi(ext, '.csv')
        raw = read_csv(src); return;
    end
    if strcmpi(ext, '.zip')
        tmp = tempname; unzip(src, tmp);
        c = onCleanup(@() rmdir(tmp, 's'));
        [raw, ref] = read_points(tmp); return;
    end
    if isfolder(src)
        d = dir(fullfile(src, '**', '*.mat'));
        assert(~isempty(d), 'read_offaxis_scans:noMat', 'No .mat export in %s.', src);
        raw = struct('label', {}, 'kind', {}, 'focus', {}, 'x', {}, 'y', {}, 'z', {}, 'I', {});
        for k = 1:numel(d)
            [r, a] = read_export(fullfile(d(k).folder, d(k).name));
            raw = [raw, r]; %#ok<AGROW>
            if ~isempty(a); ref = a; end
        end
        return;
    end
    [raw, ref] = read_export(src);
end

function [raw, ref] = read_export(file)
% Sonic Concepts FUS Metrology export: planes, and axial profiles as reference
    S = load(file);
    assert(isfield(S, 'data') && iscell(S.data), 'read_offaxis_scans:format', ...
        '%s is not a FUS Metrology export (no cell array ''data'').', file);
    raw = struct('label', {}, 'kind', {}, 'focus', {}, 'x', {}, 'y', {}, 'z', {}, 'I', {});
    ax = {};
    for k = 1:numel(S.data)
        v = S.data{k};
        x = double(v.x(:)); y = double(v.y(:)); z = double(v.z(:)); I = double(v.ISPPA.data(:));
        assert(strcmpi(char(v.ISPPA.unit), 'W/cm2'), 'read_offaxis_scans:unit', ...
            'Expected ISPPA in W/cm2 in %s, found %s.', file, char(v.ISPPA.unit));
        typ = char(v.type);
        if strcmpi(typ, 'Axial') || (range(x) == 0 && range(y) == 0)
            ax{end + 1} = struct('focus', double(v.focus), 'z', z, 'I', I); %#ok<AGROW>
            continue;
        end
        raw(end + 1) = struct('label', typ, 'kind', plane_kind(x, y, z), 'focus', double(v.focus), ...
            'x', x, 'y', y, 'z', z, 'I', I); %#ok<AGROW>
    end
    ref = [];
    if ~isempty(ax)
        z = ax{1}.z;
        same = cellfun(@(a) numel(a.z) == numel(z) && all(abs(a.z - z) < 1e-6), ax);
        assert(all(same), 'read_offaxis_scans:axialGrid', 'Axial profiles in %s differ in z.', file);
        ref = struct('depths_ep_mm', cellfun(@(a) a.focus, ax), 'z_ep_mm', z, ...
                     'I', cell2mat(cellfun(@(a) a.I, ax, 'UniformOutput', false)));
    end
end

function raw = read_csv(file)
    T = readtable(file, 'VariableNamingRule', 'preserve');
    need = {'focus_ep_mm', 'x_mm', 'y_mm', 'z_ep_mm', 'intensity_wcm2'};
    miss = setdiff(need, T.Properties.VariableNames);
    assert(isempty(miss), 'read_offaxis_scans:csvColumns', 'Plane CSV %s lacks column(s) %s.', ...
        file, strjoin(miss, ', '));
    if ismember('label', T.Properties.VariableNames)
        lab = string(T.label);
    else
        lab = strings(height(T), 1);
    end
    [g, ~, idx] = unique(table(T.focus_ep_mm, lab), 'rows');
    raw = struct('label', {}, 'kind', {}, 'focus', {}, 'x', {}, 'y', {}, 'z', {}, 'I', {});
    for k = 1:height(g)
        s = idx == k;
        x = T.x_mm(s); y = T.y_mm(s); z = T.z_ep_mm(s);
        kd = plane_kind(x, y, z);
        label = char(g{k, 2});
        if isempty(label); label = sprintf('%s %g mm setting', kd, g{k, 1}); end
        raw(end + 1) = struct('label', label, 'kind', kd, 'focus', g{k, 1}, ...
            'x', x, 'y', y, 'z', z, 'I', double(T.intensity_wcm2(s))); %#ok<AGROW>
    end
end

%% Helpers -------------------------------------------------------------------------
function k = plane_kind(x, y, z)
% Transverse: constant z. Axial plane: constant x or y (the axis lies in it).
% Volume: x, y and z all vary (several planes at one setting need a label
% column in a CSV, else they read as one volume)
    if range(z) < 1e-6
        k = 'transverse';
    elseif range(x) < 1e-6 || range(y) < 1e-6
        k = 'axial';
    else
        k = 'volume';
    end
end

function s = rescale_factor(focus, ref, meas, skip)
% Ratio of the fitted axial peak to the scans' own axial peak at this setting
    s = 1;
    if isempty(ref) || isempty(meas); return; end
    [d1, k1] = min(abs(ref.depths_ep_mm - focus));
    [d2, k2] = min(abs(meas.depths_ep_mm - focus));
    if d1 > 1 || d2 > 1
        warning('read_offaxis_scans:noReference', ['No axial profile within 1 mm of the %.1f mm ' ...
            'setting to rescale the plane to: used as is.'], focus);
        return;
    end
    p_ref  = max(ref.I(ref.z_ep_mm >= skip, k1));
    p_meas = max(meas.I(meas.z_ep_mm >= skip, k2));
    s = p_meas / p_ref;
end

function [u, v, M] = to_map(a, b, I)
    u = unique(a); v = unique(b);
    [~, iu] = ismember(a, u); [~, iv] = ismember(b, v);
    M = nan(numel(v), numel(u));
    M(sub2ind(size(M), iv, iu)) = I;
end

function v = opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
