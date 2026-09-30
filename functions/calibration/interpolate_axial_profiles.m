function out = interpolate_axial_profiles(prof, depths_ep_mm, options)
% INTERPOLATE_AXIAL_PROFILES  Axial intensity profiles at desired focal depths from measured ones
%
% Builds a profile for each desired focal depth from hydrophone profiles
% measured at other depths (READ_AXIAL_PROFILES), so a transducer can be
% calibrated at a depth that was not measured:
%
%   measured depth      - the measured column itself                ('measured')
%   between two depths  - the two neighbouring profiles are shifted so their
%                         focal positions coincide at the interpolated
%                         position, then blended linearly
%                         (1-a)*I_1 + a*I_2, a = (d - d_1)/(d_2 - d_1)  ('interpolated')
%   outside the range   - only within MaxExtrapolationMm of the nearest
%                         measured depth: that profile, shifted by the
%                         depth difference                            ('extrapolated')
%   further outside     - NaN column                                 ('unavailable')
%
% The focal position used for the alignment is the column's focal-depth
% label (Align = 'label', default: the headers are the nominal setpoints,
% which correspond to the FLHM centre, see doc_transducer_library.md) or the
% intensity peak beyond the near field (Align = 'peak', as in the legacy
% EXTRACT_REAL_INTENSITY_PROFILE). The interpolated focal position is
% interpolated linearly between the two, and equals d for 'label'.
% Shifted profiles are resampled linearly; samples shifted in from beyond
% the scan take the nearest edge value. Near-field structure is shifted with
% the focus, so it is not physical: fits skip it (NearFieldSkipMm).
%
% Use as:
%   out = interpolate_axial_profiles(prof, depths_ep_mm)
%   out = interpolate_axial_profiles(prof, depths_ep_mm, options)
%
% Input:
%   prof         - measured profiles: .depths_ep_mm [1 x N_d], .z_ep_mm [N_z x 1],
%                  .I [N_z x N_d] (READ_AXIAL_PROFILES), optional .source
%   depths_ep_mm - desired focal depths [mm from EP]
%   options      - struct, all optional:
%     .Align              - 'label' (default) | 'peak'
%     .MaxExtrapolationMm - how far beyond the measured range a profile may be
%                           shifted from the nearest one [mm] (default 0: none)
%     .NearFieldSkipMm    - ignore z < this when locating the peak [mm] (default 10)
%     .Tolerance          - a desired depth within this of a measured one uses
%                           that measurement [mm] (default 0.01)
%
% Output (out): a profile struct usable wherever measured profiles are
% (FIT_TRANSDUCER_MODEL, WRITE_AXIAL_PROFILES, ...), with
%   .depths_ep_mm, .z_ep_mm, .I  - desired depths, measured z grid, profiles [W/cm²]
%   .profile_type   - {1 x N} 'measured' | 'interpolated' | 'extrapolated' | 'unavailable'
%   .source_depths  - [N x 2] measured depths each profile was built from (NaN = none)
%   .weights        - [N x 2] their weights
%   .shift_mm       - [N x 2] shift applied to each source profile [mm]
%   .align, .source - alignment used; source of the measured profiles
%
% See also: READ_AXIAL_PROFILES, WRITE_AXIAL_PROFILES, EXTRACT_REAL_INTENSITY_PROFILE,
%           REFINE_TRANSDUCER_LIBRARY, CALIBRATION_GUI

    if nargin < 3; options = struct(); end
    align   = opt(options, 'Align', 'label');
    max_ext = opt(options, 'MaxExtrapolationMm', 0);
    skip_mm = opt(options, 'NearFieldSkipMm', 10);
    tol     = opt(options, 'Tolerance', 0.01);
    if ~any(strcmp(align, {'label', 'peak'}))
        error('interpolate_axial_profiles: unknown Align ''%s'' (label | peak).', align);
    end

    dm = prof.depths_ep_mm(:)';
    z  = prof.z_ep_mm(:);
    Im = prof.I;
    [dm, o] = sort(dm); Im = Im(:, o);
    Im(~isfinite(Im) | Im < 0) = 0;
    switch align
        case 'label'
            zf = dm;
        case 'peak'
            far = z >= skip_mm;
            zf  = arrayfun(@(k) peak_at(z(far), Im(far, k)), 1:numel(dm));
    end

    dq = depths_ep_mm(:)';
    n  = numel(dq);
    I  = nan(numel(z), n);
    type = repmat({'unavailable'}, 1, n);
    src = nan(n, 2); w = nan(n, 2); sh = nan(n, 2);
    for k = 1:n
        d = dq(k);
        [dist, j] = min(abs(dm - d));
        if dist <= tol
            I(:, k) = Im(:, j);
            type{k} = 'measured'; src(k, 1) = dm(j); w(k, 1) = 1; sh(k, 1) = 0;
        elseif d > dm(1) && d < dm(end)
            j2 = find(dm > d, 1); j1 = j2 - 1;
            a  = (d - dm(j1)) / (dm(j2) - dm(j1));
            zq = (1 - a) * zf(j1) + a * zf(j2);            % interpolated focal position
            s1 = zq - zf(j1); s2 = zq - zf(j2);
            I(:, k) = (1 - a) * shifted(z, Im(:, j1), s1) + a * shifted(z, Im(:, j2), s2);
            type{k} = 'interpolated';
            src(k, :) = dm([j1 j2]); w(k, :) = [1 - a, a]; sh(k, :) = [s1, s2];
        elseif dist <= max_ext
            s = d - dm(j);
            I(:, k) = shifted(z, Im(:, j), s);
            type{k} = 'extrapolated'; src(k, 1) = dm(j); w(k, 1) = 1; sh(k, 1) = s;
        end
    end

    out = struct('depths_ep_mm', dq, 'z_ep_mm', z, 'I', I);
    out.profile_type  = type;
    out.source_depths = src;
    out.weights       = w;
    out.shift_mm      = sh;
    out.align         = align;
    out.source        = '';
    if isfield(prof, 'source'); out.source = char(prof.source); end
end

%% ======================================================================
function y = shifted(z, I, s)
% I moved by s along z (s > 0: deeper), on the same grid
    y = interp1(z, I, z - s, 'linear', NaN);
    y = fillmissing(y, 'nearest');
    y(y < 0) = 0;
end

function zp = peak_at(z, I)
    [~, k] = max(I);
    zp = z(k);
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
