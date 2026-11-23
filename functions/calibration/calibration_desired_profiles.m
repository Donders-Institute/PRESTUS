function p = calibration_desired_profiles(tran, meas, depths, options)
% CALIBRATION_DESIRED_PROFILES  Profiles on a desired-depth grid, for prediction and fitting
%
% The profile set behind "Predict" in CALIBRATION_GUI and the Calibration
% page of the PRESTUS Nexus. Each desired depth gets its measured profile,
% one interpolated between the neighbouring measured depths, or, within
% MaxExtrapolationMm, the nearest one shifted (INTERPOLATE_AXIAL_PROFILES).
% Without measurements every depth is 'unavailable' (model only). The z axis
% runs to 1.5 x the deepest of the desired depths and the geometric focus, so
% the model can be shown beyond the scan; the profiles are NaN there.
%
% Use as:
%   p = calibration_desired_profiles(tran, meas, depths)
%   p = calibration_desired_profiles(tran, meas, depths, options)
%
% Input:
%   tran    - annular transducer struct (load_equipment_config)
%   meas    - measured profiles (READ_AXIAL_PROFILES), or [] for none
%   depths  - desired depths [mm from EP]; rounded to 0.01 mm and sorted
%   options - struct, all optional (INTERPOLATE_AXIAL_PROFILES):
%     .Align              - 'label' (default) | 'peak'
%     .MaxExtrapolationMm - default 0 (no extrapolation)
%     .NearFieldSkipMm    - default 10
%     .Tolerance          - a desired depth within this of a measured one uses
%                           the measurement [mm] (default 0.01)
%
% Output:
%   p - .depths_ep_mm [1 x N], .z_ep_mm, .I [N_z x N] (NaN where there is no
%       profile), .profile_type {1 x N}: 'measured' | 'interpolated' |
%       'extrapolated' | 'unavailable'; with meas also the provenance fields
%       of INTERPOLATE_AXIAL_PROFILES
%
% See also: INTERPOLATE_AXIAL_PROFILES, CALIBRATION_EVALUATE, CALIBRATION_FIT

arguments
    tran (1,1) struct
    meas
    depths {mustBeNumeric}
    options (1,1) struct = struct()
end

    if nargin < 4; options = struct(); end
    dep  = unique(round(double(depths(:)'), 2));
    zmax = ceil(1.5 * max([dep, tran.transducer.annular.dist_geom_ep_mm]));
    if isempty(meas)
        z = (0:0.25:zmax)';
        p = struct('depths_ep_mm', dep, 'z_ep_mm', z, 'I', nan(numel(z), numel(dep)), 'source', '');
        p.profile_type = repmat({'unavailable'}, 1, numel(dep));
        return;
    end
    io = struct('Align', opt(options, 'Align', 'label'), ...
        'MaxExtrapolationMm', opt(options, 'MaxExtrapolationMm', 0), ...
        'NearFieldSkipMm', opt(options, 'NearFieldSkipMm', 10), ...
        'Tolerance', opt(options, 'Tolerance', 0.01));
    p  = interpolate_axial_profiles(meas, dep, io);
    z  = p.z_ep_mm; dz = median(diff(z));
    p.n_z_scan = numel(z);
    if z(end) < zmax                              % model beyond the scan; no profile there
        zx = (z(end) + dz:dz:zmax)';
        p.z_ep_mm = [z; zx]; p.I = [p.I; nan(numel(zx), numel(dep))];
    end
end

function v = opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
