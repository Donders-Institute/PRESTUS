function lib = library_add_refinements(lib, depths_ep_mm, offsets_deg, blend_mm, forward_model, source, amp)
% LIBRARY_ADD_REFINEMENTS  Merge depth-specific phase offsets (and drives) into a library entry
%
% Offsets are relative to geometric steering (at the entry's
% steering_sound_speed). Existing refinements at the same depth (±0.5 mm)
% are replaced, together with their drive amplitude and any k-Wave drive calibration
% result (those belonged to the old phases); others are kept.
%
% Use as:
%   lib = library_add_refinements(lib, depths_ep_mm, offsets_deg)
%   lib = library_add_refinements(lib, depths_ep_mm, offsets_deg, blend_mm, forward_model, source)
%   lib = library_add_refinements(..., source, amp)
%
% Input:
%   depths_ep_mm - [1 x K] refined depths [mm from EP]
%   offsets_deg  - {1 x K} cell of [1 x N_elem] offsets [deg]
%   blend_mm     - fade length back to the static model (default: keep, else 5)
%   amp          - (optional) struct with the drive that gives isppa_wcm2 in
%                  water at each refined depth with the refined phases:
%                    .elem_p_amp      [1 x K] mean source pressure [Pa]
%                    .isppa_wcm2      reference ISPPA [W/cm²]
%                    .elem_amp_weight {1 x K} per-ring gains (mean 1), if fitted
%                    .z_shift_mm      [1 x K] fitted axial registration [mm]
%                                     (diagnostic, not applied at run time)
%                    .profile_type    {1 x K} profile each depth was fitted to:
%                                     'measured' | 'interpolated' | 'extrapolated'
%                                     (INTERPOLATE_AXIAL_PROFILES)
%                  Without it the new depths carry no drive (NaN) and
%                  LOAD_TRANSDUCER_FROM_LIBRARY uses the global amplitude there.
%
% Per-depth fields written under 'refinement' (one value per depths_ep_mm):
%   elem_phase_correction_deg - phase offsets [deg]
%   elem_p_amp                - drive at refinement.isppa_wcm2 [Pa] (NaN = none)
%   amp_source                - 'analytical' | 'kwave' | '' (no drive)
%   elem_amp_weight           - per-ring gains [1 x N_elem] (NaN = global gains)
%   fit_z_shift_mm            - fitted axial registration [mm] (NaN = not fitted)
%   profile_type              - profile the depth was fitted to ('' = not recorded)
%   kwave_isppa_wcm2, kwave_r2, kwave_focal_err_mm
%                             - free-water k-Wave drive calibration (NaN = not run;
%                               written by KWAVE_CALIBRATE_DRIVE)
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, REFINE_TRANSDUCER_LIBRARY, KWAVE_CALIBRATE_DRIVE

arguments
    lib
    depths_ep_mm {mustBeNumeric}
    offsets_deg cell
    blend_mm = []
    forward_model = ''
    source = ''
    amp = []
end

    if nargin < 4; blend_mm = []; end
    if nargin < 5; forward_model = ''; end
    if nargin < 6; source = ''; end
    if nargin < 7; amp = []; end

    [lib, R] = refinement_table(lib);
    ref_isppa = R.isppa_wcm2;
    if ~isempty(amp)
        amp_new = double(amp.elem_p_amp(:)');
        if isnan(ref_isppa)
            ref_isppa = double(amp.isppa_wcm2);
        elseif abs(ref_isppa - amp.isppa_wcm2) > 1e-9
            amp_new = amp_new * sqrt(ref_isppa / double(amp.isppa_wcm2));   % I ∝ p²
        end
    else
        amp_new = nan(1, numel(depths_ep_mm));
    end
    gain_new = {}; dz_new = nan(1, numel(depths_ep_mm));
    if ~isempty(amp) && isfield(amp, 'elem_amp_weight'); gain_new = amp.elem_amp_weight; end
    if ~isempty(amp) && isfield(amp, 'z_shift_mm');      dz_new   = double(amp.z_shift_mm(:)'); end
    pt_new = repmat("", 1, numel(depths_ep_mm));
    if ~isempty(amp) && isfield(amp, 'profile_type');    pt_new   = string(amp.profile_type(:)'); end

    for k = 1:numel(depths_ep_mm)
        off = mod(double(offsets_deg{k}(:)'), 360);
        m = find(abs(R.d - depths_ep_mm(k)) < 0.5, 1);
        if isempty(m); m = numel(R.d) + 1; end
        R.d(m)    = depths_ep_mm(k);
        R.off{m}  = off;
        R.amp(m)  = amp_new(k);
        R.src(m)  = string(ternary(isnan(amp_new(k)), '', 'analytical'));
        if isempty(gain_new); R.gain{m} = nan(1, numel(off)); else; R.gain{m} = double(gain_new{k}(:)'); end
        R.dz(m)   = dz_new(k);
        R.pt(m)   = pt_new(k);
        R.kw_i(m) = NaN; R.kw_r2(m) = NaN; R.kw_fe(m) = NaN;
    end
    R.isppa_wcm2 = ref_isppa;
    lib = write_refinement_table(lib, R);

    if ~isempty(blend_mm)
        lib.refinement.blend_mm = blend_mm;
    elseif ~isfield(lib.refinement, 'blend_mm')
        lib.refinement.blend_mm = 5;
    end
    if ~isempty(forward_model); lib.refinement.forward_model = forward_model; end
    if ~isempty(source);        lib.refinement.source        = char(source); end
    lib.refinement.last_updated = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
end

%% ==========================================================================
function [lib, R] = refinement_table(lib)
% Read the per-depth refinement fields into aligned vectors (missing = NaN / '').
    R = struct('d', [], 'off', {{}}, 'amp', [], 'src', strings(1, 0), 'gain', {{}}, 'dz', [], 'pt', strings(1, 0), ...
               'kw_i', [], 'kw_r2', [], 'kw_fe', [], 'isppa_wcm2', NaN);
    if ~isfield(lib, 'refinement') || ~isstruct(lib.refinement) || ...
            ~isfield(lib.refinement, 'depths_ep_mm')
        lib.refinement = struct();
        return;
    end
    r = lib.refinement;
    R.d = double(r.depths_ep_mm(:)');
    n   = numel(R.d);
    c   = r.elem_phase_correction_deg;
    if ~iscell(c)
        if n == 1; c = {c}; else; c = num2cell(c, 2); end
    end
    R.off   = reshape(c, 1, []);
    R.amp   = field_vec(r, 'elem_p_amp', n);
    R.dz    = field_vec(r, 'fit_z_shift_mm', n);
    R.gain  = cellfun(@(o) nan(1, numel(o)), R.off, 'UniformOutput', false);
    if isfield(r, 'elem_amp_weight')
        g = r.elem_amp_weight;
        if ~iscell(g); if n == 1; g = {g}; else; g = num2cell(g, 2); end; end
        if numel(g) == n; R.gain = cellfun(@(x) double(x(:)'), reshape(g, 1, []), 'UniformOutput', false); end
    end
    R.kw_i  = field_vec(r, 'kwave_isppa_wcm2', n);
    R.kw_r2 = field_vec(r, 'kwave_r2', n);
    R.kw_fe = field_vec(r, 'kwave_focal_err_mm', n);
    if isfield(r, 'amp_source') && numel(r.amp_source) == n
        R.src = string(r.amp_source(:)');
    else
        R.src = repmat("", 1, n);
        R.src(isfinite(R.amp)) = "analytical";
    end
    R.pt = repmat("", 1, n);
    if isfield(r, 'profile_type') && numel(r.profile_type) == n; R.pt = string(r.profile_type(:)'); end
    if isfield(r, 'isppa_wcm2') && ~isempty(r.isppa_wcm2); R.isppa_wcm2 = double(r.isppa_wcm2); end
end

function lib = write_refinement_table(lib, R)
    [d, o] = sort(R.d);
    r = lib.refinement;
    r.depths_ep_mm              = d;
    r.elem_phase_correction_deg = R.off(o);
    if any(cellfun(@(x) any(isfinite(x)), R.gain))
        r.elem_amp_weight = R.gain(o);
    elseif isfield(r, 'elem_amp_weight')
        r = rmfield(r, 'elem_amp_weight');
    end
    if any(isfinite(R.dz)); r.fit_z_shift_mm = R.dz(o); end
    if any(R.pt ~= ""); r.profile_type = cellstr(R.pt(o)); end
    if any(isfinite(R.amp))
        r.elem_p_amp = round(R.amp(o));
        r.isppa_wcm2 = R.isppa_wcm2;
        r.amp_source = cellstr(R.src(o));
    end
    if any(isfinite(R.kw_i))
        r.kwave_isppa_wcm2   = R.kw_i(o);
        r.kwave_r2           = R.kw_r2(o);
        r.kwave_focal_err_mm = R.kw_fe(o);
    else
        for f = {'kwave_isppa_wcm2', 'kwave_r2', 'kwave_focal_err_mm'}
            if isfield(r, f{1}); r = rmfield(r, f{1}); end
        end
    end
    lib.refinement = r;
end

function v = field_vec(r, f, n)
    v = nan(1, n);
    if isfield(r, f) && numel(r.(f)) == n; v = double(r.(f)(:)'); end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end
