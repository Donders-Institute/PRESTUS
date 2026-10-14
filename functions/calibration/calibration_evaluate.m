function [res, msg] = calibration_evaluate(tran, prof, options)
% CALIBRATION_EVALUATE  Geometric, library, fitted and k-Wave profiles on one set of depths
%
% The comparison view of CALIBRATION_GUI and the Calibration page of the
% PRESTUS Nexus. At the measured depths (Predicted false) it shows how well
% each model reproduces the hydrophone profiles (R², FIT_TRANSDUCER_MODEL
% definitions). At desired depths (Predicted true, CALIBRATION_DESIRED_PROFILES)
% it resolves the library entry and the unsaved fit as a study would
% (EVALUATE_LIBRARY_ENTRY: refinements interpolated between refined depths,
% blended back to the static model outside them) and reports the focus
% (FLHM centre beyond the near field) and R² against the measured or derived
% profile, with the full modelled profile.
%
% Use as:
%   res = calibration_evaluate(tran, prof)
%   [res, msg] = calibration_evaluate(tran, prof, options)
%
% Input:
%   tran    - annular transducer struct (load_equipment_config)
%   prof    - profiles on view: measured (READ_AXIAL_PROFILES with
%             profile_type) or desired (CALIBRATION_DESIRED_PROFILES)
%   options - struct, all optional:
%     .Predicted      - prof is a desired-depth grid (default false)
%     .Library        - current library entry, or []
%     .Fit            - CALIBRATION_FIT result, or []
%     .FitEntry       - entry the fit resolves to (default Fit.entry; the
%                       k-Wave-corrected entry after a drive check)
%     .Kwave          - k-Wave drive check (.report, .sims; KWAVE_CALIBRATE_DRIVE)
%     .FitOptions     - FIT_TRANSDUCER_MODEL options (ForwardModel, NearFieldSkipMm, ...)
%     .IsppaWcm2      - reference ISPPA for the runtime settings (default 30)
%     .MeasuredZmaxMm - end of the measured scan, for R² at desired depths (default Inf)
%
% Output:
%   res - struct with fields geometric, and if available library, fit,
%         kwave; each with .depths [1 x N], .z [N_z x 1], .Imeas, .Imod
%         [N_z x N], .r2 [1 x N], .focus [1 x N] and .settings (runtime
%         settings per depth, LIBRARY_MODEL_AT_DEPTH; [] for fits)
%   msg - summary lines (cell of char)
%
% See also: CALIBRATION_DESIRED_PROFILES, CALIBRATION_FIT, EVALUATE_LIBRARY_ENTRY,
%           FIT_TRANSDUCER_MODEL, CALIBRATION_FOCUS_BEYOND

    if nargin < 3; options = struct(); end
    predicted = logical(opt(options, 'Predicted', false));
    lib   = opt(options, 'Library', []);
    fit   = opt(options, 'Fit', []);
    kw    = opt(options, 'Kwave', []);
    isppa = opt(options, 'IsppaWcm2', 30);
    base  = opt(options, 'FitOptions', struct());
    skip  = opt(base, 'NearFieldSkipMm', 10);
    if ~isfield(prof, 'profile_type')
        prof.profile_type = repmat({'measured'}, 1, numel(prof.depths_ep_mm));
    end

    if predicted; base.NearFieldSkipMm = 0; end   % show the full modelled profile
    geo = base; geo.EvaluateOnly = true;
    rt  = base; rt.IsppaWcm2 = isppa;              % runtime path, drive at the reference ISPPA
    res = struct();
    has = any(isfinite(prof.I) & prof.I > 0, 1);  % depths with a (measured or derived) profile
    pg  = prof; pg.I(:, ~has) = 1;                % flat placeholder: the model scale is fitted to it
    res.geometric = pack(fit_transducer_model(tran, pg, geo));
    if ~isempty(lib)
        res.library = pack(evaluate_library_entry(tran, lib, prof, rt));
    end
    if ~isempty(fit)
        on_derived = any(~strcmp(fit.prof.profile_type, 'measured'));
        if predicted || on_derived
            % As a study resolves it (on the measured depths too, when the
            % fit was calibrated on desired depths)
            res.fit = pack(evaluate_library_entry(tran, opt(options, 'FitEntry', fit.entry), prof, rt));
        else
            res.fit = pack_many(fit.model);
            if ~isempty(kw); res.kwave = pack_kwave(kw.sims, kw.report, fit.model{1}.z_ep_mm, skip); end
        end
    end
    if predicted
        % The full profile is shown; R² and focus are taken beyond the near
        % field as in the fit (a near-field peak would otherwise count as
        % the focus), R² against the measured or derived profile
        zmax = opt(options, 'MeasuredZmaxMm', inf);
        for f = fieldnames(res)'
            r = res.(f{1});
            r.Imeas(:, ~has) = NaN;
            r.focus = calibration_focus_beyond(r.Imod, r.z, skip);
            r.r2 = r2_beyond(r.Imod, r.Imeas, r.z, skip, zmax);
            res.(f{1}) = r;
        end
    end
    msg = summary(res, prof, predicted, isppa);
end

function msg = summary(res, prof, predicted, isppa)
    has_lib = isfield(res, 'library');
    has_fit = isfield(res, 'fit');
    if ~predicted
        txt = sprintf('Geometric median R² %.3f', median(res.geometric.r2, 'omitnan'));
        if has_lib; txt = [txt sprintf(' | library %.3f', median(res.library.r2, 'omitnan'))]; end
        msg = {txt};
        return;
    end
    ferr = @(r) median(abs(r.focus - r.depths));
    txt = sprintf('geometric %.1f mm', ferr(res.geometric));
    if has_lib; txt = [txt sprintf(' | library %.1f mm', ferr(res.library))]; end
    if has_fit; txt = [txt sprintf(' | fit (unsaved) %.1f mm', ferr(res.fit))]; end
    pt = prof.profile_type;
    msg = {sprintf(['Predicted %d depths (profiles: %d measured, %d interpolated, %d extrapolated, ' ...
        '%d none). Median |focus − target|: %s'], numel(prof.depths_ep_mm), ...
        nnz(strcmp(pt, 'measured')), nnz(strcmp(pt, 'interpolated')), ...
        nnz(strcmp(pt, 'extrapolated')), nnz(strcmp(pt, 'unavailable')), txt)};
    src = ''; label = '';
    if has_fit && ~isempty(res.fit.settings)
        src = 'fit'; label = 'Fit (unsaved)';
    elseif has_lib && ~isempty(res.library.settings)
        src = 'library'; label = 'Library';
    end
    if ~isempty(src)
        pm = {res.(src).settings.phase_model};
        r2 = res.(src).r2;
        msg{end + 1} = sprintf(['%s at the desired depths: %d refined (within the refined span), %d blended ' ...
            'back to the static model, %d static/geometric; median R² %.3f on %d profiles ' ...
            '(phases and drive at %g W/cm² in the table).'], label, ...
            nnz(strcmp(pm, 'refined')), nnz(strcmp(pm, 'blend')), nnz(~ismember(pm, {'refined', 'blend'})), ...
            median(r2, 'omitnan'), nnz(isfinite(r2)), isppa);
    end
end

function r = pack(m)
    r = struct('depths', m.fit_depths_ep_mm(:)', 'z', m.z_ep_mm, 'Imeas', m.I_meas, ...
               'Imod', m.I_model, 'r2', m.metrics.r2(:)');
    r.focus = arrayfun(@(j) compute_flhm_center(m.I_model(:, j), m.z_ep_mm), 1:size(m.I_model, 2));
    r.settings = [];   % runtime settings per depth (EVALUATE_LIBRARY_ENTRY results only)
    if isfield(m, 'settings'); r.settings = m.settings; end
end

function r = pack_kwave(sims, rep, z, skip_mm)
% k-Wave on-axis profiles on the fit's z grid, peak-normalised beyond the near field
    r = struct('depths', [sims.depth_ep_mm], 'z', z, 'Imeas', [], 'Imod', zeros(numel(z), numel(sims)), ...
               'r2', rep.r2_kwave(:)');
    for k = 1:numel(sims)
        I = interp1(sims(k).z_ep_mm, sims(k).I_wcm2, z, 'linear', 0);
        r.Imod(:, k) = I / max(max(I(z >= skip_mm)), eps);
    end
    r.focus = arrayfun(@(j) compute_flhm_center(r.Imod(:, j), z), 1:numel(sims));
    r.settings = [];
end

function r = pack_many(models)
    r = pack(models{1});
    for k = 2:numel(models)
        q = pack(models{k});
        r.depths = [r.depths, q.depths]; r.Imeas = [r.Imeas, q.Imeas];
        r.Imod = [r.Imod, q.Imod]; r.r2 = [r.r2, q.r2]; r.focus = [r.focus, q.focus];
    end
end

function r2 = r2_beyond(Imod, Imeas, z, skip_mm, zmax)
% Profile-shape R² (peak-normalised intensity) beyond the near field and
% within the measured scan, as FIT_TRANSDUCER_MODEL; NaN without a profile
    k  = z >= skip_mm & z <= zmax;
    r2 = nan(1, size(Imod, 2));
    for j = 1:size(Imod, 2)
        yt = Imeas(k, j); ym = Imod(k, j);
        ok = isfinite(yt);
        if ~any(ok) || max(yt(ok)) <= 0; continue; end
        yt = yt(ok) / max(yt(ok)); ym = ym(ok) / max(max(ym(ok)), eps);
        r2(j) = 1 - sum((yt - ym).^2) / max(sum((yt - mean(yt)).^2), eps);
    end
end

function v = opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
