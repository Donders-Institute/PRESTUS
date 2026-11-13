function [lib, report] = library_add_conditions(lib, tran, conds, options)
% LIBRARY_ADD_CONDITIONS  Add dephased drive conditions to a library entry [EXPERIMENTAL]
%
% EXPERIMENTAL — NOT VALIDATED. Dephased (defocused) conditions drive the
% rings with explicit phases instead of geometric steering to a focal depth,
% so they cannot be stored as depth refinements. Each condition is stored
% under lib.conditions.<label> and selected in a study with
% transducer.condition: <label>.
%
% Phases of a condition:
%   absolute: phase = elem_phase_deg                         + offset(setpoint)
%   shift:    phase = geometric steering phases at setpoint
%                     + phase_shift_deg                      + offset(setpoint)
% where offset(setpoint) is the entry's hardware phase correction at the
% condition's nominal setpoint (LIBRARY_MODEL_AT_DEPTH), applied together
% with the per-ring gains of that setpoint. This assumes the correction
% fitted to focused profiles transfers to dephased drives — the assumption
% the measured dephased profiles have to confirm. Set InheritCorrection to
% false (or give no setpoint) to use the nominal phases and unit gains.
%
% The drive is chosen so that the model's free-water on-axis peak (beyond
% NearFieldSkipMm) reaches IsppaWcm2; for a dephased beam this peak can be
% a near-field or split focus. Measured profiles are used only to score the
% prediction (R², focal error), never to fit it.
%
% Use as:
%   lib = library_add_conditions(lib, tran, conds)
%   [lib, report] = library_add_conditions(lib, tran, conds, options)
%
% Input:
%   lib     - library struct (yaml.loadFile(..., 'ConvertToArray', true))
%   tran    - annular transducer struct from load_equipment_config
%   conds   - struct array, one per condition:
%               .label            - condition name (valid YAML/MATLAB field name)
%               .elem_phase_deg   - [1 x N_elem] absolute driver phases [deg], or
%               .phase_shift_deg  - [1 x N_elem] shift added to the setpoint's
%                                   steering phases [deg] (needs .setpoint_ep_mm)
%               .setpoint_ep_mm   - nominal focal setpoint the condition derives
%                                   from [mm from EP] (optional for absolute phases)
%               .profile_column   - (optional) column of options.Profile measured
%                                   for this condition
%   options - struct, all optional:
%     .IsppaWcm2         - reference ISPPA for the stored drive (default 30)
%     .InheritCorrection - apply the setpoint's phase offset and gains (default true)
%     .Profile           - measured profiles (READ_AXIAL_PROFILES) for scoring
%     .SourceCsv         - provenance string stored with each condition
%     .NearFieldSkipMm   - ignore z < this for the peak and the scores (default 10)
%     .CWater, .Rho      - medium (default 1500 m/s, 994 kg/m³)
%     .ForwardModel      - 'rayleigh' (default) | 'oneil'
%
% Output:
%   lib    - entry with lib.conditions.<label> added or replaced
%   report - table: label, r2, focal_err_mm, peak_err_mm, peak_z_ep_mm, elem_p_amp
%
% See also: LIBRARY_CONDITION, LIBRARY_INGEST_CONDITIONS, LIBRARY_MODEL_AT_DEPTH

    if nargin < 4; options = struct(); end
    o = struct('IsppaWcm2', 30, 'InheritCorrection', true, 'Profile', [], ...
               'SourceCsv', '', 'NearFieldSkipMm', 10, 'CWater', 1500, 'Rho', 994, ...
               'ForwardModel', 'rayleigh');
    fn = fieldnames(options);
    for k = 1:numel(fn); o.(fn{k}) = options.(fn{k}); end

    ann    = tran.transducer.annular;
    n_elem = ann.elem_n;
    ep_to_bowl = annular_bowl_depth_mm(ann);

    % Axial grid: the measured positions when scoring, else 0..1.5 ROC
    if ~isempty(o.Profile)
        z_ep = o.Profile.z_ep_mm(:);
    else
        z_ep = (0:0.25:1.5 * ann.curv_radius_mm)';
    end
    H  = axial_transfer_matrix(tran, z_ep + ep_to_bowl, o.CWater, o.Rho, o.ForwardModel);
    ff = z_ep >= o.NearFieldSkipMm;

    if ~isfield(lib, 'conditions') || ~isstruct(lib.conditions); lib.conditions = struct(); end

    n_c = numel(conds);
    r2 = nan(n_c, 1); fe = r2; pe = r2; pz = r2; pa = r2; labels = cell(n_c, 1);
    for c = 1:n_c
        cd = conds(c);
        label = char(cd.label);
        if ~isvarname(label)
            error('library_add_conditions:badLabel', ...
                'Condition label ''%s'' must be a valid field name (letters, digits, underscores).', label);
        end
        labels{c} = label;
        setpoint = get_or(cd, 'setpoint_ep_mm', []);
        has_abs   = isfield(cd, 'elem_phase_deg')  && ~isempty(cd.elem_phase_deg);
        has_shift = isfield(cd, 'phase_shift_deg') && ~isempty(cd.phase_shift_deg);
        if has_abs == has_shift
            error('library_add_conditions:phaseSpec', ...
                'Condition ''%s'': give exactly one of elem_phase_deg or phase_shift_deg.', label);
        end
        if has_shift && isempty(setpoint)
            error('library_add_conditions:noSetpoint', ...
                'Condition ''%s'': phase_shift_deg needs setpoint_ep_mm.', label);
        end

        % Setpoint model: geometric phases, hardware offset and gains there
        offset = zeros(1, n_elem); gains = ones(1, n_elem); geo = zeros(1, n_elem);
        if ~isempty(setpoint)
            m   = library_model_at_depth(lib, tran, setpoint, NaN);
            geo = mod(m.elem_phase_deg - m.offset_deg, 360);
            if o.InheritCorrection
                offset = m.offset_deg;
                if ~isempty(m.elem_amp_weight); gains = m.elem_amp_weight(:)'; end
            end
        end
        if has_abs
            nominal = check_n(cd.elem_phase_deg, n_elem, label, 'elem_phase_deg');
            mode = 'absolute';
        else
            nominal = mod(geo + check_n(cd.phase_shift_deg, n_elem, label, 'phase_shift_deg'), 360);
            mode = 'shift';
        end
        phase = mod(nominal + offset, 360);

        % Drive: mean element velocity giving IsppaWcm2 at the far-field peak
        I_unit = abs(H * (gains(:) .* exp(1i * phase(:) * pi/180))).^2 / (2 * o.Rho * o.CWater) * 1e-4;
        [I_pk, k_pk] = max(I_unit .* ff);
        v    = sqrt(o.IsppaWcm2 / I_pk);
        drive = v * o.Rho * o.CWater;                               % mean source pressure [Pa]

        % Score against the measured profile of this condition, if any
        col = get_or(cd, 'profile_column', []);
        if ~isempty(o.Profile) && ~isempty(col)
            yt = o.Profile.I(:, col);  ym = I_unit;
            yt = yt(ff) / max(max(yt(ff)), eps);  ym = ym(ff) / max(max(ym(ff)), eps);  zf = z_ep(ff);
            r2(c) = 1 - sum((yt - ym).^2) / max(sum((yt - mean(yt)).^2), eps);
            fe(c) = compute_flhm_center(ym, zf) - compute_flhm_center(yt, zf);
            [~, a1] = max(ym); [~, a2] = max(yt);
            pe(c) = zf(a1) - zf(a2);
        end
        pz(c) = z_ep(k_pk);
        pa(c) = drive;

        entry = struct();
        entry.status              = 'experimental';
        entry.validated           = false;
        entry.phase_mode          = mode;
        entry.nominal_phase_deg   = round(nominal, 4);
        if has_shift; entry.phase_shift_deg = round(cd.phase_shift_deg(:)', 4); end
        entry.setpoint_ep_mm      = empty_to_nan(setpoint);
        entry.inherit_correction  = logical(o.InheritCorrection && ~isempty(setpoint));
        entry.elem_phase_correction_deg = round(offset, 4);
        entry.elem_phase_deg      = round(phase, 4);
        entry.elem_amp_weight     = round(gains, 6);
        entry.elem_p_amp          = round(drive);
        entry.isppa_wcm2          = o.IsppaWcm2;
        entry.peak_z_ep_mm        = pz(c);
        entry.forward_model       = o.ForwardModel;
        entry.near_field_skip_mm  = o.NearFieldSkipMm;
        entry.predicted_r2        = r2(c);
        entry.predicted_focal_err_mm = fe(c);
        entry.predicted_peak_err_mm  = pe(c);
        entry.source_csv          = char(o.SourceCsv);
        entry.created_at          = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
        lib.conditions.(label) = entry;
    end

    report = table(labels, r2, fe, pe, pz, pa, 'VariableNames', ...
        {'label', 'r2', 'focal_err_mm', 'peak_err_mm', 'peak_z_ep_mm', 'elem_p_amp'});
end

function v = check_n(v, n, label, name)
    v = double(v(:)');
    if numel(v) ~= n
        error('library_add_conditions:phaseLength', ...
            'Condition ''%s'': %s has %d values, the transducer has %d rings.', label, name, numel(v), n);
    end
end

function v = get_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function v = empty_to_nan(v)
    if isempty(v); v = NaN; end
end
