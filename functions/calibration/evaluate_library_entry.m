function ev = evaluate_library_entry(tran, lib, prof, options)
% EVALUATE_LIBRARY_ENTRY  Model the axial profiles a study gets from a library entry, at any depths
%
% For each requested depth the entry is resolved exactly as at run time
% (LIBRARY_MODEL_AT_DEPTH: static offset, refinements interpolated between
% refined depths and blended back to the static model outside them, drive
% and gains interpolated/clamped in depth), and the on-axis profile is
% modelled with FIT_TRANSDUCER_MODEL (EvaluateOnly). Depths need not be
% measured: this is how an entry is interpolated or extrapolated to a
% desired focal depth. Where prof holds a measurement the usual metrics
% (R², FLHM-centre error, ...) are reported; elsewhere they are NaN.
%
% Unlike scoring LIBRARY_STATIC_MODEL, this includes the refinements and
% does not apply the diagnostic fit_z_shift_mm (it is not applied at run
% time either).
%
% Use as:
%   ev = evaluate_library_entry(tran, lib, prof)
%   ev = evaluate_library_entry(tran, lib, prof, options)
%
% Input:
%   tran    - transducer struct from load_equipment_config (annular)
%   lib     - library struct (yaml.loadFile) or struct() for geometric steering
%   prof    - struct: .depths_ep_mm [1 x N_d], .z_ep_mm [N_z x 1] and
%             optionally .I [N_z x N_d] measured intensity; a column that is
%             all NaN/zero (or no .I) means "not measured at this depth"
%   options - struct, all optional:
%     .IsppaWcm2 - ISPPA for the reported drive [W/cm²] (default NaN: the
%                  entry's stored reference)
%     any FIT_TRANSDUCER_MODEL scoring option (ForwardModel, NearFieldSkipMm,
%     Domain, CWater, Rho, ...); FreeParams, Prior and EvaluateOnly are set here
%
% Output (ev), shaped like a FIT_TRANSDUCER_MODEL result:
%   .fit_depths_ep_mm - [1 x N_d] depths
%   .z_ep_mm          - [N_z' x 1] axial positions (near field excluded)
%   .I_model          - [N_z' x N_d] modelled intensity, peak-normalised
%   .I_meas           - [N_z' x N_d] measured intensity, peak-normalised (NaN = not measured)
%   .measured         - [1 x N_d] logical
%   .metrics          - table per depth as FIT_TRANSDUCER_MODEL, plus
%                       model_focus_ep_mm (FLHM centre of the model); metrics
%                       that need a measurement are NaN where measured is false
%   .settings         - [1 x N_d] struct array from LIBRARY_MODEL_AT_DEPTH
%                       (phases, drive, gains, phase_model, refinement_weight)
%   .forward_model
%
% See also: LIBRARY_MODEL_AT_DEPTH, LOAD_TRANSDUCER_FROM_LIBRARY, FIT_TRANSDUCER_MODEL,
%           CALIBRATION_GUI

arguments
    tran (1,1) struct
    lib
    prof (1,1) struct
    options (1,1) struct = struct()
end

    if nargin < 4; options = struct(); end
    tran  = library_transducer(tran, lib);    % sub-ring geometry if the entry has a ring profile
    isppa = NaN;
    if isfield(options, 'IsppaWcm2') && ~isempty(options.IsppaWcm2); isppa = options.IsppaWcm2; end
    fo = rmfield_if(options, {'IsppaWcm2', 'FreeParams', 'Prior', 'EvaluateOnly', 'ChannelGroups'});
    fo.EvaluateOnly = true;
    % One group per ring: the resolved offsets and gains are applied as they
    % are, not averaged over the rings of a channel
    fo.ChannelGroups = 1:tran.transducer.annular.elem_n;

    depths = prof.depths_ep_mm(:)';
    n_d    = numel(depths);
    z      = prof.z_ep_mm(:);
    if isfield(prof, 'I') && ~isempty(prof.I); I = prof.I; else; I = nan(numel(z), n_d); end
    measured = any(isfinite(I) & I > 0, 1);
    % The model is scaled to the measurement in closed form, so an empty
    % column would give a zero model: score unmeasured depths against a flat
    % placeholder (their metrics are discarded below)
    I(:, ~measured) = 1;

    settings = arrayfun(@(d) library_model_at_depth(lib, tran, d, isppa), depths);

    % Depths that resolve to the same model (e.g. every depth on the static
    % model) are scored in one call: the transfer matrix is built once
    key = arrayfun(@(s) {sprintf('%.6f,', [s.offset_deg, s.elem_amp_weight, s.steering_sound_speed])}, settings);
    [~, ~, grp] = unique(key, 'stable');
    parts = cell(1, max(grp));
    for g = 1:max(grp)
        j = find(grp == g);
        s = settings(j(1));
        o = fo;
        o.Prior = struct('elem_phase_correction_deg', s.offset_deg, ...
            'elem_amp_weight', s.elem_amp_weight, 'steering_sound_speed', s.steering_sound_speed);
        pg = struct('depths_ep_mm', depths(j), 'z_ep_mm', z, 'I', I(:, j));
        parts{g} = struct('cols', j, 'm', fit_transducer_model(tran, pg, o));
    end

    m1 = parts{1}.m;
    ev.fit_depths_ep_mm = depths;
    ev.z_ep_mm = m1.z_ep_mm;
    ev.I_model = nan(numel(m1.z_ep_mm), n_d);
    ev.I_meas  = nan(numel(m1.z_ep_mm), n_d);
    metrics = repmat(m1.metrics(1, :), n_d, 1);
    for g = 1:numel(parts)
        j = parts{g}.cols; m = parts{g}.m;
        ev.I_model(:, j) = m.I_model;
        ev.I_meas(:, j)  = m.I_meas;
        metrics(j, :)    = m.metrics;
    end
    metrics.model_focus_ep_mm = arrayfun(@(k) compute_flhm_center(ev.I_model(:, k), ev.z_ep_mm), (1:n_d)');
    no = ~measured(:);
    for v = {'r2', 'r2_pressure', 'focal_err_mm', 'peak_err_mm', 'nearfield_excess'}
        metrics.(v{1})(no) = NaN;
    end
    ev.I_meas(:, no) = NaN;
    ev.measured      = measured;
    ev.metrics       = metrics;
    ev.settings      = settings;
    ev.forward_model = m1.forward_model;
end

function s = rmfield_if(s, f)
    s = rmfield(s, intersect(fieldnames(s), f));
end
