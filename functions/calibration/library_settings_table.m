function T = library_settings_table(tran, lib, depths_ep_mm, options)
% LIBRARY_SETTINGS_TABLE  Per-depth transducer settings a library entry resolves to, as a table / CSV
%
% For each focal depth: the phases, drive and gains a study gets from the
% entry (LIBRARY_MODEL_AT_DEPTH, the runtime path of
% LOAD_TRANSDUCER_FROM_LIBRARY), which model they come from, and the
% modelled focus (EVALUATE_LIBRARY_ENTRY). Depths need not be measured or
% refined; outside the refined span the entry blends back to its static model.
%
% Use as:
%   T = library_settings_table(tran, lib, depths_ep_mm)
%   T = library_settings_table(tran, lib, depths_ep_mm, options)
%
% Input:
%   tran         - transducer struct from load_equipment_config (annular)
%   lib          - library struct (yaml.loadFile) or struct() for geometric steering
%   depths_ep_mm - focal depths [mm from EP]
%   options      - struct, all optional:
%     .IsppaWcm2       - ISPPA for the drive [W/cm²] (default NaN: stored reference)
%     .CsvPath         - also write T to this CSV
%     .ProfileType     - {1 x N} profile each depth was calibrated on, e.g. from
%                        INTERPOLATE_AXIAL_PROFILES ('' = none), added as a column
%     ForwardModel, NearFieldSkipMm, ... - passed to EVALUATE_LIBRARY_ENTRY (focus)
%
% Output (T, one row per depth):
%   depth_ep_mm, phase_model, refinement_weight, amplitude_source,
%   isppa_wcm2, drive_pa (mean source pressure), model_focus_ep_mm,
%   [profile_type], elem_phase_deg_1..N, elem_p_amp_pa_1..N, elem_amp_weight_1..N
%
% See also: LIBRARY_MODEL_AT_DEPTH, EVALUATE_LIBRARY_ENTRY, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 4; options = struct(); end
    csv   = '';
    ptype = {};
    if isfield(options, 'CsvPath');     csv   = options.CsvPath;     options = rmfield(options, 'CsvPath'); end
    if isfield(options, 'ProfileType'); ptype = options.ProfileType; options = rmfield(options, 'ProfileType'); end
    isppa = NaN;
    if isfield(options, 'IsppaWcm2') && ~isempty(options.IsppaWcm2); isppa = options.IsppaWcm2; end

    d  = depths_ep_mm(:)';
    zmax = ceil(1.5 * max([d, tran.transducer.annular.dist_geom_ep_mm]));
    ev = evaluate_library_entry(tran, lib, struct('depths_ep_mm', d, 'z_ep_mm', (0:0.1:zmax)'), options);
    s  = ev.settings;
    n_elem = numel(s(1).elem_phase_deg);

    drive = arrayfun(@(x) empty_nan(x.drive_pa), s)';
    pamp  = nan(numel(s), n_elem);
    for k = 1:numel(s); if ~isempty(s(k).elem_p_amp); pamp(k, :) = s(k).elem_p_amp; end; end
    T = table(d(:), string({s.phase_model})', [s.refinement_weight]', string({s.amplitude_source})', ...
        repmat(isppa, numel(d), 1), drive, ev.metrics.model_focus_ep_mm, ...
        'VariableNames', {'depth_ep_mm', 'phase_model', 'refinement_weight', 'amplitude_source', ...
        'isppa_wcm2', 'drive_pa', 'model_focus_ep_mm'});
    if ~isempty(ptype); T.profile_type = string(ptype(:)); end
    T = [T, array2table(vertcat(s.elem_phase_deg), 'VariableNames', compose('elem_phase_deg_%d', 1:n_elem)), ...
            array2table(pamp, 'VariableNames', compose('elem_p_amp_pa_%d', 1:n_elem)), ...
            array2table(vertcat(s.elem_amp_weight), 'VariableNames', compose('elem_amp_weight_%d', 1:n_elem))];
    if ~isempty(csv)
        folder = fileparts(csv);
        if ~isempty(folder) && ~isfolder(folder); mkdir(folder); end
        writetable(T, csv);
    end
end

function v = empty_nan(x)
    if isempty(x); v = NaN; else; v = x; end
end
