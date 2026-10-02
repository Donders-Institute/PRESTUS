function s = library_condition(lib, tran, label, desired_intensity, source_label)
% LIBRARY_CONDITION  Phases, gains and drive of a dephased library condition [EXPERIMENTAL]
%
% Reads lib.conditions.<label> (written by LIBRARY_ADD_CONDITIONS) and
% returns the same fields as LIBRARY_MODEL_AT_DEPTH, so the loader treats a
% condition like a depth. The stored drive is scaled as sqrt(I) to
% desired_intensity. Warns that the condition is experimental unless the
% entry has been marked validated: true.
%
% Use as:
%   s = library_condition(lib, tran, label, desired_intensity)
%   s = library_condition(lib, tran, label, desired_intensity, source_label)
%
% Input:
%   lib               - library struct
%   tran              - annular transducer struct from load_equipment_config
%   label             - condition name (transducer.condition)
%   desired_intensity - target free-water ISPPA [W/cm²]; NaN = stored reference
%   source_label      - (optional) entry name for messages
%
% Output (s): as LIBRARY_MODEL_AT_DEPTH, with phase_model = 'condition',
%   amplitude_source = 'condition', plus
%   .condition       - the stored condition struct
%   .setpoint_ep_mm  - nominal focal setpoint of the condition (NaN if none)
%
% See also: LIBRARY_ADD_CONDITIONS, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 5 || isempty(source_label); source_label = '(in-memory entry)'; end
    label = char(label);
    if ~isfield(lib, 'conditions') || ~isstruct(lib.conditions) || ~isfield(lib.conditions, label)
        avail = {'(none)'};
        if isfield(lib, 'conditions') && isstruct(lib.conditions) && ~isempty(fieldnames(lib.conditions))
            avail = fieldnames(lib.conditions);
        end
        error('library_condition:unknownCondition', ...
            'Condition ''%s'' not found in %s. Available: %s', label, source_label, strjoin(avail, ', '));
    end
    cd = lib.conditions.(label);
    n_elem = tran.transducer.annular.elem_n;

    phase = double(cd.elem_phase_deg(:)');
    gains = double(cd.elem_amp_weight(:)');
    if numel(phase) ~= n_elem || numel(gains) ~= n_elem
        error('library_condition:phaseLength', ...
            'Condition ''%s'' in %s has %d phases for a %d-ring transducer (built for another geometry?).', ...
            label, source_label, numel(phase), n_elem);
    end

    if ~(isfield(cd, 'validated') && isequal(cd.validated, true))
        warning('library_condition:experimental', ...
            ['%s: dephased condition ''%s'' is EXPERIMENTAL and not validated against hydrophone ' ...
             'measurements (predicted R² %.3g). Phases assume the focused-profile hardware correction ' ...
             'transfers; the drive scales the model''s free-water peak at %.1f mm from the exit plane.'], ...
            source_label, label, double(cd.predicted_r2), double(cd.peak_z_ep_mm));
    end

    drive = double(cd.elem_p_amp);
    if ~isnan(desired_intensity)
        drive = round(drive * sqrt(desired_intensity / double(cd.isppa_wcm2)));
    end

    s = struct();
    s.offset_deg            = double(cd.elem_phase_correction_deg(:)');
    s.elem_phase_deg        = mod(phase, 360);
    s.elem_phase_rad        = s.elem_phase_deg * pi/180;
    s.steering_sound_speed  = NaN;
    s.drive_pa              = drive;
    s.elem_amp_weight       = gains;
    s.elem_p_amp            = round(drive * gains);
    s.phase_model           = 'condition';
    s.refinement_weight     = NaN;
    s.amplitude_source      = 'condition';
    s.amp_kwave_calibrated  = false;
    s.amp_kwave_nearest_depth_mm = NaN;
    s.amp_kwave_r2          = NaN;
    s.amp_kwave_focal_err_mm = NaN;
    s.condition             = cd;
    s.setpoint_ep_mm        = double(cd.setpoint_ep_mm);
end
