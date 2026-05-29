function save_elem_correction(equipment_yaml_path, opt_phases_deg, geo_phases_deg, ref_depth_ep)
% SAVE_ELEM_CORRECTION  Compute and persist per-element phase correction into the equipment YAML
%
% Derives the hardware correction as the difference between optimised and
% geometric phases at a reference depth (unwrapped before subtraction to
% avoid wrap-around artefacts), then writes it into the transducer equipment
% YAML under elem_phase_correction so it is loaded automatically by
% load_equipment_config / load_transducer_from_library.
%
% Use as:
%   save_elem_correction(equipment_yaml_path, opt_phases_deg, geo_phases_deg, ref_depth_ep)
%
% Input:
%   equipment_yaml_path - full path to the transducer equipment YAML
%                         (e.g. {equip_param.path}/{tran_serial}.yaml)
%   opt_phases_deg      - optimised element phases at the reference depth [°, 1xN]
%   geo_phases_deg      - geometric (compute_phases) element phases at same depth [°, 1xN]
%   ref_depth_ep        - reference focal depth w.r.t. exit plane [mm]
%
% The correction is stored as:
%   elem_phase_correction:
%     ref_depth_ep_mm: <value>
%     deg: [e1, e2, ...]
%
% See also: CALIBRATION_EXTRACT_CORRECTION, CALIBRATION_TRANSDUCER, COMPUTE_PHASES

arguments
    equipment_yaml_path (1,:) char
    opt_phases_deg      (1,:) {mustBeNumeric}
    geo_phases_deg      (1,:) {mustBeNumeric}
    ref_depth_ep        (1,1) {mustBeNumeric}
end

    if numel(opt_phases_deg) ~= numel(geo_phases_deg)
        error('save_elem_correction: opt_phases_deg and geo_phases_deg must have the same length.');
    end
    if ~isfile(equipment_yaml_path)
        error('save_elem_correction: equipment YAML not found: %s', equipment_yaml_path);
    end

    % Unwrap before subtracting to avoid wrap-around artefacts
    geo_unwrapped  = unwrap(geo_phases_deg * pi/180) * 180/pi;
    opt_unwrapped  = unwrap(opt_phases_deg * pi/180) * 180/pi;
    correction_deg = mod(opt_unwrapped - geo_unwrapped, 360);

    % Load existing YAML, add/overwrite the correction field, write back
    d = yaml.loadFile(equipment_yaml_path, 'ConvertToArray', true);
    d.elem_phase_correction.ref_depth_ep_mm = ref_depth_ep;
    d.elem_phase_correction.deg             = correction_deg;
    yaml.dumpFile(equipment_yaml_path, d);

    fprintf('Element correction (%d elements, ref depth %.2f mm) written to: %s\n', ...
        numel(correction_deg), ref_depth_ep, equipment_yaml_path);

end
