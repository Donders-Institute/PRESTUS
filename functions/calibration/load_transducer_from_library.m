function parameters = load_transducer_from_library(combo_name, focal_distance_ep, desired_intensity, equip_param, library_path)
% LOAD_TRANSDUCER_FROM_LIBRARY  Load calibrated transducer parameters from library
%
% Reads the calibration library YAML for the given equipment combination and
% the geometry from equip_param (config_equipment.yaml). Reconstructs per-
% element phases from the fitted parametric model and looks up (or scales)
% the amplitude for the desired intensity.
%
% This is the PRESTUS equivalent of how BabelBrain loads a device YAML at
% runtime: geometry + fitted model -> ready-to-use transducer parameters.
% The difference is that BabelBrain trusts the manufacturer phase table
% directly, while PRESTUS uses phases fitted by calibration_transducer.
%
% Use as:
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param)
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param, library_path)
%
% Input:
%   combo_name         - equipment combination key, e.g.
%                        'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'
%   focal_distance_ep  - desired focal distance from exit plane [mm]
%   desired_intensity  - target ISPPA [W/cm²]
%   equip_param        - struct from yaml.loadFile('config_equipment.yaml')
%   library_path       - (optional) path to transducer library folder;
%                        defaults to config/transducer/ under PRESTUS root
%
% Output:
%   parameters - struct with:
%     .transducer          — geometry + calibrated phases and amplitude
%     .library_meta        — provenance from the library YAML meta section
%     .library_focal_used  — actual calibrated focal depth used [mm]
%
% See also: UPDATE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER

    if nargin < 5 || isempty(library_path)
        library_path = fullfile(get_PRESTUSpath(), 'config', 'transducer');
    end

    %% Load calibration library entry
    yaml_path = fullfile(library_path, [combo_name, '.yaml']);
    if ~isfile(yaml_path)
        error('load_transducer_from_library: library file not found:\n  %s', yaml_path);
    end
    lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);

    %% Load geometry from equipment config
    tran_serial = lib.meta.tran_serial;
    if ~isfield(equip_param.trans, tran_serial)
        error('load_transducer_from_library: transducer ''%s'' not found in equip_param.', tran_serial);
    end
    tran = equip_param.trans.(tran_serial);
    parameters.transducer = tran.prestus.transducer;

    %% Find closest calibrated focal depth
    cal = lib.calibration.focal_depths;
    if isempty(fieldnames(cal))
        error(['load_transducer_from_library: no calibrated focal depths in library ' ...
            'for combo ''%s''. Run calibration_standalone first.'], combo_name);
    end

    depth_keys   = fieldnames(cal);
    depth_values = cellfun(@(k) str2double(strrep(strrep(k, 'f', ''), 'p', '.')), depth_keys);

    [dist, idx]    = min(abs(depth_values - focal_distance_ep));
    used_key       = depth_keys{idx};
    focal_used     = depth_values(idx);

    if dist > 5
        warning(['load_transducer_from_library: closest calibrated depth is %.1f mm ' ...
            '(requested %.1f mm). Consider running calibration at this depth.'], ...
            focal_used, focal_distance_ep);
    end

    entry  = cal.(used_key);
    n_elem = tran.prestus.transducer.annular.elem_n;

    %% Reconstruct per-element phases from parametric model
    phase_start_rad = entry.phase_start_deg / 180 * pi;

    switch entry.precession_mode
        case 'linear'
            phase_step_rad = entry.phase_step_deg / 180 * pi;
            elem_phase_rad = mod(phase_start_rad + (0:n_elem-1) * phase_step_rad, 2*pi);
        case 'monotonic'
            increments_rad = entry.phase_increments_deg / 180 * pi;
            elem_phase_rad = mod(phase_start_rad + [0, cumsum(increments_rad)], 2*pi);
        otherwise
            error('load_transducer_from_library: unknown precession_mode ''%s''.', ...
                entry.precession_mode);
    end

    parameters.transducer.annular.elem_phase_rad = elem_phase_rad;
    parameters.transducer.annular.elem_phase_deg = elem_phase_rad / pi * 180;

    %% Look up amplitude for desired intensity
    amp_map     = entry.amplitude_scaling;
    int_keys    = fieldnames(amp_map);
    int_values  = cellfun(@(k) str2double(strrep(strrep(k, 'i', ''), 'p', '.')), int_keys);

    [int_dist, amp_idx] = min(abs(int_values - desired_intensity));
    int_used            = int_values(amp_idx);

    if int_dist / desired_intensity > 0.05
        % Scale analytically: amp ∝ sqrt(I) since I ∝ v² ∝ amp²
        base_amp = amp_map.(int_keys{amp_idx});
        elem_amp = round(base_amp * sqrt(desired_intensity / int_used));
        warning(['load_transducer_from_library: scaling amplitude from %.1f to %.1f W/cm² ' ...
            '(%.0f -> %d). Run calibration at this intensity for an exact value.'], ...
            int_used, desired_intensity, base_amp, elem_amp);
    else
        elem_amp = amp_map.(int_keys{amp_idx});
    end

    parameters.transducer.annular.elem_amp = repmat(elem_amp, 1, n_elem);

    %% Focal distance bookkeeping
    parameters.transducer.focal_distance_ep   = focal_distance_ep;
    parameters.transducer.focal_distance_bowl = focal_distance_ep + ...
        (tran.prestus.transducer.annular.curv_radius_mm - ...
         tran.prestus.transducer.annular.dist_geom_ep_mm);
    parameters.transducer.set_intensity_w_per_cm2 = desired_intensity;

    %% Provenance
    parameters.library_meta       = lib.meta;
    parameters.library_focal_used = focal_used;

    fprintf('Library loaded: %s | focal %.1f mm (calibrated %.1f mm) | amp %d @ %.1f W/cm²\n', ...
        combo_name, focal_distance_ep, focal_used, elem_amp(1), desired_intensity);
end
