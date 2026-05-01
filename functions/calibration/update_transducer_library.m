function update_transducer_library(parameters, opt_params_raw, precession_mode, library_path)
% UPDATE_TRANSDUCER_LIBRARY  Write calibration results into a transducer library YAML
%
% Adds or updates the focal-depth entry for the current calibration run in the
% combo's library YAML. Called automatically by save_parametric_model when
% parameters.calibration.library_combo is set.
%
% Use as:
%   update_transducer_library(parameters, opt_params_raw, precession_mode)
%   update_transducer_library(parameters, opt_params_raw, precession_mode, library_path)
%
% Input:
%   parameters      - PRESTUS config after calibration; must contain:
%                       calibration.library_combo    — combo name string
%                       calibration.desired_focal_distance_ep [mm]
%                       calibration.desired_intensity [W/cm²]
%                       transducer.annular.elem_amp
%   opt_params_raw  - raw optimisation vector from perform_global_search
%   precession_mode - 'linear' or 'monotonic'
%   library_path    - (optional) path to transducer library folder;
%                     defaults to config/transducer/ under PRESTUS root
%
% See also: SAVE_PARAMETRIC_MODEL, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 4 || isempty(library_path)
        library_path = fullfile(get_PRESTUSpath(), 'config', 'transducer');
    end

    combo_name = parameters.calibration.library_combo;
    yaml_path  = fullfile(library_path, [combo_name, '.yaml']);

    if ~isfile(yaml_path)
        error('update_transducer_library: library file not found:\n  %s\n', yaml_path);
    end

    lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);

    focal_ep    = parameters.calibration.desired_focal_distance_ep;
    intensity   = parameters.calibration.desired_intensity;
    elem_amp    = double(parameters.transducer.annular.elem_amp(1));
    n_elem      = parameters.transducer.annular.elem_n;

    % YAML struct field keys cannot start with digits; encode focal depth as
    % 'fXXX' with decimal point replaced by 'p' (e.g. 85.0 -> f85p0).
    focal_key = sprintf('f%s', strrep(sprintf('%.1f', focal_ep), '.', 'p'));
    % Intensity key: 'iXXX' with decimal point replaced by 'p'.
    int_key   = sprintf('i%s', strrep(sprintf('%.1f', intensity), '.', 'p'));

    %% Build focal-depth entry
    entry.precession_mode = precession_mode;
    entry.phase_start_deg = opt_params_raw(1) / pi * 180;

    switch precession_mode
        case 'linear'
            entry.phase_step_deg = opt_params_raw(2) / pi * 180;
        case 'monotonic'
            entry.phase_increments_deg = opt_params_raw(2:n_elem) / pi * 180;
    end

    entry.velocity_m_s = opt_params_raw(end);
    entry.validated    = false;

    %% Preserve or initialise amplitude scaling map
    if isfield(lib.calibration.focal_depths, focal_key)
        existing = lib.calibration.focal_depths.(focal_key);
        if isfield(existing, 'amplitude_scaling')
            entry.amplitude_scaling = existing.amplitude_scaling;
        else
            entry.amplitude_scaling = struct();
        end
        % Preserve validated flag if previously set to true
        if isfield(existing, 'validated') && existing.validated
            entry.validated = true;
        end
    else
        entry.amplitude_scaling = struct();
    end

    entry.amplitude_scaling.(int_key) = elem_amp;

    %% Update library
    lib.calibration.focal_depths.(focal_key) = entry;
    lib.meta.calibration_date = datestr(now, 'yyyy-mm-dd');

    yaml.dumpFile(yaml_path, lib);
    fprintf('Library updated: %s | focal %s mm | %s W/cm² -> amp %d\n', ...
        combo_name, num2str(focal_ep), num2str(intensity), elem_amp);
end
