function calibration_setup(parameters, equip_param)
% CALIBRATION_SETUP  Run CALIBRATION_TRANSDUCER for every equipment combination in the calibration config
%
% For each entry of calibration.combinations: locate its hydrophone profile
% CSV (equipment YAML char_data_path, relative to calibration.path_input_axial),
% keep the requested focal depths within the transducer's focal range, and
% calibrate all depths and intensities in one CALIBRATION_PIPELINE_START call
% (locally or as one HPC job per combination).
%
% Use as:
%   calibration_setup(parameters, equip_param)
%
% Input:
%   parameters  - PRESTUS config with a calibration section (config_calibration.yaml)
%   equip_param - struct from load_equipment_config
%
% See also: CALIBRATION_TRANSDUCER, CALIBRATION_PIPELINE_START, CALIBRATE_LIBRARY_ENTRY

if ~isfolder(parameters.calibration.path_output)
    mkdir(parameters.calibration.path_output);
end

N_i = length(parameters.calibration.combinations);

% Coerce list-of-list fields to cell-of-cell if yaml ConvertToArray
% flattened them into numeric matrices.
%
% Ambiguous case: when every combo's list has exactly one value (e.g.
% N_i combos each requesting a single focal depth), ConvertToArray
% flattens "[[85],[85]]" into a plain 1xN_i (or N_ix1) vector that is
% structurally indistinguishable from "one combo's list of N_i values".
% Since N_i (from `combinations`) is known, a flattened vector whose length
% equals N_i is treated as one scalar per combo.
for coerce_field = {'focal_depths_wrt_exit_plane', 'desired_intensities'}
    f = coerce_field{1};
    val = parameters.calibration.(f);
    if ~iscell(val)
        if N_i > 1 && isvector(val) && numel(val) == N_i
            parameters.calibration.(f) = arrayfun(@(v) {v}, val(:)', 'UniformOutput', false);
        else
            parameters.calibration.(f) = ...
                cellfun(@num2cell, num2cell(val, 2), 'UniformOutput', false);
        end
    elseif ~iscell(val{1})
        parameters.calibration.(f) = ...
            cellfun(@num2cell, val, 'UniformOutput', false);
    end
end

% If a list-of-list field has only one entry, replicate it for all combinations
for coerce_field = {'focal_depths_wrt_exit_plane', 'desired_intensities'}
    f = coerce_field{1};
    if iscell(parameters.calibration.(f)) && length(parameters.calibration.(f)) == 1 && N_i > 1
        parameters.calibration.(f) = repmat(parameters.calibration.(f), 1, N_i);
    end
end

for i = 1:N_i
    combo_name = char(parameters.calibration.combinations{i});
    combo      = equip_param.combos.(combo_name);
    tran       = equip_param.trans.(char(combo.tran_serial));
    fprintf('Equipment: transducer %s and driving system %s\n', tran.name, combo.ds.name);

    % Hydrophone profiles (columns = focal depths from the exit plane)
    csv_path = combo.char_data_path;
    if ~isempty(parameters.calibration.path_input_axial)
        csv_path = fullfile(parameters.calibration.path_input_axial, csv_path);
    end

    % Former per-element correction stored in the equipment YAML: used as
    % the prior the fits are regularised towards
    cal_i = parameters;
    cal_i.calibration.elem_phase_correction_deg = [];
    if isfield(tran, 'elem_phase_correction') && ~isempty(tran.elem_phase_correction)
        cal_i.calibration.elem_phase_correction_deg = tran.elem_phase_correction.deg(:)';
    end

    % Requested depths ([] = every measured depth) within the focal range
    depths = cell2mat(parameters.calibration.focal_depths_wrt_exit_plane{i});
    if ~isempty(depths)
        out = depths < tran.min_foc | depths > tran.max_foc;
        if any(out)
            warning('calibration_setup:outOfRange', ...
                'Focus %s mm outside the range [%.2f, %.2f] mm. Skipping.', ...
                mat2str(depths(out)), tran.min_foc, tran.max_foc);
            depths = depths(~out);
            if isempty(depths); continue; end
        end
    end
    intensities = cell2mat(parameters.calibration.desired_intensities{i});

    fprintf('Calibrating %s | F %s mm | I %s W/cm²\n', combo_name, ...
        mat2str(round(depths, 2)), mat2str(intensities));
    calibration_pipeline_start(csv_path, combo_name, intensities, depths, cal_i, i);
end
end
