function save_optimized_values(parameters)
% SAVE_OPTIMIZED_VALUES  Save optimised transducer phases and amplitude to a CSV file
%
% Writes or updates a calibration CSV where rows correspond to intensity
% targets and columns to focal depths. Each cell stores a string of
% element phases and amplitude.
%
% Use as:
%   save_optimized_values(parameters)
%
% Input:
%   parameters - PRESTUS config; uses calibration.path_output_profiles,
%                calibration.filename_calibrated_CSV,
%                calibration.desired_focal_distance_ep [mm],
%                calibration.desired_intensity [W/cm²],
%                calibration.equipment_name, and transducer.annular phases/amplitude
%
% See also: CALIBRATION_TRANSDUCER, SET_REAL_PHASES

arguments
    parameters (1,1) struct
end
    
    disp('Saving optimized values to CSV file...');

    output_file_path = fullfile(parameters.calibration.path_output_profiles, ...
        parameters.calibration.filename_calibrated_CSV);
    
    % Extract and round phases and amplitudes
    opt_phases = round(parameters.transducer.annular.elem_phase_deg, 2);
    elem_p_amp = double(parameters.transducer.annular.elem_p_amp(1));

    fprintf('CSV file can be found here: %s \n', output_file_path);
    
    % Check if the CSV file exists
    if isfile(output_file_path)
        % Read existing data
        virtual_data = readcell(output_file_path, 'Delimiter', ',');
        prestus_foci_wrt_exit_plane = round(cell2mat(virtual_data(1, 2:end)), 2);
        prestus_int = cell2mat(virtual_data(2:end, 1));

        % Find or add the focal distance column
        col_index_foc = find(prestus_foci_wrt_exit_plane == parameters.calibration.desired_focal_distance_ep, 1);
        if isempty(col_index_foc)
            col_index_foc = size(virtual_data, 2) + 1;
            virtual_data{1, col_index_foc} = parameters.calibration.desired_focal_distance_ep;
        else
            col_index_foc = col_index_foc + 1; % Adjust for header row
        end

        % Find or add the desired intensity row
        row_index_int = find(prestus_int == parameters.calibration.desired_intensity, 1);
        if isempty(row_index_int)
            row_index_int = size(virtual_data, 1) + 1;
            virtual_data{row_index_int, 1} = parameters.calibration.desired_intensity;
        else
            row_index_int = row_index_int + 1; % Adjust for header row
        end

        % Update optimized values
        virtual_data{row_index_int, col_index_foc} = mat2str([opt_phases, elem_p_amp]);

        % Handle missing values
        mask = cellfun(@(x) isempty(x) || isa(x, 'missing'), virtual_data);
        virtual_data(mask) = {[]};

        % Sort rows by desired intensity
        first_col = [virtual_data{2:end, 1}];
        [~, sortIdxCol] = sort(first_col);
        virtual_data_sorted = [virtual_data(1, :); virtual_data(1 + sortIdxCol, :)];

        % Sort columns by focal distance
        first_row = [virtual_data_sorted{1, 2:end}];
        [~, sortIdxRow] = sort(first_row);
        virtual_data_sorted = [virtual_data_sorted(:, 1), virtual_data_sorted(:, 1 + sortIdxRow)];

        % Write sorted data back to file
        writecell(virtual_data_sorted, output_file_path, 'FileType', 'text');
    else
        % Create a new file
        virtual_data = {
            'Desired Intensity [W/cm^2] - Focus wrt Exit Plane [mm]', parameters.calibration.desired_focal_distance_ep;
            parameters.calibration.desired_intensity, mat2str([opt_phases, elem_p_amp])
        };

        writecell(virtual_data, output_file_path, 'FileType', 'text');
    end

    % Save optimized transducer parameters to a YAML file for easy 
    % integration in PRESTUS config file
    if mod(parameters.calibration.desired_focal_distance_ep,1) == 0
        yaml_file = sprintf('%s-F%.0fmm-I%.0fwpercm2.yaml', parameters.calibration.equipment_name, ...
            parameters.calibration.desired_focal_distance_ep, parameters.calibration.desired_intensity);
    else
        yaml_file = sprintf('%s-F%.1fmm-I%.0fwpercm2.yaml', parameters.calibration.equipment_name, ...
            parameters.calibration.desired_focal_distance_ep, parameters.calibration.desired_intensity);
    end
    yaml_path = fullfile(parameters.calibration.path_output_profiles, yaml_file);

    parameters.transducer.focal_distance_ep = parameters.calibration.desired_focal_distance_ep;
    parameters.transducer.tx_calibration.isppa_wcm2   = parameters.calibration.desired_intensity;
    parameters.transducer.tx_calibration.phase_source = 'library';
    parameters.transducer.tx_calibration.isppa_source = 'library';

    data = struct();
    data.meta.tran_serial = parameters.calibration.equipment_name;
    data.meta.created_at  = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    try
        data.meta.prestus_version = prestus_version();
    catch
        data.meta.prestus_version = 'unknown';
    end
    % Calibration provenance
    cal = parameters.calibration;
    if isfield(cal, 'phase_mode');       data.meta.calibration.phase_mode       = cal.phase_mode;       end
    if isfield(cal, 'forward_model');    data.meta.calibration.forward_model    = cal.forward_model;    end
    if isfield(cal, 'opt_regularization_lambda')
        data.meta.calibration.lambda = cal.opt_regularization_lambda;
    end
    if isfield(cal, 'n_depths_used');    data.meta.calibration.n_depths_used    = cal.n_depths_used;    end
    if isfield(cal, 'min_fit_error');    data.meta.calibration.min_fit_error    = cal.min_fit_error;    end
    if isfield(cal, 'opt_velocity');     data.meta.calibration.opt_velocity     = cal.opt_velocity;     end
    if isfield(cal, 'calibration_type'); data.meta.calibration.calibration_type = cal.calibration_type; end
    data.transducer = parameters.transducer;
    yaml.dumpFile(yaml_path, data);

    fprintf('Transducer parameters with optimized values saved to YAML: %s \n', yaml_path);

end