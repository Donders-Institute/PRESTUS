function save_optimized_values(prestus_path, focus_wrt_exit_plane, desired_intensity, opt_param, ds_serial, tran_serial, save_in_general_folder)
    % Save optimized phases and amplitude values to a CSV file.
    %
    % Arguments:
    % - prestus_path: Path to the PRESTUS file for saving optimized values.
    % - focus_wrt_exit_plane: Target focal distance with respect to the exit plane [mm].
    % - desired_intensity: Target intensity for optimization [W/cm^2].
    % - opt_param: Structure containing optimized parameters, including transducer values.
    % - ds_serial: Serial number of the data system.
    % - tran_serial: Serial number of the transducer.
    % - save_in_general_folder: Option to save data in the general PRESTUS output folder.
    
    disp('Saving optimized values to CSV file...');
    
    % Extract and round phases and amplitudes
    opt_phases = round(opt_param.transducer.source_phase_deg, 2);
    source_amp = double(opt_param.transducer.source_phase_deg(1));

    if save_in_general_folder
        fprintf('CSV file can be found here: %s \n', prestus_path);
        
        % Check if the CSV file exists
        if isfile(prestus_path)
            % Read existing data
            virtual_data = readcell(prestus_path);
            prestus_foci_wrt_exit_plane = round(cell2mat(virtual_data(1, 2:end)), 2);
            prestus_int = cell2mat(virtual_data(2:end, 1));
    
            % Find or add the focal distance column
            col_index_foc = find(prestus_foci_wrt_exit_plane == focus_wrt_exit_plane, 1);
            if isempty(col_index_foc)
                col_index_foc = size(virtual_data, 2) + 1;
                virtual_data{1, col_index_foc} = focus_wrt_exit_plane;
            else
                col_index_foc = col_index_foc + 1; % Adjust for header row
            end
    
            % Find or add the desired intensity row
            row_index_int = find(prestus_int == desired_intensity, 1);
            if isempty(row_index_int)
                row_index_int = size(virtual_data, 1) + 1;
                virtual_data{row_index_int, 1} = desired_intensity;
            else
                row_index_int = row_index_int + 1; % Adjust for header row
            end
    
            % Update optimized values
            virtual_data{row_index_int, col_index_foc} = mat2str([opt_phases, source_amp]);
    
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
            writecell(virtual_data_sorted, prestus_path, 'FileType', 'text');
        else
            % Create a new file
            virtual_data = {
                'Desired Intensity [W/cm^2] - Focus wrt Exit Plane [mm]', focus_wrt_exit_plane;
                desired_intensity, mat2str([opt_phases, source_amp])
            };
    
            writecell(virtual_data, prestus_path, 'FileType', 'text');
        end
    else
        % Create the local path
        [~, name, ext] = fileparts(prestus_path);
        local_path = fullfile(opt_param.output_location, strcat(name, ext));
        fprintf('CSV file can be found here: %s \n', local_path);
        
        % Create a new file
        virtual_data = {
            'Desired Intensity [W/cm^2] - Focus wrt Exit Plane [mm]', focus_wrt_exit_plane;
            desired_intensity, mat2str([opt_phases, source_amp])
        };

        writecell(virtual_data, local_path, 'FileType', 'text');
    end

    % Save optimized transducer parameters to a YAML file for easy 
    % integration in PRESTUS config file
    yaml_file = sprintf('%s-%s-F%.0fmm-I%.0fwpercm2.yaml', tran_serial, ds_serial, focus_wrt_exit_plane, desired_intensity);
    yaml_path = fullfile(opt_param.output_location, yaml_file);

    opt_param.transducer.set_focus_wrt_exit_plane_mm = focus_wrt_exit_plane;
    opt_param.transducer.set_intensity_w_per_cm2 = desired_intensity;

    % Wrap transducer parameters in a parent structure for YAML
    data = struct('transducer', opt_param.transducer);
    yaml.dumpFile(yaml_path, data);

    fprintf('Transducer parameters with optimized values saved to YAML: %s \n', yaml_path);

end