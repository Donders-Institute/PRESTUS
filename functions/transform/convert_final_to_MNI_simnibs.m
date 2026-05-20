function convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters, options)
% CONVERT_FINAL_TO_MNI_SIMNIBS  Transform one or more images to MNI space using SimNIBS
%
% Calls the SimNIBS subject2mni CLI tool to warp NIfTI image(s) from
% subject-specific (T1 conform) space to MNI space. An optional
% LD_LIBRARY_PATH export is prepended when parameters.hpc.ld_library_path
% is set. Because SimNIBS appends '_MNI' to the output filename by default,
% each result is renamed to the corresponding path_to_output_img.
%
% When called with cell-array inputs (batch mode) on Unix/macOS, all
% subject2mni processes are launched simultaneously as background jobs and
% the function waits until every output file appears (timeout: 1 hour).
% On Windows or for a single file, execution is sequential.
%
% Use as:
%   convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
%   convert_final_to_MNI_simnibs({in1,in2,...}, m2m_folder, {out1,out2,...}, parameters)
%   convert_final_to_MNI_simnibs(..., 'interpolation_order', 0)
%
% Input:
%   path_to_input_img  - path or cell array of paths to input NIfTI(s) in subject space
%   m2m_folder         - path to SimNIBS m2m folder with registration data (toMNI subdirectory)
%   path_to_output_img - desired output path(s) in MNI space (.nii.gz)
%   parameters         - PRESTUS config; uses startup.simnibs_bin_path and hpc.ld_library_path
%   options.interpolation_order - resampling order passed to subject2mni
%                                 (0=nearest, 1=linear; default: 1)
%
% See also: SIMULATION_NIFTI, CONVERT_FINAL_TO_MNI_MATLAB, SUBJECT2MNI_COORDS_LDFIX

    arguments
        path_to_input_img
        m2m_folder string
        path_to_output_img
        parameters struct
        options.interpolation_order = 1
    end

    % Normalise scalar inputs to cell arrays so the rest of the code is uniform
    if ischar(path_to_input_img) || isstring(path_to_input_img)
        input_list  = {char(path_to_input_img)};
        output_list = {char(path_to_output_img)};
    else
        input_list  = path_to_input_img;
        output_list = path_to_output_img;
    end
    n = numel(input_list);

    if isfield(parameters.hpc, 'ld_library_path') && ~isempty(parameters.hpc.ld_library_path) && ~ispc
        ld_command = sprintf('export LD_LIBRARY_PATH="%s"; ', parameters.hpc.ld_library_path);
    else
        ld_command = '';
    end
    subject2mni_bin = fullfile(parameters.startup.simnibs_bin_path, 'subject2mni');

    % SimNIBS appends '_MNI' to the output stem; track the actual filenames it writes
    simnibs_names = cell(1, n);
    for k = 1:n
        out = output_list{k};
        if ~matches(out, wildcardPattern + '_MNI.nii.gz')
            simnibs_names{k} = strrep(out, '.nii.gz', '_MNI.nii.gz');
        else
            simnibs_names{k} = out;
        end
    end

    if n > 1 && ~ispc
        % Batch mode: launch all conversions simultaneously as background processes
        for k = 1:n
            cmd = sprintf('%s"%s" --in "%s" --out "%s" --m2mpath "%s" --interpolation_order %d', ...
                ld_command, subject2mni_bin, input_list{k}, output_list{k}, ...
                m2m_folder, options.interpolation_order);
            system([cmd ' &']);
        end
        % Poll until every expected output file exists (max 1 hour)
        t0 = tic;
        while toc(t0) < 3600
            if all(cellfun(@isfile, simnibs_names)); break; end
            pause(5);
        end
        if ~all(cellfun(@isfile, simnibs_names))
            warning('PRESTUS:MNITimeout', ...
                'convert_final_to_MNI_simnibs: timed out waiting for subject2mni outputs.');
        end
    else
        % Single file or Windows: sequential
        for k = 1:n
            system(sprintf('%s"%s" --in "%s" --out "%s" --m2mpath "%s" --interpolation_order %d', ...
                ld_command, subject2mni_bin, input_list{k}, output_list{k}, ...
                m2m_folder, options.interpolation_order));
        end
    end

    % Rename each SimNIBS output to the requested destination path
    for k = 1:n
        if ~strcmp(simnibs_names{k}, output_list{k}) && isfile(simnibs_names{k})
            movefile(simnibs_names{k}, output_list{k});
        end
    end
end
