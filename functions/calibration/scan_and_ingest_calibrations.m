function scan_and_ingest_calibrations(calibration_output_folder, library_path)
% SCAN_AND_INGEST_CALIBRATIONS  Scan a folder and ingest all calibration YAMLs
%
% Finds all *-F*mm-I*wpercm2.yaml files in calibration_output_folder, groups
% them by combo_name (the prefix before -F), and calls
% update_transducer_library for each combo.
%
% Use as:
%   scan_and_ingest_calibrations()
%   scan_and_ingest_calibrations(calibration_output_folder)
%   scan_and_ingest_calibrations(calibration_output_folder, library_path)
%
% Input:
%   calibration_output_folder - folder containing per-run YAML files
%                               (default: current directory)
%   library_path              - folder for the transducer library
%                               (default: config/transducer/ under PRESTUS root)
%
% Output:
%   (none) — calls update_transducer_library for each combo found
%
% See also: UPDATE_TRANSDUCER_LIBRARY, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 1 || isempty(calibration_output_folder)
        calibration_output_folder = pwd();
    end
    if nargin < 2 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'transducer');
    end

    %% Scan for all calibration files
    files = dir(fullfile(calibration_output_folder, '*-F*mm-I*wpercm2.yaml'));
    if isempty(files)
        fprintf('scan_and_ingest_calibrations: no calibration files found in:\n  %s\n', ...
            calibration_output_folder);
        return;
    end

    %% Group by combo_name (prefix before -F)
    combo_map = struct();
    for fi = 1:numel(files)
        fname = files(fi).name;
        tok   = regexp(fname, '^(.+)-F[0-9.]+mm-I[0-9.]+wpercm2\.yaml$', 'tokens', 'once');
        if isempty(tok)
            warning('scan_and_ingest_calibrations: cannot parse combo from ''%s'', skipping.', fname);
            continue;
        end
        key = matlab.lang.makeValidName(tok{1});
        if ~isfield(combo_map, key)
            combo_map.(key) = tok{1};
        end
    end

    combo_keys = fieldnames(combo_map);
    if isempty(combo_keys)
        fprintf('scan_and_ingest_calibrations: no valid combo files found.\n');
        return;
    end

    fprintf('scan_and_ingest_calibrations: found %d combo(s) in %s\n', ...
        numel(combo_keys), calibration_output_folder);

    %% Process each combo
    for ci = 1:numel(combo_keys)
        combo_name = combo_map.(combo_keys{ci});

        % Count how many depths for this combo
        depth_files = dir(fullfile(calibration_output_folder, [combo_name '-F*mm-I*wpercm2.yaml']));
        depth_toks  = regexp({depth_files.name}, '-F([0-9.]+)mm-I', 'tokens', 'once');
        depths_str  = cellfun(@(t) t{1}, depth_toks(~cellfun(@isempty, depth_toks)), 'UniformOutput', false);
        unique_d    = unique(cellfun(@str2double, depths_str));

        fprintf('  [%d/%d] %s — %d depth(s): [%s] mm\n', ...
            ci, numel(combo_keys), combo_name, numel(unique_d), num2str(unique_d(:)', '%.0f '));

        update_transducer_library(combo_name, calibration_output_folder, library_path);
    end

    fprintf('scan_and_ingest_calibrations: done.\n');
end
