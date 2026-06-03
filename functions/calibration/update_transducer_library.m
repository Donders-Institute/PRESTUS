function update_transducer_library(combo_name, calibration_output_folder, library_path)
% UPDATE_TRANSDUCER_LIBRARY  Build or update the global model in the library YAML
%
% Scans calibration_output_folder for per-run YAML files matching
% {combo_name}-F*mm-I*wpercm2.yaml, groups them by depth, and builds or
% updates the global_model entry in the library YAML for this combo.
%
% Use as:
%   update_transducer_library(combo_name)
%   update_transducer_library(combo_name, calibration_output_folder)
%   update_transducer_library(combo_name, calibration_output_folder, library_path)
%
% Input:
%   combo_name                - equipment combination key, e.g.
%                               'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'
%   calibration_output_folder - folder containing per-run YAML files
%                               (default: current directory)
%   library_path              - folder for the transducer library
%                               (default: config/transducer/ under PRESTUS root)
%
% Output:
%   (none) — writes/updates {library_path}/{combo_name}.yaml
%
% See also: SCAN_AND_INGEST_CALIBRATIONS, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 2 || isempty(calibration_output_folder)
        calibration_output_folder = pwd();
    end
    if nargin < 3 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'transducer');
    end

    %% Scan for matching calibration files
    pattern = fullfile(calibration_output_folder, [combo_name '-F*mm-I*wpercm2.yaml']);
    files   = dir(pattern);
    if isempty(files)
        error('update_transducer_library: no calibration files found matching:\n  %s', pattern);
    end

    %% Parse each file: extract depth, intensity, phases, amplitude
    n_files    = numel(files);
    depths_all = zeros(1, n_files);
    intens_all = zeros(1, n_files);
    phases_all = cell(1, n_files);
    amps_all   = zeros(1, n_files);
    tran_serial_found = '';

    for fi = 1:n_files
        fname = files(fi).name;
        fpath = fullfile(calibration_output_folder, fname);
        d     = yaml.loadFile(fpath, 'ConvertToArray', true);

        % Parse depth and intensity from filename
        tok = regexp(fname, '-F([0-9.]+)mm-I([0-9.]+)wpercm2\.yaml$', 'tokens', 'once');
        if isempty(tok)
            warning('update_transducer_library: cannot parse depth/intensity from ''%s'', skipping.', fname);
            continue;
        end
        depths_all(fi) = str2double(tok{1});
        intens_all(fi) = str2double(tok{2});

        tran = d.transducer;
        phases_all{fi} = tran.annular.elem_phase_deg(:)';
        amps_all(fi)   = double(tran.annular.elem_amp(1));

        % Try to extract tran_serial from meta if present
        if isempty(tran_serial_found) && isfield(d, 'meta') && isfield(d.meta, 'tran_serial')
            tran_serial_found = d.meta.tran_serial;
        end
        % Fallback: name field in transducer or combo_name itself
        if isempty(tran_serial_found) && isfield(tran, 'name')
            tran_serial_found = tran.name;
        end
    end

    % Remove any skipped entries
    valid = depths_all > 0;
    depths_all = depths_all(valid);
    intens_all = intens_all(valid);
    phases_all = phases_all(valid);
    amps_all   = amps_all(valid);

    if isempty(depths_all)
        error('update_transducer_library: no valid calibration files found for ''%s''.', combo_name);
    end

    % Determine tran_serial: check combo_name for underscore pattern
    if isempty(tran_serial_found)
        if contains(combo_name, '_')
            % Try to match first segment against equipment config serials
            parts = strsplit(combo_name, '_');
            % Greedily try longest prefix that matches equipment
            tran_serial_found = combo_name; % fallback
            try
                eq = load_equipment_config();
                for p = numel(parts):-1:1
                    candidate = strjoin(parts(1:p), '_');
                    if isfield(eq.trans, candidate)
                        tran_serial_found = candidate;
                        break;
                    end
                end
            catch
                tran_serial_found = parts{1};
            end
        else
            tran_serial_found = combo_name;
        end
    end

    %% Group by depth: for each unique depth, select one intensity for phases
    unique_depths = unique(depths_all);
    n_depths      = numel(unique_depths);

    depth_phase_map = cell(1, n_depths);  % one phase array per depth
    amp_struct      = struct();           % intensity -> amplitude mapping

    for di = 1:n_depths
        d_val = unique_depths(di);
        mask  = depths_all == d_val;
        % For phase: use first file at this depth (intensity doesn't change phases)
        idx_list = find(mask);
        depth_phase_map{di} = phases_all{idx_list(1)};
    end

    % Amplitude scaling: all intensity/amplitude pairs across all files
    all_intens_unique = unique(intens_all);
    for ii = 1:numel(all_intens_unique)
        i_val = all_intens_unique(ii);
        mask  = intens_all == i_val;
        % Average amplitude if multiple files at same intensity (should be same)
        amp_val = round(mean(amps_all(mask)));
        key     = sprintf('i%s', strrep(num2str(i_val), '.', 'p'));
        amp_struct.(key) = amp_val;
    end

    %% Monotonicity check
    if numel(unique_depths) > 1
        n_elem = numel(depth_phase_map{1});
        for elem_i = 1:n_elem
            elem_phases = cellfun(@(p) p(elem_i), depth_phase_map);
            diffs = diff(elem_phases);
            if ~all(diffs >= 0) && ~all(diffs <= 0)
                warning('update_transducer_library:nonMonotonic', ...
                    'Element %d phases are non-monotonic across depths [%s] deg. ' ...
                    'Check calibration quality.', ...
                    elem_i, num2str(elem_phases, '%.1f '));
            end
        end
    end

    %% Load existing library YAML (if present), merge
    if ~exist(library_path, 'dir')
        mkdir(library_path);
    end
    lib_yaml = fullfile(library_path, [combo_name '.yaml']);

    if isfile(lib_yaml)
        lib = yaml.loadFile(lib_yaml, 'ConvertToArray', true);
    else
        lib = struct();
    end

    %% Build/update meta
    today_str = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    if ~isfield(lib, 'meta') || isempty(lib.meta)
        lib.meta = struct();
        lib.meta.created_at = today_str;
    end
    lib.meta.tran_serial   = tran_serial_found;
    lib.meta.last_updated  = today_str;
    try
        [~, ver_str]            = prestus_version();
        lib.meta.prestus_version = ver_str;
    catch
        lib.meta.prestus_version = 'unknown';
    end

    %% Merge depths into existing global_model (if any)
    if isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm') && ...
            ~isempty(lib.global_model.depths_ep_mm)
        existing_depths  = lib.global_model.depths_ep_mm(:)';
        existing_phases  = lib.global_model.elem_phase_deg;  % cell array
        if ~iscell(existing_phases)
            existing_phases = num2cell(existing_phases, 2);
        end

        % Merge: update existing depths, append new ones
        all_depths   = existing_depths;
        all_phases   = existing_phases;
        for di = 1:n_depths
            d_val = unique_depths(di);
            [~, match] = min(abs(all_depths - d_val));
            if abs(all_depths(match) - d_val) < 0.5
                all_phases{match} = depth_phase_map{di};
            else
                all_depths(end+1)   = d_val; %#ok<AGROW>
                all_phases{end+1}   = depth_phase_map{di}; %#ok<AGROW>
            end
        end
        % Sort by depth
        [all_depths, sort_idx] = sort(all_depths);
        all_phases = all_phases(sort_idx);

        lib.global_model.depths_ep_mm  = all_depths;
        lib.global_model.elem_phase_deg = all_phases;

        % Merge amplitude_scaling
        if isfield(lib.global_model, 'amplitude_scaling')
            existing_amp = lib.global_model.amplitude_scaling;
            new_keys = fieldnames(amp_struct);
            for ki = 1:numel(new_keys)
                existing_amp.(new_keys{ki}) = amp_struct.(new_keys{ki});
            end
            lib.global_model.amplitude_scaling = existing_amp;
        else
            lib.global_model.amplitude_scaling = amp_struct;
        end
    else
        % Build fresh global_model
        lib.global_model.depths_ep_mm   = unique_depths;
        lib.global_model.elem_phase_deg = depth_phase_map;
        lib.global_model.amplitude_scaling = amp_struct;
        lib.global_model.source = 'update_transducer_library';
    end

    %% Write updated YAML
    yaml.dumpFile(lib_yaml, lib);
    fprintf('Library updated: %s\n  %d depths: [%s] mm\n', ...
        lib_yaml, numel(lib.global_model.depths_ep_mm), ...
        num2str(lib.global_model.depths_ep_mm(:)', '%.1f '));
end
