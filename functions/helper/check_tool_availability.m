function checks = check_tool_availability(parameters, options)
% CHECK_TOOL_AVAILABILITY  Check that the external tools the requested stages need exist
%
% Derives the required tools from the enabled modules and settings, checks
% each one, and stops with a single error that lists every missing tool. A
% stage that would otherwise fail hours into a job (or fail silently) fails
% before any computation starts.
%
% Required tools, by condition:
%   JVM                    always (YAML loader, run_command, job submission)
%   sbatch / qsub          Context "submit" with platform slurm / qsub
%   SimNIBS charm          layered medium without final_tissues.nii.gz in the
%                          m2m folder, or io.overwrite_simnibs set
%   SimNIBS subject2mni    NIfTI stage with io.save_MNI, layered medium and
%                          analysis.mni_warp_method 'simnibs'; also needs m2m/toMNI
%   SimNIBS mni2subject_coords  transducer placement from MNI coordinates
%                          (placement.mode 'mni', or 'manual' with
%                          placement.manual.coordinate_space 'mni'); needs m2m/toMNI
%   k-Wave binary          acoustic stages with code_type cpp_gpu
%                          (kspaceFirstOrder-CUDA), cpp_cpu (kspaceFirstOrder-OMP)
%                          or cpp_mlx (kspaceFirstOrder-MLX)
%   GPU                    Context "run", acoustic stages with code_type
%                          matlab_gpu or cpp_gpu. Another process holding more
%                          than half of a GPU's memory fails the check; more
%                          than 5% raises a warning (prestus:gpuInUse).
%
% Use as:
%   check_tool_availability(parameters)                     % inside the job
%   check_tool_availability(parameters, Context="submit")  % before HPC submission
%   checks = check_tool_availability(...)
%
% Input:
%   parameters      - PRESTUS config (after load_parameters)
%   options.Context - "run" (default): the node that runs the simulation;
%                     checks everything. "submit": the node that submits the
%                     job; skips the GPU check and adds the scheduler check.
%
% Output:
%   checks - struct array with fields name, ok, detail (one entry per check)
%
% See also: PRESTUS_PIPELINE, PRESTUS_PIPELINE_START, RUN_COMMAND

    arguments
        parameters (1,1) struct
        options.Context (1,1) string {mustBeMember(options.Context, ["run", "submit"])} = "run"
    end

    checks = struct('name', {}, 'ok', {}, 'detail', {});
    gpu_warnings = {};

    % ---- JVM ----
    checks(end+1) = make_check('JVM', usejava('jvm') == 1, ...
        'MATLAB started without a JVM (needed by the YAML loader, run_command and job submission). Start MATLAB with Java, e.g. another MATLAB version or node.');

    % ---- scheduler (submitting node only) ----
    platform = get_field(parameters, {'platform'}, 'matlab');
    if options.Context == "submit"
        switch platform
            case 'slurm'
                checks(end+1) = make_check('sbatch', ~isempty(find_on_path('sbatch')), ...
                    'sbatch is not on PATH; jobs cannot be submitted from this node.');
            case 'qsub'
                checks(end+1) = make_check('qsub', ~isempty(find_on_path('qsub')), ...
                    'qsub is not on PATH; jobs cannot be submitted from this node.');
        end
    end

    % ---- SimNIBS ----
    medium      = char(get_field(parameters, {'simulation', 'medium'}, ''));
    is_layered  = strcmp(medium, 'layered');
    simnibs_bin = char(get_field(parameters, {'startup', 'simnibs_bin_path'}, ''));
    m2m = '';
    seg_root = char(get_field(parameters, {'path', 'seg'}, ''));
    if ~isempty(seg_root) && isfield(parameters, 'subject_id')
        m2m = fullfile(seg_root, sprintf('m2m_sub-%03d', parameters.subject_id));
    end

    simnibs_tools = {};   % {tool, reason, needs_toMNI}
    if is_layered && ~isempty(m2m)
        seg_missing = ~isfile(fullfile(m2m, 'final_tissues.nii.gz'));
        if seg_missing || logical(get_field(parameters, {'io', 'overwrite_simnibs'}, false))
            simnibs_tools(end+1, :) = {'charm', 'segmentation', false};
        end
    end
    nifti_on  = module_on(parameters, 'run_nifti_creation');
    mni_on    = logical(get_field(parameters, {'io', 'save_MNI'}, true));
    warp      = char(get_field(parameters, {'analysis', 'mni_warp_method'}, 'simnibs'));
    if nifti_on && mni_on && is_layered && strcmpi(warp, 'simnibs')
        simnibs_tools(end+1, :) = {'subject2mni', 'MNI conversion (io.save_MNI)', true};
    end
    if module_on(parameters, 'run_transducer_placement')
        pmode = lower(char(get_field(parameters, {'placement', 'mode'}, '')));
        cspace = lower(char(get_field(parameters, {'placement', 'manual', 'coordinate_space'}, '')));
        if strcmp(pmode, 'mni') || (strcmp(pmode, 'manual') && strcmp(cspace, 'mni'))
            simnibs_tools(end+1, :) = {'mni2subject_coords', 'placement from MNI coordinates', true};
        end
    end

    needs_toMNI = false;
    for k = 1:size(simnibs_tools, 1)
        tool = simnibs_tools{k, 1};
        name = ['SimNIBS ' tool];
        if isempty(simnibs_bin)
            checks(end+1) = make_check(name, false, sprintf( ...
                'startup.simnibs_bin_path is empty; %s needs %s.', simnibs_tools{k, 2}, tool)); %#ok<AGROW>
        else
            checks(end+1) = make_check(name, isfile(fullfile(simnibs_bin, tool)), sprintf( ...
                '%s not found in %s (startup.simnibs_bin_path); %s needs it.', ...
                tool, simnibs_bin, simnibs_tools{k, 2})); %#ok<AGROW>
        end
        needs_toMNI = needs_toMNI || simnibs_tools{k, 3};
    end
    if needs_toMNI && ~isempty(m2m)
        checks(end+1) = make_check('m2m toMNI', isfolder(fullfile(m2m, 'toMNI')), sprintf( ...
            'No toMNI folder in %s; MNI transforms need the SimNIBS segmentation warps.', m2m));
    end

    % ---- k-Wave binary and GPU ----
    acoustic_on = module_on(parameters, 'run_acoustic_sims') || ...
                  module_on(parameters, 'run_water_baseline') || ...
                  module_on(parameters, 'run_posthoc_water_sims');
    code_type = char(get_field(parameters, {'simulation', 'code_type'}, ''));
    if acoustic_on
        binaries = struct('cpp_gpu', 'kspaceFirstOrder-CUDA', ...
                          'cpp_cpu', 'kspaceFirstOrder-OMP', ...
                          'cpp_mlx', 'kspaceFirstOrder-MLX');
        if isfield(binaries, code_type)
            bin_name = binaries.(code_type);
            if ispc && ~strcmp(code_type, 'cpp_mlx')
                bin_name = [bin_name '.exe'];
            end
            if exist('getkWavePath', 'file')
                bin_path = fullfile(getkWavePath('binaries'), bin_name);
                checks(end+1) = make_check(bin_name, isfile(bin_path), sprintf( ...
                    'code_type %s needs %s, which is missing.', code_type, bin_path));
            else
                checks(end+1) = make_check(bin_name, false, ...
                    'k-Wave is not on the MATLAB path (getkWavePath not found).');
            end
        end

        if options.Context == "run" && any(strcmp(code_type, {'matlab_gpu', 'cpp_gpu'}))
            [gpu_check, gpu_warnings] = check_gpu(code_type);
            checks(end+1) = gpu_check;
        end
    end

    % ---- report ----
    for k = 1:numel(gpu_warnings)
        warning('prestus:gpuInUse', '%s', gpu_warnings{k});
    end
    failed = checks(~[checks.ok]);
    if ~isempty(failed)
        lines = arrayfun(@(c) sprintf('  - %s: %s', c.name, c.detail), failed, 'UniformOutput', false);
        error('prestus:missingTools', ...
            'Required tools are not available (context: %s):\n%s', ...
            options.Context, strjoin(lines, newline));
    end
    fprintf('[tool check] %s: %s ok\n', options.Context, strjoin({checks.name}, ', '));
end

% -------------------------------------------------------------------------
function [chk, warnings_out] = check_gpu(code_type)
% A visible GPU is required. nvidia-smi reports memory without creating a
% CUDA context in this process; memory already in use means another process
% shares the GPU.
    warnings_out = {};
    smi = find_on_path('nvidia-smi');
    if ~isempty(smi)
        [status, out] = run_command({smi, '--query-gpu=index,memory.used,memory.total', ...
            '--format=csv,noheader,nounits'});
        rows = strsplit(strtrim(out), newline);
        rows = rows(~cellfun(@isempty, rows));
        if status ~= 0 || isempty(rows)
            chk = make_check('GPU', false, sprintf( ...
                'code_type %s needs a GPU, but nvidia-smi found none (exit %d): %s', ...
                code_type, status, strtrim(out)));
            return
        end
        blocked = {};
        for r = 1:numel(rows)
            v = str2double(strsplit(rows{r}, ','));
            if numel(v) ~= 3 || any(isnan(v)); continue; end
            if v(2) > 0.5 * v(3)
                blocked{end+1} = sprintf('GPU %d has %.1f of %.1f GB in use by another process', ...
                    v(1), v(2)/1024, v(3)/1024); %#ok<AGROW>
            elseif v(2) > 0.05 * v(3)
                warnings_out{end+1} = sprintf( ...
                    ['GPU %d already has %.1f of %.1f GB in use by another process. ' ...
                     'The simulation may run out of GPU memory.'], v(1), v(2)/1024, v(3)/1024); %#ok<AGROW>
            end
        end
        if ~isempty(blocked)
            chk = make_check('GPU', false, sprintf( ...
                ['%s on node %s. Exclude the node (hpc.exclude) or ask the cluster ' ...
                 'admins to clear the GPU.'], strjoin(blocked, '; '), host_name()));
        else
            chk = make_check('GPU', true, '');
        end
        return
    end
    % No nvidia-smi: fall back to MATLAB's device count (Parallel Computing Toolbox).
    if exist('gpuDeviceCount', 'file') && gpuDeviceCount > 0
        chk = make_check('GPU', true, '');
    else
        chk = make_check('GPU', false, sprintf( ...
            'code_type %s needs a GPU, but neither nvidia-smi nor gpuDeviceCount found one.', code_type));
    end
end

function h = host_name()
    h = strtrim(getenv('HOSTNAME'));
    if isempty(h)
        h = strtrim(getenv('SLURMD_NODENAME'));
    end
    if isempty(h)
        h = 'unknown';
    end
end

function c = make_check(name, ok, detail)
    c = struct('name', name, 'ok', logical(ok), 'detail', detail);
end

function tf = module_on(parameters, name)
% Modules missing from the config count as enabled, so the check errs on the
% side of requiring the tool. Text values such as 'auto' (run_water_baseline)
% count as enabled unless they say off; logical() on a char array returned a
% vector, which made the || in the caller fail when run_acoustic_sims was 0.
    v = get_field(parameters, {'modules', name}, true);
    if ischar(v) || isstring(v)
        tf = ~any(strcmpi(strtrim(char(v)), {'0', 'off', 'false', 'no', 'never'}));
    else
        tf = any(logical(v(:)));
    end
end

function v = get_field(s, path, default)
    v = default;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}) || isempty(s.(path{k}))
            return
        end
        s = s.(path{k});
    end
    v = s;
end

function p = find_on_path(name)
    p = '';
    dirs = strsplit(getenv('PATH'), pathsep);
    for k = 1:numel(dirs)
        cand = fullfile(dirs{k}, name);
        if ~isempty(dirs{k}) && isfile(cand)
            p = cand;
            return
        end
    end
end
