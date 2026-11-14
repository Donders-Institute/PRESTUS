function runs = prestus_study_start(study_config, options)
% PRESTUS_STUDY_START  Run the PRESTUS pipeline for every row of a study manifest
%
% A study is a study config (YAML) whose study.manifest points to a table
% with one row per run: a subject, an optional condition label, and any
% per-run config values (placement, target, depth, intensity, transducer,
% ...). See STUDY_READ_MANIFEST for the manifest format and
% STUDY_ROW_PARAMETERS for the merge order.
%
% Each run is started with PRESTUS_PIPELINE_START, so the platform of the
% run's config applies: 'matlab' runs the rows one after the other in this
% session; 'slurm'/'qsub' submit one job per row (or a job graph for
% uncertainty / multi-ISPPA / async rows). The per-subject automatic group
% report is switched off for study runs (modules.generate_group_report = 0,
% applied after the manifest columns; only 'Overrides' can re-enable it); instead the study report
% (PRESTUS_GROUP_REPORT_START on the study config) is built once at the end:
% directly on 'matlab', or as a dependent HPC job ('afterany', so it also
% runs when some subjects fail).
%
% Use as:
%   prestus_study_start(study_config)
%   prestus_study_start(study_config, 'DryRun', true)
%   runs = prestus_study_start(study_config, 'Subjects', [1 3], 'Conditions', "M1")
%   prestus_study_start(study_config, 'Rows', 1, 'Platform', 'matlab', ...
%       'Overrides', {'grid.resolution_mm', 2})              % debug one row
%
% Input:
%   study_config - path to the study YAML (must set study.manifest and path.sim)
%   options      - name-value pairs:
%     'Subjects'    numeric subject IDs to run (default: all included rows)
%     'Conditions'  condition labels to run (default: all)
%     'Rows'        manifest row numbers to run (default: all)
%     'DryRun'      load and check every row, print the plan, run nothing
%                   (default false)
%     'Platform'    override the platform of every row ('matlab','slurm','qsub')
%     'GroupReport' build the study report after the runs (default true)
%     'StopOnError' stop at the first failing row on 'matlab' (default false)
%     'Overrides'   N×2 cell {key, value} applied to every row after the
%                   manifest columns, e.g. {'grid.resolution_mm', 2} for a
%                   quick test run
%
% Output:
%   runs - table of the selected runs with the resolved settings, platform,
%          job_id and status ('planned', 'done', 'failed', 'submitted',
%          'config error')
%
% A plan / launch log is written to <path.sim>/group/<study>/log/.
%
% See also: STUDY_SETUP, STUDY_READ_MANIFEST, STUDY_ROW_PARAMETERS,
%           PRESTUS_PIPELINE_START, PRESTUS_GROUP_REPORT_START

arguments
    study_config {mustBeTextScalar}
    options.Subjects (1,:) double = []
    options.Conditions (1,:) string = string.empty(1, 0)
    options.Rows (1,:) double = []
    options.DryRun (1,1) logical = false
    options.Platform {mustBeTextScalar} = ''
    options.GroupReport (1,1) logical = true
    options.StopOnError (1,1) logical = false
    options.Overrides cell = cell(0, 2)
end

    study = study_setup(study_config);
    runs = study.runs(study.runs.include, :);
    if ~isempty(options.Subjects);   runs = runs(ismember(runs.subject_id, options.Subjects), :); end
    if ~isempty(options.Conditions); runs = runs(ismember(runs.condition, options.Conditions), :); end
    if ~isempty(options.Rows);       runs = runs(ismember(runs.row, options.Rows), :); end
    n = height(runs);

    fprintf('========================================\n');
    fprintf('PRESTUS STUDY: %s\n', study.name);
    fprintf('Manifest: %s (%d row(s), %d selected)\n', study.manifest_path, height(study.manifest), n);
    fprintf('Outputs:  %s\n', char(study.parameters.path.sim));
    fprintf('========================================\n\n');
    if n == 0
        fprintf('No runs selected — nothing to do.\n');
        return
    end

    % Study runs: one study report at the end instead of one per subject.
    extra = [{'modules.generate_group_report', 0}; options.Overrides];
    if ~isempty(options.Platform)
        extra(end+1, :) = {'platform', char(options.Platform)};
    end

    %% Load every row first: config errors surface before anything runs
    params = cell(n, 1);
    status = repmat("planned", n, 1);
    message = strings(n, 1);
    platform = strings(n, 1);
    [placement_mode, target, transducer, focal_distance_ep, target_isppa] = deal(strings(n, 1));
    for i = 1:n
        try
            p = study_row_parameters(study.config, runs(i, :), extra);
            assert(strcmp(p.io.output_affix, runs.affix(i)), 'prestus_study_start:affixMismatch', ...
                'Resolved affix "%s" differs from the planned "%s".', p.io.output_affix, runs.affix(i));
            if ~isfield(p, 'platform') || strcmp(p.platform, 'auto')
                try
                    p.platform = hpc_detect_system();
                catch ME_detect
                    % No scheduler here: fine for a dry run, an error for a real run
                    if ~options.DryRun; rethrow(ME_detect); end
                    p.platform = 'auto';
                end
            end
            params{i} = p;
            platform(i) = string(p.platform);
            d = study_run_descriptors(p);
            placement_mode(i) = d.placement_mode;
            target(i)         = d.target_label;
            transducer(i)     = d.transducer_name;
            focal_distance_ep(i) = num2str_list(d.focal_distance_ep_mm);
            target_isppa(i)   = num2str_list(d.target_isppa_wcm2);
        catch ME
            status(i) = "config error";
            message(i) = string(ME.message);
            fprintf(2, 'Row %d (%s): %s\n', runs.row(i), runs.label(i), ME.message);
        end
    end
    runs = addvars(runs, platform, placement_mode, target, transducer, focal_distance_ep, target_isppa, ...
        'NewVariableNames', {'platform', 'placement_mode', 'target', 'transducer', ...
                             'focal_distance_ep_mm', 'target_isppa_wcm2'});
    runs.job_id  = strings(n, 1);
    runs.status  = status;
    runs.message = message;

    print_plan(runs);
    log_dir = fullfile(study.out_dir, 'log');
    if ~isfolder(log_dir); mkdir(log_dir); end
    stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

    if options.DryRun
        write_log(runs, fullfile(log_dir, sprintf('plan_%s.tsv', stamp)));
        n_err = sum(runs.status == "config error");
        fprintf('\nDry run: %d run(s) ready, %d with config errors. Nothing was started.\n', n - n_err, n_err);
        return
    end
    if any(runs.status == "config error") && options.StopOnError
        write_log(runs, fullfile(log_dir, sprintf('launch_%s.tsv', stamp)));
        error('prestus_study_start:configErrors', 'Config errors in %d row(s); nothing was started.', ...
            sum(runs.status == "config error"));
    end

    %% Run / submit
    hpc_ids = {};
    for i = 1:n
        if runs.status(i) == "config error"; continue; end
        fprintf('\n──── Run %d/%d: %s (row %d, affix "%s") ────\n', i, n, runs.label(i), ...
            runs.row(i), runs.affix(i));
        try
            job_id = prestus_pipeline_start(params{i});
            if isempty(job_id)
                if strcmp(params{i}.platform, 'matlab')
                    runs.status(i) = "done";
                else
                    runs.status(i) = "submitted";   % job graph without a single final ID
                end
            else
                runs.status(i) = "submitted";
                runs.job_id(i) = string(num2str(job_id));
                hpc_ids{end+1} = job_id; %#ok<AGROW>
            end
            runs.platform(i) = string(params{i}.platform);
        catch ME
            runs.status(i) = "failed";
            runs.message(i) = string(ME.message);
            fprintf(2, 'Run %s failed: %s\n', runs.label(i), ME.message);
            if options.StopOnError
                write_log(runs, fullfile(log_dir, sprintf('launch_%s.tsv', stamp)));
                rethrow(ME);
            end
        end
    end
    params = []; %#ok<NASGU>

    %% Study report
    if options.GroupReport
        on_hpc = any(ismember(runs.platform, ["slurm", "qsub"]) & runs.status == "submitted");
        if on_hpc
            if isempty(study.config)
                warning('prestus_study_start:noReportJob', ...
                    'Study report job needs the study config path; run prestus_group_report_start later.');
            else
                submit_report_job(study, runs, hpc_ids, log_dir);
            end
        else
            try
                prestus_group_report_start(study.config);
            catch ME
                warning('prestus_study_start:reportFailed', 'Study report failed: %s', ME.message);
            end
        end
    end

    write_log(runs, fullfile(log_dir, sprintf('launch_%s.tsv', stamp)));
    fprintf('\nStudy %s: %d done, %d submitted, %d failed, %d config error(s).\n', study.name, ...
        sum(runs.status == "done"), sum(runs.status == "submitted"), ...
        sum(runs.status == "failed"), sum(runs.status == "config error"));
end

% =========================================================================
function print_plan(runs)
    fprintf('%-4s %-22s %-24s %-9s %-10s %-14s %-14s %-9s %-7s %s\n', 'row', 'run', 'affix', ...
        'platform', 'placement', 'target', 'transducer', 'depth_ep', 'Isppa', 'status');
    for i = 1:height(runs)
        fprintf('%-4d %-22s %-24s %-9s %-10s %-14s %-14s %-9s %-7s %s\n', runs.row(i), ...
            clip(runs.label(i), 22), clip(runs.affix(i), 24), runs.platform(i), ...
            clip(runs.placement_mode(i), 10), clip(runs.target(i), 14), clip(runs.transducer(i), 14), ...
            runs.focal_distance_ep_mm(i), runs.target_isppa_wcm2(i), runs.status(i));
    end
end

function s = clip(s, w)
    s = char(s);
    if numel(s) > w; s = [s(1:w-1) '~']; end
end

function s = num2str_list(v)
    if isempty(v) || ~isnumeric(v)
        s = "";
    else
        s = string(strjoin(compose('%.4g', v(:)'), ' '));
    end
end

function write_log(runs, file)
    t = removevars(runs, {'overrides', 'config'});
    writetable(t, file, 'FileType', 'text', 'Delimiter', '\t');
    fprintf('Study log: %s\n', file);
end

function submit_report_job(study, runs, hpc_ids, log_dir)
% Dependent CPU job that builds the study report after all subject jobs.
    p = study.parameters;
    platform = char(runs.platform(find(ismember(runs.platform, ["slurm", "qsub"]), 1)));
    ids = cellfun(@(x) str2double(string(x)), hpc_ids);
    ids = ids(isfinite(ids));
    if isempty(ids)
        fprintf(['No subject job IDs to wait for; run prestus_group_report_start(''%s'') ' ...
                 'once the jobs have finished.\n'], study.config);
        return
    end
    if numel(ids) < sum(runs.status == "submitted")
        fprintf(['Note: some rows submit job graphs without a single final job ID; the study ' ...
                 'report job may start before they finish. Re-run prestus_group_report_start then.\n']);
    end

    p.subject_id = 0;
    p.io.output_affix = '';
    p.hpc.job_name = regexprep(['PS_grp_' study.name], '^(.{1,20}).*$', '$1');
    p.hpc.timelimit = '01:00:00';
    p.hpc.memorylimit = 16;
    p.hpc.depend_job_ids = ids;
    p.hpc.depend_type = 'afterany';
    % The report job runs no solver: no GPU request (hpc_submit_job derives
    % it from the solvers and their compute backends).
    p.simulation.code_type = 'matlab_cpu';
    p.acoustic_backend = 'kwave';
    p.thermal_backend = 'kwave';
    % Snellius jobs must name one of its GPU partitions (hpc_submit_job
    % derives cores/memory/time limits from it), so keep it there.
    if ~(isfield(p.hpc, 'name') && strcmp(p.hpc.name, 'snellius'))
        for f = {'gpu', 'n_gpu', 'partition'}
            if isfield(p.hpc, f{1}); p.hpc = rmfield(p.hpc, f{1}); end
        end
    end

    hpc_dir = fullfile(study.out_dir, 'log_hpc');
    if ~isfolder(hpc_dir); mkdir(hpc_dir); end
    prestus_path = get_prestus_path;
    m_file = sprintf('study_report_%s', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    fid = fopen(fullfile(hpc_dir, [m_file '.m']), 'w');
    fprintf(fid, 'addpath(''%s'');\n', matlab_escape_quotes(fullfile(prestus_path, 'functions', 'helper')));
    fprintf(fid, 'safe_addpath(''%s'');\n', matlab_escape_quotes(fullfile(prestus_path, 'functions')));
    fprintf(fid, 'safe_addpath(''%s'');\n', matlab_escape_quotes(fullfile(prestus_path, 'external')));
    fprintf(fid, 'prestus_group_report_start(''%s'');\n', matlab_escape_quotes(strrep(study.config, '''', '''''')));
    fclose(fid);

    job_id = hpc_submit_job(platform, m_file, p, hpc_dir);
    fprintf('Study report job %s submitted (after %d subject job(s)); logs in %s\n', ...
        num2str(job_id), numel(ids), hpc_dir);
    fprintf('Study log dir: %s\n', log_dir);
end
