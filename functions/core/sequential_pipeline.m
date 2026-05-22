function sequential_pipeline(parameters, options)
% SEQUENTIAL_PIPELINE  Dispatch follow-up simulations that chain thermal state
%
% After a base PRESTUS run completes, call this function to pop the next config
% from options.sequential_configs, wire up the adopted heatmap / CEM43 files
% from the just-completed run, set sensible cache-affix defaults, and recurse
% via prestus_pipeline_start.
%
% Affix defaults applied for each follow-up run:
%   output_affix        — base output_affix + '_seq<N>' (drives all output filenames)
%   preproc_affix       — reused from base run (head geometry unchanged)
%   acoustic_cache_affix — NOT defaulted; fresh run unless user sets it
%   thermal_cache_affix — base output_affix + '_seq<N>' to avoid "already done"
%
% A multi-run summary HTML report is generated automatically after the last
% sequential run completes (requires generate_sequential_report on the path).
%
% Intermediate data cleanup:
%   Set options.sequential_cleanup_intermediate = true (default: false) to
%   delete each run's per-run NIfTIs (nii/) and images (img/) after the
%   integrated sequential report has been generated.  Cache files (including
%   heating timeseries .mat) and top-level outputs (reports, CSVs) are kept.
%
% Use as:
%   sequential_pipeline(parameters, options)
%
% where parameters is the just-completed run's parameters struct and options
% contains at minimum options.sequential_configs.
%
% See also: PRESTUS_PIPELINE_START, GENERATE_SEQUENTIAL_REPORT

arguments
    parameters (1,1) struct
    options    (1,1) struct
end

    sequential_configs = options.sequential_configs;
    fields   = fieldnames(sequential_configs);
    numbers  = cellfun(@(x) sscanf(x, 'config_%d'), fields);
    [~, minIdx]   = min(numbers);
    lowestField   = fields{minIdx};

    sequential_parameters = sequential_configs.(lowestField);
    sequential_configs    = rmfield(sequential_configs, lowestField);

    % ---- wire up adopted thermal maps from the completed run ----
    sequential_parameters.io.adopted_heatmap = fullfile(parameters.io.dir_nii_T1w, ...
        sprintf('sub-%03d_%s_T1w%s_heating_end.nii.gz', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    sequential_parameters.io.adopted_cem43 = fullfile(parameters.io.dir_nii_T1w, ...
        sprintf('sub-%03d_%s_T1w%s_CEM43_end.nii.gz', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    sequential_parameters.io.adopted_cem43_iso = fullfile(parameters.io.dir_nii_T1w, ...
        sprintf('sub-%03d_%s_T1w%s_CEM43_iso_end.nii.gz', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    % ---- affix defaults ----
    seq_suffix = ['_seq', num2str(min(numbers))];

    % output_affix drives all output filenames; default to base + seq index
    % so follow-up runs get their own NIfTIs, reports, and CSVs.
    if ~isfield(sequential_parameters.io, 'output_affix') || ...
            isempty(sequential_parameters.io.output_affix)
        sequential_parameters.io.output_affix = [parameters.io.output_affix, seq_suffix];
    end

    if ~isfield(sequential_parameters.io, 'preproc_affix')
        sequential_parameters.io.preproc_affix = parameters.io.output_affix;
    end
    % acoustic_cache_affix intentionally NOT defaulted: fresh run by default.
    % Set it explicitly in the sequential config to reuse prior acoustics.
    if ~isfield(sequential_parameters.io, 'thermal_cache_affix')
        sequential_parameters.io.thermal_cache_affix = ...
            [parameters.io.output_affix, seq_suffix];
    end

    % ---- accumulate run configs for the sequential report ----
    if ~isfield(options, 'sequential_report_runs')
        all_run_params = [{parameters}, {sequential_parameters}];
    else
        all_run_params = [options.sequential_report_runs, {sequential_parameters}];
    end

    fprintf('Running subsequent heating simulation on %s\n', lowestField);

    is_last_sequential = isempty(fieldnames(sequential_configs));
    if ~is_last_sequential
        options.sequential_configs     = sequential_configs;
        options.sequential_report_runs = all_run_params;
    else
        options = rmfield(options, 'sequential_configs');
        if isfield(options, 'sequential_report_runs')
            options = rmfield(options, 'sequential_report_runs');
        end
    end

    prestus_pipeline_start(sequential_parameters, options);

    % ---- generate multi-run summary after the chain completes ----
    if is_last_sequential && numel(all_run_params) > 1
        report_ok = false;
        try
            generate_sequential_report(all_run_params);
            report_ok = true;
        catch ME_rep
            warning('prestus_pipeline:sequentialReport', ...
                'Sequential report generation failed: %s', ME_rep.message);
        end

        % ---- optional per-run NIfTI / image cleanup ----
        % Only runs after a successful report so integrated outputs exist first.
        % Cache (including heating timeseries .mat) is always retained.
        do_cleanup = report_ok && ...
                     isfield(options, 'sequential_cleanup_intermediate') && ...
                     options.sequential_cleanup_intermediate;
        if do_cleanup
            for ri = 1:numel(all_run_params)
                p = all_run_params{ri};
                if isfield(p.io, 'dir_output')
                    base = p.io.dir_output;
                else
                    continue
                end
                for subdir = {fullfile(base, 'nii'), fullfile(base, 'img')}
                    d = subdir{1};
                    if isfolder(d)
                        fprintf('Removing intermediate outputs: %s\n', d);
                        rmdir(d, 's');
                    end
                end
            end
        end
    end
end
