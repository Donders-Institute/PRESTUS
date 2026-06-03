function report = hpc_job_report(job_id, hpc_type, parameters, log_dir, carbon_params, cost_params)
% HPC_JOB_REPORT  Query scheduler for post-job resource usage and carbon metrics
%
% After a job finishes, queries sacct (SLURM) for actual CPU time, wall time,
% peak memory (MaxRSS), and GPU allocation, then estimates the CO2e footprint
% using the Green Algorithms framework and the monetary cost using cluster
% billing rates. Prints a formatted summary and writes a report file to
% log_dir for use in tuning future simulations.
%
% Use as:
%   report = hpc_job_report(job_id, hpc_type, parameters, log_dir)
%   report = hpc_job_report(job_id, hpc_type, parameters, log_dir, carbon_params)
%   report = hpc_job_report(job_id, hpc_type, parameters, log_dir, carbon_params, cost_params)
%
% Input:
%   job_id        - scheduler job ID (numeric for SLURM, char for qsub)
%   hpc_type      - 'slurm' or 'qsub'
%   parameters    - PRESTUS config struct; used for requested resource fields:
%                     .hpc.memorylimit  — requested memory in GB
%                     .hpc.timelimit    — requested wall time string
%                     .subject_id       — used in report file name
%   log_dir       - directory where the report file is written
%   carbon_params - (optional) struct overriding emission constants:
%                     .pue                  — Power Usage Effectiveness (default 1.67)
%                     .carbon_intensity     — gCO2e/kWh (default 270, Netherlands)
%                     .tdp_per_core_w       — CPU TDP per core in W (default 12)
%                     .tdp_per_gpu_w        — GPU TDP in W (default 300)
%                     .mem_power_w_per_gb   — W per GB RAM (default 0.3725)
%   cost_params   - (optional) struct overriding billing rates:
%                     .cpu_sbu_per_core_hour — SBUs per CPU core-hour (default 1)
%                     .gpu_sbu_per_gpu_hour  — SBUs per GPU-hour; auto-set from
%                                              parameters.hpc.partition for Snellius:
%                                              gpu_a100→128, gpu_h100→192 (default 0)
%                     .eur_per_1000_sbu      — EUR per 1000 SBUs (default 16, SURF 2025)
%                     .billing_profile       — label string for the report
%                   When parameters.hpc.name = 'snellius', Snellius rates are applied
%                   automatically unless cost_params overrides them.
%                   For other clusters (e.g. DCCN), cost fields will be NaN unless
%                   cost_params is supplied.
%
% Output:
%   report - struct with fields:
%     Scheduler info:
%       .job_id          — job ID string
%       .job_state       — final state (COMPLETED, FAILED, TIMEOUT, ...)
%       .exit_code       — exit code string
%     Actual usage:
%       .elapsed_s       — wall time in seconds
%       .cpu_core_s      — CPU core-seconds consumed
%       .alloc_cpus      — number of allocated CPU cores
%       .max_rss_gb      — peak memory in GB (from MaxRSS across all steps)
%       .n_gpus          — number of GPUs allocated (0 if none)
%     Requested resources:
%       .requested_mem_gb   — memory requested (GB)
%       .requested_timelimit — wall time limit string
%     Carbon metrics (Green Algorithms framework):
%       .energy_cpu_kwh     — energy used by CPUs
%       .energy_gpu_kwh     — energy used by GPUs
%       .energy_mem_kwh     — energy used by memory
%       .energy_total_kwh   — total energy including PUE overhead
%       .co2e_g             — CO2e in grams
%       .co2e_kg            — CO2e in kg
%     Monetary cost (cluster billing):
%       .cpu_sbu            — CPU SBUs consumed (core-hours × rate)
%       .gpu_sbu            — GPU SBUs consumed (GPU-hours × rate)
%       .total_sbu          — total SBUs consumed
%       .cost_eur           — estimated cost in EUR (NaN if rates unknown)
%       .billing_profile    — label identifying the billing profile used
%     Resource efficiency:
%       .memory_efficiency  — max_rss / requested_mem (NaN if MaxRSS unavailable)
%       .time_efficiency    — elapsed / requested wall time (NaN if unavailable)
%
% Notes:
%   - Only SLURM (sacct) is supported for automated metric collection.
%     For qsub, usage fields are set to NaN with a warning.
%   - MaxRSS is taken as the maximum across all job steps (batch, extern, etc.)
%     to capture peak memory reliably.
%   - Carbon defaults match the Donders/DCCN cluster (Netherlands grid, 2026).
%     Override via carbon_params for other sites.
%   - Snellius billing rates (SURF 2025): CPU 1 SBU/core-hour, A100 128 SBU/GPU-hour,
%     H100 192 SBU/GPU-hour, EUR 16 per 1000 SBU.
%     Source: https://servicedesk.surf.nl/wiki/spaces/WIKI/pages/30660209
%   - Cite: Lannelongue et al. 2021 https://doi.org/10.1002/advs.202100707
%
% See also: HPC_WAIT_FOR_COMPLETION, HPC_SUBMIT_JOB, HPC_JOB_INFO

% ---- Carbon constants (Green Algorithms framework defaults for Donders) ----
cp.pue                = 1.67;   % implied from cluster report: 40% overhead
cp.carbon_intensity   = 270;    % gCO2e/kWh, Netherlands electricity 2026
cp.tdp_per_core_w     = 12;     % W per CPU core (AMD EPYC / Intel Xeon class)
cp.tdp_per_gpu_w      = 300;    % W per GPU (A100 / V100 nominal)
cp.mem_power_w_per_gb = 0.3725; % W/GB RAM (Green Algorithms default)

if nargin >= 5 && ~isempty(carbon_params)
    fields = fieldnames(carbon_params);
    for i = 1:numel(fields)
        cp.(fields{i}) = carbon_params.(fields{i});
    end
end

% ---- Billing rates ----
% Snellius SURF 2025: https://servicedesk.surf.nl/wiki/spaces/WIKI/pages/30660209
% CPU: 1 SBU per core-hour. GPU: A100=128, H100=192 SBU per GPU-hour.
% EUR 16 per 1000 SBU (SURF rates 2025).
bp.cpu_sbu_per_core_hour = NaN;  % unknown unless profile matched below
bp.gpu_sbu_per_gpu_hour  = NaN;
bp.eur_per_1000_sbu      = NaN;
bp.billing_profile       = 'unknown';

% Auto-detect Snellius profile from parameters
if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'name') && ...
        strcmpi(parameters.hpc.name, 'snellius')
    bp.cpu_sbu_per_core_hour = 1;
    bp.eur_per_1000_sbu      = 16;
    if isfield(parameters.hpc, 'partition')
        switch lower(char(parameters.hpc.partition))
            case 'gpu_a100'
                bp.gpu_sbu_per_gpu_hour = 128;
                bp.billing_profile      = 'Snellius gpu_a100 (SURF 2025)';
            case 'gpu_h100'
                bp.gpu_sbu_per_gpu_hour = 192;
                bp.billing_profile      = 'Snellius gpu_h100 (SURF 2025)';
            otherwise
                bp.gpu_sbu_per_gpu_hour = 0;
                bp.billing_profile      = 'Snellius CPU (SURF 2025)';
        end
    else
        bp.gpu_sbu_per_gpu_hour = 0;
        bp.billing_profile      = 'Snellius (SURF 2025)';
    end
end

if nargin >= 6 && ~isempty(cost_params)
    fields = fieldnames(cost_params);
    for i = 1:numel(fields)
        bp.(fields{i}) = cost_params.(fields{i});
    end
end

% ---- Initialise report struct ----
job_id_str = sprintf('%.0f', job_id);

report.job_id           = job_id_str;
report.job_state        = 'UNKNOWN';
report.exit_code        = '';
report.elapsed_s        = NaN;
report.cpu_core_s       = NaN;
report.alloc_cpus       = NaN;
report.max_rss_gb       = NaN;
report.n_gpus           = 0;
report.requested_mem_gb = NaN;
report.requested_timelimit = '';
report.energy_cpu_kwh   = NaN;
report.energy_gpu_kwh   = NaN;
report.energy_mem_kwh   = NaN;
report.energy_total_kwh = NaN;
report.co2e_g           = NaN;
report.co2e_kg          = NaN;
report.memory_efficiency = NaN;
report.time_efficiency   = NaN;
report.cpu_sbu           = NaN;
report.gpu_sbu           = NaN;
report.total_sbu         = NaN;
report.cost_eur          = NaN;
report.billing_profile   = bp.billing_profile;

% ---- Requested resources from parameters ----
if isfield(parameters, 'hpc')
    if isfield(parameters.hpc, 'memorylimit')
        report.requested_mem_gb = parameters.hpc.memorylimit;
    end
    if isfield(parameters.hpc, 'timelimit')
        report.requested_timelimit = char(parameters.hpc.timelimit);
    end
end

% ---- Query scheduler ----
switch lower(hpc_type)
    case 'slurm'
        report = query_slurm(report, job_id_str, cp, bp);
    otherwise
        warning('PRESTUS:hpc_job_report:unsupportedScheduler', ...
            'hpc_job_report: post-job metrics not available for %s (only SLURM/sacct supported).', hpc_type);
end

% ---- Compute efficiency ratios ----
if ~isnan(report.max_rss_gb) && ~isnan(report.requested_mem_gb) && report.requested_mem_gb > 0
    report.memory_efficiency = report.max_rss_gb / report.requested_mem_gb;
end

if ~isnan(report.elapsed_s) && ~isempty(report.requested_timelimit)
    requested_s = parse_timelimit_to_seconds(report.requested_timelimit);
    if ~isnan(requested_s) && requested_s > 0
        report.time_efficiency = report.elapsed_s / requested_s;
    end
end

% ---- Print and write report ----
print_report(report, cp, bp);
if nargin >= 4 && ~isempty(log_dir) && isfolder(log_dir)
    write_report_file(report, parameters, log_dir, cp, bp);
end

end

% ==========================================================================
% LOCAL FUNCTIONS
% ==========================================================================

function report = query_slurm(report, job_id_str, cp, bp)
% Query main job allocation for timing and CPU info
cmd_alloc = sprintf( ...
    'sacct -j %s --allocations --format=State,ExitCode,ElapsedRaw,CPUTimeRAW,AllocCPUS,AllocTRES --noheader --parsable2 2>/dev/null', ...
    job_id_str);
[status, out] = system(cmd_alloc);

if status == 0 && ~isempty(strtrim(out))
    lines = strtrim(strsplit(out, '\n'));
    lines = lines(~cellfun(@isempty, lines));
    if ~isempty(lines)
        parts = strsplit(lines{1}, '|');
        if numel(parts) >= 5
            report.job_state  = strtrim(parts{1});
            report.exit_code  = strtrim(parts{2});
            report.elapsed_s  = str2double(strtrim(parts{3}));
            report.cpu_core_s = str2double(strtrim(parts{4}));
            report.alloc_cpus = str2double(strtrim(parts{5}));
        end
        if numel(parts) >= 6
            report.n_gpus = parse_gpu_count(strtrim(parts{6}));
        end
    end
end

% MaxRSS lives in job steps, not the allocation — take max across all steps
cmd_rss = sprintf( ...
    'sacct -j %s --format=MaxRSS --noheader --parsable2 2>/dev/null | grep -v "^$"', ...
    job_id_str);
[status_r, out_r] = system(cmd_rss);

if status_r == 0 && ~isempty(strtrim(out_r))
    rss_vals = strtrim(strsplit(out_r, '\n'));
    rss_vals = rss_vals(~cellfun(@isempty, rss_vals));
    rss_bytes = cellfun(@parse_rss_to_gb, rss_vals);
    rss_bytes = rss_bytes(~isnan(rss_bytes));
    if ~isempty(rss_bytes)
        report.max_rss_gb = max(rss_bytes);
    end
end

% ---- Carbon calculation ----
if ~isnan(report.elapsed_s) && ~isnan(report.alloc_cpus)
    elapsed_h = report.elapsed_s / 3600;
    n_cpus    = report.alloc_cpus;
    n_gpus    = report.n_gpus;
    mem_gb    = report.max_rss_gb;
    if isnan(mem_gb), mem_gb = 0; end

    report.energy_cpu_kwh = elapsed_h * n_cpus * cp.tdp_per_core_w / 1000;
    report.energy_gpu_kwh = elapsed_h * n_gpus * cp.tdp_per_gpu_w  / 1000;
    report.energy_mem_kwh = elapsed_h * mem_gb  * cp.mem_power_w_per_gb / 1000;

    energy_compute        = report.energy_cpu_kwh + report.energy_gpu_kwh + report.energy_mem_kwh;
    report.energy_total_kwh = energy_compute * cp.pue;
    report.co2e_g           = report.energy_total_kwh * cp.carbon_intensity;
    report.co2e_kg          = report.co2e_g / 1000;

    % ---- Monetary cost (SBU-based billing) ----
    if ~isnan(bp.cpu_sbu_per_core_hour)
        report.cpu_sbu = elapsed_h * n_cpus * bp.cpu_sbu_per_core_hour;
    end
    if ~isnan(bp.gpu_sbu_per_gpu_hour)
        report.gpu_sbu = elapsed_h * n_gpus * bp.gpu_sbu_per_gpu_hour;
    end
    cpu_sbu_val = report.cpu_sbu; if isnan(cpu_sbu_val), cpu_sbu_val = 0; end
    gpu_sbu_val = report.gpu_sbu; if isnan(gpu_sbu_val), gpu_sbu_val = 0; end
    if ~isnan(bp.eur_per_1000_sbu) && (~isnan(report.cpu_sbu) || ~isnan(report.gpu_sbu))
        report.total_sbu = cpu_sbu_val + gpu_sbu_val;
        report.cost_eur  = report.total_sbu * bp.eur_per_1000_sbu / 1000;
    end
end
end

% --------------------------------------------------------------------------
function n = parse_gpu_count(tres_str)
% Parse GPU count from AllocTRES string, e.g. "cpu=16,gres/gpu=1,mem=128G"
n = 0;
tok = regexp(tres_str, 'gres/gpu[^,=]*=(\d+)', 'tokens', 'once');
if ~isempty(tok)
    n = str2double(tok{1});
end
end

% --------------------------------------------------------------------------
function gb = parse_rss_to_gb(rss_str)
% Convert sacct MaxRSS string (e.g. "12345678K", "1234M", "1G") to GB
gb = NaN;
rss_str = strtrim(rss_str);
if isempty(rss_str) || strcmp(rss_str, '0'), return; end
tok = regexp(rss_str, '^([\d.]+)([KMGTP]?)$', 'tokens', 'once');
if isempty(tok), return; end
val  = str2double(tok{1});
unit = tok{2};
switch upper(unit)
    case 'K', gb = val / (1024^2);
    case 'M', gb = val / 1024;
    case 'G', gb = val;
    case 'T', gb = val * 1024;
    otherwise, gb = val / (1024^3); % assume bytes if no unit
end
end

% --------------------------------------------------------------------------
function s = parse_timelimit_to_seconds(tl_str)
% Convert SLURM time strings (D-HH:MM:SS, HH:MM:SS, MM:SS) to seconds
s = NaN;
tl_str = strtrim(char(tl_str));
% D-HH:MM:SS
tok = regexp(tl_str, '^(\d+)-(\d+):(\d+):(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    s = str2double(tok{1})*86400 + str2double(tok{2})*3600 + ...
        str2double(tok{3})*60   + str2double(tok{4});
    return;
end
% HH:MM:SS
tok = regexp(tl_str, '^(\d+):(\d+):(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    s = str2double(tok{1})*3600 + str2double(tok{2})*60 + str2double(tok{3});
    return;
end
% MM:SS
tok = regexp(tl_str, '^(\d+):(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    s = str2double(tok{1})*60 + str2double(tok{2});
end
end

% --------------------------------------------------------------------------
function print_report(r, cp, bp)
sep = repmat('─', 1, 50);
fprintf('\n%s\n', sep);
fprintf('  POST-JOB REPORT  (job %s)\n', r.job_id);
fprintf('%s\n', sep);

% Status
state_tag = r.job_state;
fprintf('  State:          %s  (exit: %s)\n', state_tag, r.exit_code);

% Timing
if ~isnan(r.elapsed_s)
    fprintf('  Wall time:      %s  (of %s requested)\n', ...
        seconds_to_hms(r.elapsed_s), r.requested_timelimit);
    if ~isnan(r.time_efficiency)
        fprintf('  Time used:      %.1f%%\n', r.time_efficiency * 100);
    end
end

% CPU / GPU
if ~isnan(r.alloc_cpus)
    fprintf('  CPUs allocated: %.0f\n', r.alloc_cpus);
end
if r.n_gpus > 0
    fprintf('  GPUs allocated: %.0f\n', r.n_gpus);
end

% Memory
if ~isnan(r.max_rss_gb)
    fprintf('  Peak memory:    %.2f GB  (of %.0f GB requested)\n', ...
        r.max_rss_gb, r.requested_mem_gb);
    if ~isnan(r.memory_efficiency)
        fprintf('  Memory used:    %.1f%%', r.memory_efficiency * 100);
        if r.memory_efficiency < 0.2
            fprintf('  ← consider reducing --mem to ~%.0f GB', ...
                ceil(r.max_rss_gb * 1.5));
        end
        fprintf('\n');
    end
end

% Carbon
fprintf('%s\n', sep);
fprintf('  CARBON ESTIMATE  [Green Algorithms, CI=%.0f gCO2/kWh, PUE=%.2f]\n', ...
    cp.carbon_intensity, cp.pue);
if ~isnan(r.co2e_kg)
    fprintf('  Energy (total): %.4f kWh\n', r.energy_total_kwh);
    fprintf('    CPUs:         %.4f kWh\n', r.energy_cpu_kwh);
    if r.n_gpus > 0
        fprintf('    GPUs:         %.4f kWh\n', r.energy_gpu_kwh);
    end
    if ~isnan(r.max_rss_gb) && r.max_rss_gb > 0
        fprintf('    Memory:       %.4f kWh\n', r.energy_mem_kwh);
    end
    fprintf('  CO2e:           %.4f kg  (%.1f g)\n', r.co2e_kg, r.co2e_g);
else
    fprintf('  (insufficient data for carbon estimate)\n');
end
fprintf('  Cite: Lannelongue et al. (2021) Adv Sci https://doi.org/10.1002/advs.202100707\n');

% Cost
fprintf('%s\n', sep);
if ~isnan(r.cost_eur)
    fprintf('  COST ESTIMATE  [%s]\n', r.billing_profile);
    if ~isnan(r.cpu_sbu)
        fprintf('    CPU SBUs:     %.1f\n', r.cpu_sbu);
    end
    if ~isnan(r.gpu_sbu) && r.gpu_sbu > 0
        fprintf('    GPU SBUs:     %.1f\n', r.gpu_sbu);
    end
    fprintf('    Total SBUs:   %.1f\n',  r.total_sbu);
    fprintf('    Cost:         EUR %.2f  (@ EUR %.0f / 1000 SBU)\n', ...
        r.cost_eur, bp.eur_per_1000_sbu);
else
    fprintf('  COST ESTIMATE  (no billing profile — pass cost_params to enable)\n');
end
fprintf('%s\n\n', sep);
end

% --------------------------------------------------------------------------
function write_report_file(report, parameters, log_dir, cp, bp)
subj_str = sprintf('sub-%03d', parameters.subject_id);
fname = fullfile(log_dir, sprintf('%s_job%s_report.txt', subj_str, report.job_id));

fid = fopen(fname, 'w');
if fid == -1
    warning('PRESTUS:hpc_job_report:writeError', ...
        'hpc_job_report: could not write report to %s', fname);
    return;
end

fprintf(fid, 'job_id: %s\n',         report.job_id);
fprintf(fid, 'job_state: %s\n',      report.job_state);
fprintf(fid, 'exit_code: %s\n',      report.exit_code);
fprintf(fid, 'subject_id: %s\n',     subj_str);
fprintf(fid, 'timestamp: %s\n',      datestr(now, 'yyyy-mm-dd HH:MM:SS'));

fprintf(fid, '\n# Actual usage\n');
fprintf(fid, 'elapsed_s: %.0f\n',    report.elapsed_s);
fprintf(fid, 'elapsed_hms: %s\n',    seconds_to_hms(report.elapsed_s));
fprintf(fid, 'cpu_core_s: %.0f\n',   report.cpu_core_s);
fprintf(fid, 'alloc_cpus: %.0f\n',   report.alloc_cpus);
fprintf(fid, 'n_gpus: %.0f\n',       report.n_gpus);
fprintf(fid, 'max_rss_gb: %.4f\n',   report.max_rss_gb);

fprintf(fid, '\n# Requested resources\n');
fprintf(fid, 'requested_mem_gb: %.0f\n',    report.requested_mem_gb);
fprintf(fid, 'requested_timelimit: %s\n',   report.requested_timelimit);

fprintf(fid, '\n# Resource efficiency\n');
fprintf(fid, 'memory_efficiency: %.4f\n',   report.memory_efficiency);
fprintf(fid, 'time_efficiency: %.4f\n',     report.time_efficiency);

fprintf(fid, '\n# Carbon estimate (Green Algorithms)\n');
fprintf(fid, '# cite: https://doi.org/10.1002/advs.202100707\n');
fprintf(fid, 'carbon_intensity_g_per_kwh: %.0f\n', cp.carbon_intensity);
fprintf(fid, 'pue: %.2f\n',                        cp.pue);
fprintf(fid, 'tdp_per_core_w: %.1f\n',             cp.tdp_per_core_w);
fprintf(fid, 'tdp_per_gpu_w: %.1f\n',              cp.tdp_per_gpu_w);
fprintf(fid, 'mem_power_w_per_gb: %.4f\n',         cp.mem_power_w_per_gb);
fprintf(fid, 'energy_cpu_kwh: %.6f\n',     report.energy_cpu_kwh);
fprintf(fid, 'energy_gpu_kwh: %.6f\n',     report.energy_gpu_kwh);
fprintf(fid, 'energy_mem_kwh: %.6f\n',     report.energy_mem_kwh);
fprintf(fid, 'energy_total_kwh: %.6f\n',   report.energy_total_kwh);
fprintf(fid, 'co2e_g: %.4f\n',             report.co2e_g);
fprintf(fid, 'co2e_kg: %.6f\n',            report.co2e_kg);

fprintf(fid, '\n# Cost estimate (cluster billing)\n');
fprintf(fid, 'billing_profile: %s\n',      report.billing_profile);
fprintf(fid, 'cpu_sbu_per_core_hour: %s\n', format_nan(bp.cpu_sbu_per_core_hour, '%.4f'));
fprintf(fid, 'gpu_sbu_per_gpu_hour: %s\n',  format_nan(bp.gpu_sbu_per_gpu_hour,  '%.4f'));
fprintf(fid, 'eur_per_1000_sbu: %s\n',      format_nan(bp.eur_per_1000_sbu,       '%.2f'));
fprintf(fid, 'cpu_sbu: %s\n',               format_nan(report.cpu_sbu,   '%.2f'));
fprintf(fid, 'gpu_sbu: %s\n',               format_nan(report.gpu_sbu,   '%.2f'));
fprintf(fid, 'total_sbu: %s\n',             format_nan(report.total_sbu, '%.2f'));
fprintf(fid, 'cost_eur: %s\n',              format_nan(report.cost_eur,  '%.4f'));

fclose(fid);
fprintf('  Report written: %s\n', fname);
end

% --------------------------------------------------------------------------
function hms = seconds_to_hms(s)
if isnan(s)
    hms = 'N/A';
    return;
end
s   = round(s);
h   = floor(s / 3600);
m   = floor(mod(s, 3600) / 60);
sec = mod(s, 60);
hms = sprintf('%02d:%02d:%02d', h, m, sec);
end

% --------------------------------------------------------------------------
function s = format_nan(val, fmt)
if isnan(val)
    s = 'NaN';
else
    s = sprintf(fmt, val);
end
end
