function out = run_stats(action, varargin)
% RUN_STATS  Collect a run's simulation size, duration and peak memory
%
% The acoustic and thermal stages report what they simulated; at the end of
% a successful run the pipeline appends one JSON line to
% <path.sim>/run_stats.jsonl. The Nexus resource estimate (and its headless
% form, python -m prestus_nexus.estimate) reads that file. Nothing is sent
% anywhere; this is independent of telemetry consent.
%
% Use as:
%   run_stats('reset')
%   run_stats('add', 'acoustic_s', 812.4)        % sums over repeated calls
%   run_stats('set', 'acoustic_nt', 4200)
%   rec = run_stats('write', parameters, total_s) % append; returns the record
%   rec = run_stats('get')
%
% Input:
%   action     - 'reset' | 'set' | 'add' | 'get' | 'write'
%   varargin   - name/value pairs ('set', 'add'); parameters and total run
%                time [s] ('write')
%
% Output:
%   out - the record struct ('get', 'write'); empty otherwise
%
% Record fields: timestamp, platform, hpc_name, hpc_partition, host,
% n_cpu, gpu_model, gpu_mem_gb, code_type, precision, acoustic_backend, medium,
% resolution_mm, acoustic_voxels, acoustic_dims, acoustic_pml, axisymmetric,
% acoustic_nt, acoustic_s, thermal_voxels, thermal_steps, thermal_s, total_s,
% peak_ram_gb, and stage_<label>_s per timed stage (LOG_TIMER labels, e.g.
% stage_preproc_s, stage_freefield_baseline_s; summed over repeats), and the
% modelled skull properties (SKULL_PROPERTY_STATS: skull_c/_rho/_atten
% _p10/_median/_p90, skull_path_mm, skull_path_*_median), and queue_wait_s
% (seconds the SLURM job waited; SLURM only). The acoustic fields describe the grid the solver ran on (the
% 2D plane of an axisymmetric run; PML not included in voxels or dims).
% n_cpu is the job's CPU allocation under SLURM, else the machine's cores.
% A stage that did not simulate (cache reuse, dry run, stage off) leaves its
% fields out.
%
% Peak RAM is the MATLAB process high-water mark, or the job cgroup peak
% (which includes child solvers) under SLURM; Linux only, NaN elsewhere.
%
% See also: LOG_MEMORY_USAGE, PRESTUS_PIPELINE

persistent S
if isempty(S), S = struct(); end
out = [];

switch action
    case 'reset'
        S = struct();
    case 'set'
        for k = 1:2:numel(varargin)
            S.(varargin{k}) = varargin{k+1};
        end
    case 'add'
        for k = 1:2:numel(varargin)
            if isfield(S, varargin{k})
                S.(varargin{k}) = S.(varargin{k}) + varargin{k+1};
            else
                S.(varargin{k}) = varargin{k+1};
            end
        end
    case 'get'
        out = S;
    case 'write'
        out = S;
        S = struct();   % a nested or follow-up run starts clean
        % Only runs that simulated something are useful for estimates.
        if ~isfield(out, 'acoustic_s') && ~isfield(out, 'thermal_s')
            return
        end
        try
            out = build_record(out, varargin{1}, varargin{2});
            append_record(out, varargin{1});
        catch me
            % Bookkeeping must never fail a finished run.
            fprintf('[run_stats] not recorded: %s\n', me.message);
        end
    otherwise
        error('run_stats:badAction', 'Unknown action ''%s''.', action);
end
end

%% ------------------------------------------------------------------------
function rec = build_record(S, parameters, total_s)
rec = struct();
rec.timestamp        = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
rec.platform         = char(get_nonempty_field(parameters, {'platform'}, 'matlab'));
rec.hpc_name         = char(get_nonempty_field(parameters, {'hpc', 'name'}, ''));
rec.hpc_partition    = char(get_nonempty_field(parameters, {'hpc', 'partition'}, ''));
rec.host             = host_name();
rec.n_cpu            = cpu_count();
rec.code_type        = char(get_nonempty_field(parameters, {'simulation', 'code_type'}, ''));
rec.gpu_model        = '';
if contains(rec.code_type, 'gpu')
    [rec.gpu_model, rec.gpu_mem_gb] = gpu_name();
end
rec.precision        = char(get_nonempty_field(parameters, {'simulation', 'precision'}, ''));
rec.acoustic_backend = char(get_nonempty_field(parameters, {'acoustic_backend'}, 'kwave'));
rec.medium           = char(get_nonempty_field(parameters, {'simulation', 'medium'}, ''));
rec.resolution_mm    = get_nonempty_field(parameters, {'grid', 'resolution_mm'}, NaN);
fn = fieldnames(S);
for k = 1:numel(fn)
    rec.(fn{k}) = S.(fn{k});
end
rec.total_s     = total_s;
rec.peak_ram_gb = peak_ram_gb();
qw = queue_wait_s();
if isfinite(qw)
    rec.queue_wait_s = qw;
end
end

function append_record(rec, parameters)
sim_dir = char(get_nonempty_field(parameters, {'path', 'sim'}, ''));
if isempty(sim_dir) || ~isfolder(sim_dir)
    return
end
fid = fopen(fullfile(sim_dir, 'run_stats.jsonl'), 'a');
if fid < 0
    return
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', jsonencode(rec));
end

function name = host_name()
name = getenv('SLURMD_NODENAME');
if isempty(name)
    [~, name] = system('hostname');
end
name = strtrim(name);
end

function n = cpu_count()
n = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isfinite(n)
    n = feature('numcores');
end
end

function [name, mem_gb] = gpu_name()
name = '';
mem_gb = NaN;
if canUseGPU()   % false without a GPU or the Parallel Computing Toolbox
    g = gpuDevice();
    name = g.Name;
    mem_gb = double(g.TotalMemory) / 1e9;   % a MIG slice reports its own memory
end
end

function s = queue_wait_s()
% Seconds the SLURM job waited in the queue (StartTime - SubmitTime of
% `scontrol show job`); NaN outside SLURM or when scontrol is unavailable.
s = NaN;
id = getenv('SLURM_JOB_ID');
if isempty(id)
    return
end
try
    [status, out] = run_command({'scontrol', 'show', 'job', id, '-o'});
    if status ~= 0
        return
    end
    t = regexp(out, 'SubmitTime=(\S+).*?StartTime=(\S+)', 'tokens', 'once');
    if isempty(t)
        return
    end
    fmt = 'yyyy-MM-dd''T''HH:mm:ss';
    s = seconds(datetime(t{2}, 'InputFormat', fmt) - datetime(t{1}, 'InputFormat', fmt));
catch
    % A job not yet started reports StartTime=Unknown; the wait is then unknown.
    s = NaN;
end
end

function gb = peak_ram_gb()
gb = NaN;
if ~isunix || ismac || ~isfile('/proc/self/status')
    return
end
tok = regexp(fileread('/proc/self/status'), 'VmHWM:\s+(\d+)\s*kB', 'tokens', 'once');
if ~isempty(tok)
    gb = str2double(tok{1}) / 1024^2;
end
if ~isempty(getenv('SLURM_JOB_ID'))
    % The job cgroup also holds child solvers (k-Wave C++ binary, BabelBrain).
    % (evalc: log_memory_usage prints a line; keep the run log clean. It
    % returns NaN where cgroups are not readable.)
    [~, mem] = evalc('log_memory_usage(''run_stats'')');
    gb = max(gb, mem.cg_peak_gb);
end
end
