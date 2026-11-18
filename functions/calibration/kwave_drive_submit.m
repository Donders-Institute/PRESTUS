function job_id = kwave_drive_submit(script, fit_mat, sim_config)
% KWAVE_DRIVE_SUBMIT  Submit a k-Wave drive check script as a scheduler job
%
% The free-water k-Wave drive check (a script that calls KWAVE_CALIBRATE_DRIVE
% and saves the corrected entry, written by the PRESTUS Nexus calibration page)
% submits one simulation job per refined depth and waits for each. Run as a
% job itself, that waiting happens on a compute node, not on a login node or
% the user's machine. The job needs no GPU (the simulations request their own),
% 8 GB of memory, and a wall time of the per-simulation limit (hpc.timelimit of
% sim_config, else config_default) times the number of refined depths, plus one
% hour for queueing and saving. Prints 'Job "<name>" (ID: <id>) submitted
% successfully', which the Nexus Run tab follows.
%
% Use as:
%   job_id = kwave_drive_submit(script, fit_mat)
%   job_id = kwave_drive_submit(script, fit_mat, sim_config)
%
% Input:
%   script     - path of the drive-check .m script; the job runs it from its folder
%   fit_mat    - .mat with the unsaved fit (variable fit, fit.entry.refinement.depths_ep_mm)
%   sim_config - (optional) study config YAML the check simulates with ('' = none):
%                its hpc block (partition and limits of the simulations) and platform
%
% Output:
%   job_id     - scheduler job ID
%
% See also: KWAVE_CALIBRATE_DRIVE, HPC_SUBMIT_JOB, HPC_DETECT_SYSTEM

arguments
    script     (1,:) char
    fit_mat    (1,:) char
    sim_config (1,:) char = ''
end

[folder, name, ext] = fileparts(script);
if ~strcmp(ext, '.m') || ~isfile(script)
    error('kwave_drive_submit:noScript', 'Drive-check script not found: %s', script);
end
s = load(fit_mat, 'fit');
n_depths = numel(s.fit.entry.refinement.depths_ep_mm);

dflt = yaml.loadFile(fullfile(get_prestus_path(), 'config', 'config_default.yaml'), 'ConvertToArray', true);
hpc = dflt.hpc;
platform = 'auto';
if ~isempty(sim_config)
    c = yaml.loadFile(sim_config, 'ConvertToArray', true);
    if isfield(c, 'hpc') && isstruct(c.hpc)
        for f = fieldnames(c.hpc)'
            hpc.(f{1}) = c.hpc.(f{1});
        end
    end
    if isfield(c, 'platform') && ~isempty(c.platform)
        platform = char(c.platform);
    end
end
if strcmp(platform, 'matlab')
    error('kwave_drive_submit:localPlatform', ...
        ['platform is ''matlab'' in %s, so the check would simulate in the job itself without a GPU. ' ...
         'Set platform to ''auto'', ''slurm'' or ''qsub''.'], sim_config);
end
hpc_type = platform;
if ~ismember(hpc_type, {'slurm', 'qsub'})
    hpc_type = hpc_detect_system();
end
if ~ismember(hpc_type, {'slurm', 'qsub'})
    error('kwave_drive_submit:noScheduler', 'No SLURM or qsub scheduler on this machine.');
end
if isfield(hpc, 'name') && strcmp(char(hpc.name), 'snellius')
    error('kwave_drive_submit:snellius', ...
        ['Snellius jobs need a GPU partition; this waiting job runs none. Submit %s as a CPU job ' ...
         'yourself (matlab -batch %s in %s).'], script, name, folder);
end

% The waiting job: no GPU, little memory, the simulations' time one after another
per_sim = walltime_s(hpc.timelimit);
total = n_depths * per_sim + 3600;
p = struct();
p.subject_id = 0;
p.hpc = hpc;
p.hpc.gpu = '';
p.hpc.partition = '';
p.hpc.memorylimit = 8;
p.hpc.timelimit = sprintf('%d:%02d:%02d', floor(total / 3600), floor(mod(total, 3600) / 60), mod(total, 60));
p.hpc.job_name = 'PS_kwave_drive';
p.acoustic_backend = 'kwave';
p.simulation.code_type = 'matlab_cpu';       % the job itself runs no solver: no GPU request
p.modules.run_heating_sims = 0;

fprintf('[kwave_drive_submit] %d refined depths x %s per simulation + 1 h: wall time %s, 8 GB, no GPU.\n', ...
    n_depths, char(hpc.timelimit), p.hpc.timelimit);
job_id = hpc_submit_job(hpc_type, name, p, folder);   % prints 'Job "..." (ID: ...) submitted successfully'
end

function s = walltime_s(t)
% '[D-]HH:MM:SS' (also MM:SS or minutes) -> seconds
t = strtrim(char(string(t)));
days = 0;
tok = regexp(t, '^(\d+)-(.*)$', 'tokens', 'once');
if ~isempty(tok)
    days = str2double(tok{1});
    t = tok{2};
end
parts = str2double(strsplit(t, ':'));
switch numel(parts)
    case 1, s = parts * 60;
    case 2, s = parts(1) * 60 + parts(2);
    otherwise, s = parts(1) * 3600 + parts(2) * 60 + parts(3);
end
s = s + days * 86400;
end
