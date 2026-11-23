function hpc_scratch_finalize(scratch_dir, created_root)
% HPC_SCRATCH_FINALIZE  Remove a job's node-local scratch folder
%
% Called through the onCleanup object of HPC_SCRATCH_SETUP, so it runs when
% the pipeline returns or errors. Removes the folder, including .h5 files a
% failed k-Wave run left behind. Problems are reported as warnings, so they
% never hide the pipeline's own error.
%
% Use as:
%   hpc_scratch_finalize(scratch_dir)
%   hpc_scratch_finalize(scratch_dir, created_root)
%
% Input:
%   scratch_dir  - job-specific scratch folder (io.dir_scratch)
%   created_root - per-job scratch root that HPC_SCRATCH_SETUP created
%                  (e.g. /scratch/$USER/$SLURM_JOB_NAME); removed only if it
%                  is then empty. '' or omitted: the root is left alone.
%
% See also: HPC_SCRATCH_SETUP

arguments
    scratch_dir {mustBeTextScalar}
    created_root = ''
end

if nargin < 2
    created_root = '';
end
try
    if isfolder(scratch_dir)
        rmdir(scratch_dir, 's');
        fprintf('Removed node-local scratch: %s\n', scratch_dir);
    end
    [~, ~] = rmdir(fileparts(scratch_dir));   % prestus_<jobid>, only if now empty
    if ~isempty(created_root)
        [~, ~] = rmdir(created_root);         % per-job root created by setup, only if now empty
    end
catch me
    warning('prestus:scratch', 'Could not remove scratch folder %s: %s', scratch_dir, me.message);
end
end
