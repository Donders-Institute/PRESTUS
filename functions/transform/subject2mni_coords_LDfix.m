function coords_mni = subject2mni_coords_LDfix(coords_subject, subdir, parameters, transformation_type)
% SUBJECT2MNI_COORDS_LDFIX  Transform subject RAS+ coordinates to MNI space via SimNIBS CLI (LD_LIBRARY_PATH hotfix)
%
% Writes coords_subject to a temporary CSV in SimNIBS Generic format, calls the
% subject2mni_coords command-line tool, and reads back the result. An optional
% LD_LIBRARY_PATH is prefixed on the child process's environment with
% hpc.ld_library_path, else the install's libraries (SIMNIBS_ENV_OVERRIDES) — this hotfix is required on some HPC
% systems where SimNIBS bash scripts do not inherit the library path correctly
% (see: https://github.com/simnibs/simnibs/issues/106). Has a large per-call
% overhead due to the external process call; do not call inside tight loops.
%
% Use as:
%   coords_mni = subject2mni_coords_LDfix(coords_subject, subdir, parameters)
%   coords_mni = subject2mni_coords_LDfix(coords_subject, subdir, parameters, transformation_type)
%
% Input:
%   coords_subject       - [Nx3] coordinates in subject RAS+ space [mm]
%   subdir               - path to SimNIBS m2m subject directory (e.g. 'm2m_ernie/')
%   parameters           - PRESTUS config; uses hpc.ld_library_path and
%                         startup.simnibs_bin_path
%   transformation_type - transform type: 'nonl' | '12dof' | '6dof'
%                         (optional, default: 'nonl')
%
% Output:
%   coords_mni - [Nx3] coordinates in MNI space [mm]
%
% See also: MNI2SUBJECT_COORDS_LDFIX, TRANSFORM_COORDINATES, CONVERT_FINAL_TO_MNI_SIMNIBS

if nargin < 4
    transformation_type = 'nonl';
end
if ~any(strcmp(transformation_type, {'nonl', '12dof', '6dof'}))
    error('transformation_type must be nonl, 12dof, or 6dof')
end
assert(size(coords_subject, 2) == 3, 'coords_subject must be in Nx3 format');
assert(exist(subdir, 'dir') == 7, ['Could not find directory ' subdir])

% Write CSV file in SimNIBS format
fn_in = [tempname,'.csv'];
fid = fopen( fn_in, 'w' );
for i = 1 : size(coords_subject, 1)
    fprintf(...
        fid, 'Generic, %f, %f, %f\n', ...
        coords_subject(i,1), coords_subject(i, 2),  coords_subject(i, 3));
end
fclose( fid );
fn_out = [tempname,'.csv'];

% include LD fix. run_command's Env option sets LD_LIBRARY_PATH natively
% on the child process -- no `export ...;` shell syntax needed at all.
% (hpc.ld_library_path, else the selected install's libraries: SIMNIBS_ENV_OVERRIDES)
env_overrides = simnibs_env_overrides(parameters);

% Run subject2mni_coords
simnibs_bin = parameters.startup.simnibs_bin_path;
assert(~isempty(simnibs_bin) && ~strcmp(simnibs_bin, '/'), ...
    'parameters.startup.simnibs_bin_path must be set to the SimNIBS bin directory (not empty or ''/'')')
subject2mni_bin = fullfile(simnibs_bin, 'subject2mni_coords');
% run_command invokes subject2mni_coords directly via argv, no shell at
% all -- subdir (the m2m folder, which may come from a shared/untrusted
% data directory) reaches the process as a literal argument regardless
% of its contents.
[status, result] = run_command({subject2mni_bin, '-m', subdir, '-s', fn_in, '-o', fn_out, '-t', transformation_type}, ...
    'Env', env_overrides);

% Check if call was successful
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running subject2mni_coords:\n %s',result)
end

% Read output
coords_mni = csvread(fn_out, 0, 1);
% SimNIBS returns [0 0 0] for a point it cannot map (outside the warp field):
% NaN, so callers cannot silently use the coordinate origin
unmapped = all(coords_mni == 0, 2);
if any(unmapped)
    coords_mni(unmapped, :) = NaN;
    warning('PRESTUS:coords:noMapping', 'subject -> MNI: %d of %d points could not be mapped (SimNIBS returned 0 0 0); returned as NaN.', ...
        nnz(unmapped), numel(unmapped));
end
delete(fn_in);
delete(fn_out);

end
