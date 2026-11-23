function coords_sub = mni2subject_coords_LDfix(coords_mni, subdir, parameters, transformation_type)
% MNI2SUBJECT_COORDS_LDFIX  Transform MNI coordinates to subject space via SimNIBS CLI (LD_LIBRARY_PATH hotfix)
%
% Writes coords_mni to a temporary CSV in SimNIBS Generic format, calls the
% mni2subject_coords command-line tool, and reads back the result. An optional
% LD_LIBRARY_PATH is prefixed on the child process's environment with
% hpc.ld_library_path, else the install's libraries (SIMNIBS_ENV_OVERRIDES) — this hotfix is required on some HPC
% systems where SimNIBS bash scripts do not inherit the library path correctly
% (see: https://github.com/simnibs/simnibs/issues/106). Has a large per-call
% overhead due to the external process call; do not call inside tight loops.
%
% Use as:
%   coords_sub = mni2subject_coords_LDfix(coords_mni, subdir, parameters)
%   coords_sub = mni2subject_coords_LDfix(coords_mni, subdir, parameters, transformation_type)
%
% Input:
%   coords_mni          - [Nx3] coordinates in MNI space [mm]
%   subdir              - path to SimNIBS m2m subject directory (e.g. 'm2m_ernie/')
%   parameters          - PRESTUS config; uses hpc.ld_library_path and
%                         startup.simnibs_bin_path
%   transformation_type - 'nonl' (charm's nonlinear warp; the only transform
%                         SimNIBS 4 writes)
%                         (optional, default: 'nonl')
%
% Output:
%   coords_sub - [Nx3] coordinates in subject RAS+ space [mm]
%
% See also: SUBJECT2MNI_COORDS_LDFIX, TRANSFORM_COORDINATES, CONVERT_FINAL_TO_MNI_SIMNIBS

if nargin < 4
    transformation_type = 'nonl';
end
if ~strcmp(transformation_type, 'nonl')
    error('prestus:transform:simnibs4Only', ...
        ['Transform ''%s'' needs the affine files of a SimNIBS 3 (headreco) segmentation. ' ...
         'PRESTUS uses charm segmentations (SimNIBS 4 and later), which hold only the ' ...
         'nonlinear warp: use ''nonl''.'], transformation_type);
end
assert(size(coords_mni, 2) == 3, 'coords_mni must be in Nx3 format');
assert(exist(subdir, 'dir') == 7, ['Could not find directory ' subdir])

% Write CSV file in SimNIBS format
fn_in = [tempname,'.csv'];
fid = fopen( fn_in, 'w' );
for i = 1 : size(coords_mni, 1)
    fprintf(...
        fid, 'Generic, %f, %f, %f\n', ...
        coords_mni(i,1), coords_mni(i, 2),  coords_mni(i, 3));
end
fclose( fid );
fn_out = [tempname,'.csv'];

% include LD fix (Unix only; Windows does not use LD_LIBRARY_PATH).
% run_command's Env option sets LD_LIBRARY_PATH natively on the child
% process -- no `export ...;` shell syntax needed at all.
% (hpc.ld_library_path, else the selected install's libraries: SIMNIBS_ENV_OVERRIDES)
env_overrides = simnibs_env_overrides(parameters);

% run_command invokes mni2subject_coords directly via argv, no shell at
% all -- subdir (the m2m folder, which may come from a shared/untrusted
% data directory) reaches the process as a literal argument regardless
% of its contents.
mni2sub_bin = fullfile(parameters.startup.simnibs_bin_path, 'mni2subject_coords');
[status, result] = run_command({mni2sub_bin, '-m', subdir, '-s', fn_in, '-o', fn_out, '-t', transformation_type}, ...
    'Env', env_overrides);

% Check if call was successful
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running mni2subject_coords:\n %s',result)
end

% Read output
coords_sub = csvread(fn_out, 0, 1);
% SimNIBS returns [0 0 0] for a point it cannot map (outside the warp field):
% NaN, so callers cannot silently use the coordinate origin
unmapped = all(coords_sub == 0, 2);
if any(unmapped)
    coords_sub(unmapped, :) = NaN;
    warning('PRESTUS:coords:noMapping', 'MNI -> subject: %d of %d points could not be mapped (SimNIBS returned 0 0 0); returned as NaN.', ...
        nnz(unmapped), numel(unmapped));
end
delete(fn_in);
delete(fn_out);

end
