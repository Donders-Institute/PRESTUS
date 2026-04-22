function coords_sub = mni2subject_coords_LDfix(coords_mni, subdir, parameters, transformation_type)
% MNI2SUBJECT_COORDS_LDFIX  Transform MNI coordinates to subject space via SimNIBS CLI (LD_LIBRARY_PATH hotfix)
%
% Writes coords_mni to a temporary CSV in SimNIBS Generic format, calls the
% mni2subject_coords command-line tool, and reads back the result. An optional
% LD_LIBRARY_PATH export is prepended to the system call when
% parameters.hpc.ld_library_path is set — this hotfix is required on some HPC
% systems where SimNIBS bash scripts do not inherit the library path correctly
% (see: https://github.com/simnibs/simnibs/issues/106). Has a large per-call
% overhead due to the system() call; do not call inside tight loops.
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
%   transformation_type - transform type: 'nonl' | '12dof' | '6dof'
%                         (optional, default: 'nonl')
%
% Output:
%   coords_sub - [Nx3] coordinates in subject RAS+ space [mm]
%
% See also: SUBJECT2MNI_COORDS_LDFIX, TRANSFORM_COORDINATES, CONVERT_FINAL_TO_MNI_SIMNIBS

if nargin < 4
    transformation_type = 'nonl';
end
if ~any(strcmp(transformation_type, {'nonl', '12dof', '6dof'}))
    error('transformation_type must be nonl, 12dof, or 6dof')
end
assert(size(coords_mni, 2) == 3, 'coords_subject must be in Nx3 format');
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

% include LD fix (Unix only; Windows does not use LD_LIBRARY_PATH)
if isfield(parameters.hpc,'ld_library_path') && ~isempty(parameters.hpc.ld_library_path) && ~ispc
    ld_command = sprintf('export LD_LIBRARY_PATH="%s"; ', parameters.hpc.ld_library_path);
else
    ld_command = '';
end

% Run mni2subject_coords (quote paths for spaces)
mni2sub_bin = fullfile(parameters.startup.simnibs_bin_path, 'mni2subject_coords');
[status,result] = system(sprintf('%s"%s" -m "%s" -s "%s" -o "%s" -t %s', ...
    ld_command, mni2sub_bin, subdir, fn_in, fn_out, transformation_type));


% system([simnibs_cli_call('mni2subject_coords')...
%                          ' -m "' subdir '" -s ' fn_in ' -o ' fn_out ...
%                          ' -t ' transformation_type]);

% Check if call was successefull
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running mni2subject_coords:\n %s',result)
end

% Read output
coords_sub = csvread(fn_out, 0, 1);
delete(fn_in);
delete(fn_out);

end
