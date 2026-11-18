function files = export_kplan_plan(t1_file, trans_pos_ras, focus_pos_ras, out_dir, name, opts)
% EXPORT_KPLAN_PLAN  Write a transducer pose and its images for k-Plan
%
% k-Plan counts world coordinates from the first voxel of its planning image
% along the image axes (world = 0-based voxel x voxel size), checked on plans
% 50-770 and 69-928. For an oblique T1, k-Plan would resample the image onto
% its world axes instead. So the images are written with the same voxels and
% the rotation removed from the header (voxel 0 at world 0), and the pose is
% given in that frame. Plans made from them import back exactly with the
% Nexus (Import trajectory). Distances and angles are unchanged.
%
% Files in out_dir:
%   <name>_T1w.nii          T1, rotation removed from the header
%   <name>_ct.nii           CT or pseudo-CT on the T1 grid (opts.ct_file), likewise
%   <name>_transducer.kps   transducer pose: bowl apex, beam along -z (k-Plan TPM v1.2.0 layout)
%   <name>_kplan_entry.txt  the values to check or type into k-Plan
%
% Use as:
%   files = export_kplan_plan(t1_file, trans_pos_ras, focus_pos_ras, out_dir, name)
%   files = export_kplan_plan(..., ct_file=ct, kplan_offset_mm=10.82, transducer_name='...')
%
% Input:
%   t1_file       - T1 NIfTI with an sform, voxel axes (obliquely) R, A, S
%   trans_pos_ras - [1x3] bowl apex in subject RAS [mm]
%   focus_pos_ras - [1x3] target in subject RAS [mm]
%   out_dir       - output folder (created)
%   name          - file name stem, e.g. 'sub-001'
%   opts.ct_file          - CT or pseudo-CT on the T1 grid ('' = none); k-Plan's planning image
%   opts.kplan_offset_mm  - k-Plan's apex -> exit plane distance for the device, if known;
%                           the entry sheet then gives k-Plan's focal distance
%   opts.transducer_name  - k-Plan transducer name for the entry sheet
%
% Output:
%   files - struct with fields t1, ct ('' if none), kps, entry
%
% See also: EXPORT_BRAINSIGHT_TRAJECTORY, EXPORT_SLICER_TRANSFORM, EXPORT_PLACEMENT

arguments
    t1_file       (1,:) char
    trans_pos_ras (1,3) double
    focus_pos_ras (1,3) double
    out_dir       (1,:) char
    name          (1,:) char
    opts.ct_file         (1,:) char = ''
    opts.kplan_offset_mm double = []
    opts.transducer_name (1,:) char = ''
end

F = local_frame(t1_file);
if ~isfolder(out_dir), mkdir(out_dir); end
files = struct('t1', fullfile(out_dir, [name '_T1w.nii']), 'ct', '', ...
    'kps', fullfile(out_dir, [name '_transducer.kps']), 'entry', fullfile(out_dir, [name '_kplan_entry.txt']));

local_write_nii(niftiread(t1_file), niftiinfo(t1_file), F, files.t1);
if ~isempty(opts.ct_file)
    Fc = local_frame(opts.ct_file);
    if isequal(Fc.dims, F.dims) && max(abs(Fc.affine(:) - F.affine(:))) < 1e-3
        files.ct = fullfile(out_dir, [name '_ct.nii']);
        local_write_nii(single(niftiread(opts.ct_file)), niftiinfo(opts.ct_file), F, files.ct);
    else
        warning('PRESTUS:kplan:ctGrid', 'CT %s is not on the T1 grid; not exported.', opts.ct_file);
    end
end

% Pose in k-Plan world [mm]
apex_w   = F.ras2world(trans_pos_ras);
target_w = F.ras2world(focus_pos_ras);
d_bowl = norm(target_w - apex_w);
if d_bowl < eps
    error('PRESTUS:kplan:degenerate', 'trans_pos_ras and focus_pos_ras are identical: no beam axis.');
end
beam = (target_w - apex_w) / d_bowl;              % apex -> target
z_t  = -beam;                                     % k-Plan: beam along -z of the transducer
a = [0 0 1]; if abs(dot(a, z_t)) > 0.9, a = [0 1 0]; end
x_t = cross(a, z_t); x_t = x_t / norm(x_t);
y_t = cross(z_t, x_t);
T_kplan = [x_t 0; y_t 0; z_t 0; apex_w / 1000 1];  % row-vector form, metres
local_write_kps(files.kps, T_kplan, name);

% Entry sheet
fid = fopen(files.entry, 'w');
if fid == -1, error('PRESTUS:kplan:ioError', 'Cannot write %s', files.entry); end
cleanup = onCleanup(@() fclose(fid));   %#ok<NASGU>
fprintf(fid, 'k-Plan plan for %s (PRESTUS export_kplan_plan, %s)\n\n', name, char(datetime('now')));
if ~isempty(files.ct)
    fprintf(fid, 'Planning image (primary)  %s_ct.nii   as CT\n', name);
    fprintf(fid, 'Secondary image           %s_T1w.nii  registration: identity (same voxels and header)\n', name);
else
    fprintf(fid, 'Images                    %s_T1w.nii  (no CT on the T1 grid was exported; add the planning CT in k-Plan)\n', name);
end
if ~isempty(opts.transducer_name)
    fprintf(fid, 'Transducer                %s\n', opts.transducer_name);
end
fprintf(fid, 'Transducer position       import %s_transducer.kps (label %s)\n', name, name);
fprintf(fid, '  or enter by hand:\n');
fprintf(fid, '  bowl apex  [mm, world]  %9.3f %9.3f %9.3f\n', apex_w);
fprintf(fid, '  beam dir   [world]      %9.5f %9.5f %9.5f   (apex -> target; transducer -z)\n', beam);
fprintf(fid, 'Target                    [mm, world]  %9.3f %9.3f %9.3f\n', target_w);
fprintf(fid, 'Apex to target            %.2f mm\n', d_bowl);
if ~isempty(opts.kplan_offset_mm)
    fprintf(fid, 'Focal distance            %.2f mm from k-Plan''s exit plane (%.2f mm in front of its bowl apex)\n', ...
        d_bowl - opts.kplan_offset_mm, opts.kplan_offset_mm);
else
    fprintf(fid, 'Focal distance            apex to target minus k-Plan''s apex-to-exit-plane distance for this transducer\n');
end
fprintf(fid, ['\nWorld coordinates count from the first voxel of these images (world = voxel index x voxel size).\n' ...
    'Upload these images, not the originals: for an oblique T1, k-Plan would otherwise use other world coordinates.\n']);
end

% ── Frame ────────────────────────────────────────────────────────────────

function F = local_frame(nifti_file)
info = niftiinfo(nifti_file);
if info.raw.sform_code > 0
    A = [info.raw.srow_x; info.raw.srow_y; info.raw.srow_z; 0 0 0 1];
else
    A = info.Transform.T';
end
zooms = sqrt(sum(A(1:3,1:3).^2, 1));
R = A(1:3,1:3) ./ zooms;
[U, ~, V] = svd(R);
if any(diag(U * V') < cosd(45))
    error('PRESTUS:kplan:notRAS', '%s: voxel axes are not (obliquely) R, A, S; reorient it first.', nifti_file);
end
F.dims  = info.ImageSize(1:3);
F.zooms = zooms;
F.affine = A;
F.affine_kplan = [diag(zooms) zeros(3,1); 0 0 0 1];
F.ras2world = @(r) ((A \ [r(:); 1])' * [eye(3); 0 0 0]) .* zooms;   % 0-based voxel x voxel size
end

function local_write_nii(img, info, F, filename)
% same voxels; header without rotation: voxel 0 at world 0, axes = world axes
info.Transform = affine3d(F.affine_kplan');
info.TransformName = 'Sform';
info.Qfactor = 1;
info.Description = 'PRESTUS k-Plan export: rotation removed';
info.Datatype = class(img);
info.BitsPerPixel = 8 * numel(typecast(cast(0, class(img)), 'uint8'));
info.MultiplicativeScaling = 1;
info.AdditiveOffset = 0;
[d, n] = fileparts(filename);
niftiwrite(img, fullfile(d, n), info, 'Compressed', false);
end

% ── k-Plan transducer position file (.kps) ───────────────────────────────

function local_write_kps(filename, T_kplan, label)
% HDF5 as written by k-Plan TPM v1.2.0: root attributes (fixed-length ASCII,
% uint64 versions), group '/1' with transform_label and position_transform
% (single, 1 x 4 x 4 in HDF5 order). Accepted by k-Plan (plan 69-928).
if isfile(filename), delete(filename); end
stamp = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm''Z'''));
h5create(filename, '/1/position_transform', [4 4 1], 'Datatype', 'single', ...
    'ChunkSize', [4 4 1], 'Deflate', 4);
h5write(filename, '/1/position_transform', single(reshape(T_kplan', [4 4 1])));
local_str(filename, '/', 'application_name', 'k-Plan');
local_str(filename, '/', 'created_by',       'PRESTUS');
local_str(filename, '/', 'creation_date',    stamp);
local_str(filename, '/', 'file_description', 'k-Plan Transducer Position');
local_str(filename, '/', 'file_type',        'k-Plan Transducer Position');
local_str(filename, '/', 'modified_date',    stamp);
local_u64(filename, '/', 'major_version', 1);
local_u64(filename, '/', 'minor_version', 3);
local_u64(filename, '/', 'number_transforms', 1);
local_str(filename, '/1', 'transform_label', label);
end

function local_str(file, loc, name, value)
fid = H5F.open(file, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
obj = H5O.open(fid, loc, 'H5P_DEFAULT');
tid = H5T.copy('H5T_C_S1');
H5T.set_size(tid, numel(value) + 1);
H5T.set_strpad(tid, 'H5T_STR_NULLTERM');
sid = H5S.create('H5S_SCALAR');
aid = H5A.create(obj, name, tid, sid, 'H5P_DEFAULT');
H5A.write(aid, tid, [value char(0)]);
H5A.close(aid); H5S.close(sid); H5T.close(tid); H5O.close(obj); H5F.close(fid);
end

function local_u64(file, loc, name, value)
fid = H5F.open(file, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
obj = H5O.open(fid, loc, 'H5P_DEFAULT');
sid = H5S.create_simple(1, 1, []);
aid = H5A.create(obj, name, 'H5T_STD_U64LE', sid, 'H5P_DEFAULT');
H5A.write(aid, 'H5ML_DEFAULT', uint64(value));
H5A.close(aid); H5S.close(sid); H5O.close(obj); H5F.close(fid);
end
