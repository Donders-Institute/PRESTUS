function out_path = export_slicer_transform(trans_pos_ras, focus_pos_ras, out_path)
% EXPORT_SLICER_TRANSFORM  Write the transducer pose as a 3D Slicer ITK transform (.tfm)
%
% Writes the pose of EXPORT_BRAINSIGHT_TRAJECTORY (location = acoustic focus,
% third axis = unit vector from the focus to the bowl apex) as an ITK affine
% transform in the convention BabelBrain reads Slicer trajectories with
% (ConvMatTransform.read_itk_affine_transform + itk_to_BSight): converting
% the file back gives that pose. The PRESTUS Nexus reads it the same way.
%
% Use as:
%   out_path = export_slicer_transform(trans_pos_ras, focus_pos_ras, out_path)
%
% Input:
%   trans_pos_ras - [1x3] bowl apex in subject RAS [mm]
%   focus_pos_ras - [1x3] acoustic focus in subject RAS [mm]
%   out_path      - file to write (typically *.tfm)
%
% Output:
%   out_path      - path written (same as input)
%
% See also: EXPORT_BRAINSIGHT_TRAJECTORY, EXPORT_KPLAN_PLAN, EXPORT_PLACEMENT

arguments
    trans_pos_ras (1,3) double
    focus_pos_ras (1,3) double
    out_path      (1,:) char
end

beam = trans_pos_ras - focus_pos_ras;
if norm(beam) < eps
    error('PRESTUS:slicer:degenerate', 'trans_pos_ras and focus_pos_ras are identical: no beam axis.');
end
z = beam / norm(beam);                       % focus -> transducer, as the Brainsight export
if abs(z(3)) < 0.9, ref = [0 0 1]; else, ref = [1 0 0]; end
x = cross(z, ref); x = x / norm(x);
y = cross(z, x);

B = [x(:) y(:) z(:) focus_pos_ras(:); 0 0 0 1];   % Brainsight/BabelBrain pose
Bn = B;  Bn(1:3,1:3) = -B(1:3,1:3);                % itk_to_BSight negates the rotation
L = diag([-1 -1 1 1]);                             % RAS <-> LPS
itk = L / Bn * L;                                  % inverse of itk_to_BSight
M = itk(1:3,1:3);  t = itk(1:3,4);

out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end
fid = fopen(out_path, 'w');
if fid == -1
    error('PRESTUS:slicer:ioError', 'Cannot write transform file:\n  %s', out_path);
end
cleanup = onCleanup(@() fclose(fid));   %#ok<NASGU>
fprintf(fid, '#Insight Transform File V1.0\n#Transform 0\n');
fprintf(fid, 'Transform: AffineTransform_double_3_3\n');
fprintf(fid, 'Parameters:%s\n', sprintf(' %.12g', [reshape(M', 1, []) t']));   % rotation row-major, then translation (as BabelBrain)
fprintf(fid, 'FixedParameters: 0 0 0\n');
end
