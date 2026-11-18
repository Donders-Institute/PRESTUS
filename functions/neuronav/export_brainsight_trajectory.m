function out_path = export_brainsight_trajectory(trans_pos_ras, focus_pos_ras, target_name, out_path)
% EXPORT_BRAINSIGHT_TRAJECTORY  Write a Brainsight version-14 trajectory file.
%
% Exports a tab-delimited Brainsight trajectory file (version 15 as BabelBrain writes it, NIfTI:S:<sform code>
% coordinate system) that can be imported into Brainsight for neuronavigation-guided
% transducer verification or placement confirmation.
%
% The target location (Loc.X/Y/Z) is set to focus_pos_ras — the acoustic focus
% point in RAS mm, matching what Brainsight displays as the neuronavigation target.
% The rotation matrix encodes the transducer beam orientation:
%   m2 (toolZ) = beam axis, unit vector pointing from acoustic focus → transducer rear
%   m0 (toolX) = lateral axis 1 (arbitrary, perpendicular to beam)
%   m1 (toolY) = lateral axis 2 (cross product of toolZ and toolX)
%
% Coordinate system note:
%   "NIfTI:S:<code>" names the T1's sform and its NIfTI code (1 Scanner,
%   2 Aligned, 3 Talairach, 4 MNI), i.e. the RAS mm space PRESTUS calls *_ras.
%   Pass t1_file to take the code from the image; default 'Scanner' (charm's
%   T1). BabelBrain writes 'Aligned' for every image and does not read it back.
%
% Use as:
%   out_path = export_brainsight_trajectory(trans_pos_ras, focus_pos_ras, target_name, out_path)
%   out_path = export_brainsight_trajectory(..., t1_file=t1)
%
% Input:
%   trans_pos_ras - [1×3] bowl rear-centre position in RAS space [mm]
%                  (the PRESTUS trans_pos_ras from any placement mode)
%   focus_pos_ras - [1×3] acoustic focus position in RAS space [mm]
%                  (the PRESTUS focus_pos_ras from any placement mode)
%   target_name   - label string for the trajectory row in Brainsight
%                  (spaces and underscores are fine; avoid tabs)
%   out_path      - full path to write the output file (typically *.txt)
%
% Output:
%   out_path  - path to the written file (same as input)
%
% See also: PREPROC_TRANSDUCER_PLACEMENT, POSITION_TRANSDUCER_PLANTUS,
%           POSITION_TRANSDUCER_TUS_ENTRY

arguments
    trans_pos_ras (1,3) double
    focus_pos_ras (1,3) double
    target_name   (1,:) char
    out_path      (1,:) char
    opts.t1_file  (1,:) char = ''
end

% Coordinate system label from the T1's sform code (NIfTI xform codes)
cs = 'NIfTI:S:Scanner';
if ~isempty(opts.t1_file)
    codes = {'Scanner', 'Aligned', 'Talairach', 'MNI'};
    sc = niftiinfo(opts.t1_file).raw.sform_code;
    if sc >= 1 && sc <= 4
        cs = ['NIfTI:S:' codes{sc}];
    else
        warning('PRESTUS:brainsight:noSform', '%s has no sform; labelled NIfTI:S:Scanner.', opts.t1_file);
    end
end

% Beam axis: unit vector from acoustic focus toward transducer rear centre.
% This matches the BBOP convention for toolZ (target → scalp entry direction).
beam_vec  = trans_pos_ras - focus_pos_ras;
beam_norm = norm(beam_vec);
if beam_norm < eps
    error('PRESTUS:brainsight:degenerate', ...
        'trans_pos_ras and focus_pos_ras are identical — cannot define a beam axis.');
end
toolZ = beam_vec / beam_norm;

% Build an orthonormal lateral frame (toolX, toolY) perpendicular to toolZ.
% Use Gram-Schmidt with a reference vector; pick the cardinal axis least
% collinear with toolZ to avoid near-degenerate cross products.
if abs(toolZ(3)) < 0.9
    ref = [0, 0, 1];
else
    ref = [1, 0, 0];
end
toolX = cross(toolZ, ref);
toolX = toolX / norm(toolX);
toolY = cross(toolZ, toolX);   % already unit length (both inputs are unit)

% Create parent directory if needed
out_dir = fileparts(out_path);
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

fid = fopen(out_path, 'w');
if fid == -1
    error('PRESTUS:brainsight:ioError', 'Cannot write trajectory file:\n  %s', out_path);
end
cleanup = onCleanup(@() fclose(fid));   %#ok<NASGU>  closed on function exit

fprintf(fid, '# Version: 15\n');
fprintf(fid, '# Coordinate system: %s\n', cs);
fprintf(fid, '# Created by: PRESTUS\n');
fprintf(fid, '# Units: millimetres, degrees, milliseconds, and microvolts\n');
fprintf(fid, '# Encoding: UTF-8\n');
fprintf(fid, '# Notes: Each column is delimited by a tab character.\n');
fprintf(fid, '# Target Name\tLoc. X\tLoc. Y\tLoc. Z\tm0n0\tm0n1\tm0n2\tm1n0\tm1n1\tm1n2\tm2n0\tm2n1\tm2n2\n');
fprintf(fid, '%s\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n', ...
    target_name, ...
    focus_pos_ras(1), focus_pos_ras(2), focus_pos_ras(3), ...
    toolX(1), toolX(2), toolX(3), ...
    toolY(1), toolY(2), toolY(3), ...
    toolZ(1), toolZ(2), toolZ(3));

fprintf('[export_brainsight_trajectory] Written: %s\n', out_path);

end
