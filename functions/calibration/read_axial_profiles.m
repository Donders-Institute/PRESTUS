function prof = read_axial_profiles(csv_path)
% READ_AXIAL_PROFILES  Read a hydrophone axial-profile CSV (one column per focal depth)
%
% Format (doc_transducer_library.md, Appendix): a title row, then a row whose
% columns 2..end are the focal depths [mm from EP], then one row per axial
% position with the distance from the exit plane [mm] in column 1 and the
% measured intensity [W/cm²] per focal depth. Ragged/NaN columns and rows are
% dropped; negative and non-finite intensities are set to 0.
%
% Use as:
%   prof = read_axial_profiles(csv_path)
%
% Output:
%   prof.depths_ep_mm - [1 x N_d] focal depths [mm from EP]
%   prof.z_ep_mm      - [N_z x 1] axial positions [mm from EP]
%   prof.I            - [N_z x N_d] intensity [W/cm²]
%   prof.source       - csv_path
%
% See also: FIT_TRANSDUCER_MODEL, REFINE_TRANSDUCER_LIBRARY

    raw = readmatrix(csv_path, 'NumHeaderLines', 1);
    depths = round(raw(1, 2:end), 2);
    z      = raw(2:end, 1);
    I      = raw(2:end, 2:end);
    ok = isfinite(depths);  depths = depths(ok); I = I(:, ok);
    rows = isfinite(z);     z = z(rows);         I = I(rows, :);
    I(~isfinite(I) | I < 0) = 0;
    [z, o] = sort(z); I = I(o, :);
    prof = struct('depths_ep_mm', depths, 'z_ep_mm', z, 'I', I, 'source', csv_path);
end
