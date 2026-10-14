function m = time_reversal_pose_metrics(parameters, affix, target_ras, apex_ras, radius_mm)
% TIME_REVERSAL_POSE_METRICS  Target intensity and focus position of one pose-verification run
%
% Reads the outputs of a finished forward run (affix) and measures them
% against the target, which for a run with corrected aim is not its planned
% focus (so TARGETING_ANALYSIS's focus_error does not apply):
%
%   intensity_target  mean pulse-average intensity (Ipa, else Itp) within
%                     radius_mm of the target, T1w NIfTI [W/cm^2]
%   peak_*_mm         intensity peak in the -3 dB focal volume, as offset from
%                     the target: total, axial (along apex -> target, > 0
%                     deeper) and lateral (targeting CSV)
%   centroid_*_mm     the same for the centroid of the -3 dB focal volume
%
% Missing outputs give NaN and found = false.
%
% Use as:
%   m = time_reversal_pose_metrics(parameters, affix, target_ras, apex_ras, radius_mm)
%
% Input:
%   parameters - PRESTUS config (subject_id, simulation.medium, output paths)
%   affix      - io.output_affix of the run
%   target_ras, apex_ras - [1x3] target and bowl apex [RAS mm]
%   radius_mm  - radius of the target sphere [mm]
%
% Output:
%   m - struct (see above), plus intensity_quantity ('Ipa' | 'Itp')
%
% See also: POSE_WORKFLOW_PIPELINE, TARGETING_ANALYSIS

arguments
    parameters (1,1) struct
    affix      (1,:) char
    target_ras (1,3) {mustBeNumeric}
    apex_ras   (1,3) {mustBeNumeric}
    radius_mm  (1,1) {mustBeNumeric, mustBePositive}
end

    m = struct('found', false, 'intensity_target', NaN, 'intensity_quantity', "", ...
        'peak_error_mm', NaN, 'peak_axial_mm', NaN, 'peak_lateral_mm', NaN, ...
        'centroid_error_mm', NaN, 'centroid_axial_mm', NaN, 'centroid_lateral_mm', NaN);
    out = get_output_dir(parameters);
    prefix = sprintf('sub-%03d_%s', parameters.subject_id, char(parameters.simulation.medium));
    u = (target_ras - apex_ras) / norm(target_ras - apex_ras);

    for q = {'Ipa', 'Itp'}
        nii = fullfile(out, 'nii', sprintf('%s_T1w%s_%s.nii.gz', prefix, affix, q{1}));
        if ~isfile(nii), continue; end
        hdr = niftiinfo(nii);
        vol = double(niftiread(hdr));
        v0  = ras_to_grid(target_ras, hdr);
        vox = hdr.PixelDimensions(1:3);
        r   = ceil(radius_mm ./ vox);
        [i, j, k] = ndgrid(-r(1):r(1), -r(2):r(2), -r(3):r(3));
        keep = (i * vox(1)).^2 + (j * vox(2)).^2 + (k * vox(3)).^2 <= radius_mm^2;
        idx  = round(v0(:)') + [i(keep), j(keep), k(keep)];
        ok   = all(idx >= 1 & idx <= size(vol), 2);
        if any(ok)
            m.intensity_target   = mean(vol(sub2ind(size(vol), idx(ok, 1), idx(ok, 2), idx(ok, 3))));
            m.intensity_quantity = string(q{1});
            m.found = true;
        end
        break
    end

    tcsv = fullfile(out, sprintf('%s%s_targeting.csv', prefix, affix));
    if isfile(tcsv)
        T = readtable(tcsv);
        for nm = {'peak', 'centroid'}
            col = 'peak_ras';
            if strcmp(nm{1}, 'centroid'), col = 'focal_centroid_ras'; end
            if ismember([col '_x'], T.Properties.VariableNames)
                d = [T.([col '_x'])(1), T.([col '_y'])(1), T.([col '_z'])(1)] - target_ras;
                a = dot(d, u);
                m.([nm{1} '_error_mm'])   = norm(d);
                m.([nm{1} '_axial_mm'])   = a;
                m.([nm{1} '_lateral_mm']) = norm(d - a * u);
            end
        end
    end
end
