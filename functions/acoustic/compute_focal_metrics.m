function fm = compute_focal_metrics(acoustic_Ipa, peak_Isppa, region_mask, trans_pos, res_mm)
% COMPUTE_FOCAL_METRICS  −3 dB and −6 dB focal zone metrics per ITRUSST reporting
%
% Computes the size and position of the acoustic focal zone at two standard
% thresholds, following Martin et al. (2025, Brain Stimulation) ITRUSST
% reporting guidelines:
%
%   −3 dB  I ≥ Isppa_peak/2    (half-intensity; called "−3 dB" because
%                                10·log10(0.5) = −3 dB)
%   −6 dB  I ≥ Isppa_peak/4    (= P ≥ Pmax/2; half-pressure amplitude;
%                                called "−6 dB" because 10·log10(0.25) = −6 dB
%                                and 20·log10(0.5) = −6 dB — both consistent)
%
% ITRUSST requires reporting for free-field characterisation:
%   - Position of spatial-peak pressure amplitude (relative to transducer)
%   - −3 dB and −6 dB focal volumes
%   - Position of centre of −3 dB focal volume (relative to transducer)
%
% Use as:
%   fm = compute_focal_metrics(acoustic_Ipa, peak_Isppa, region_mask, ...
%                              trans_pos, res_mm)
%
% Input:
%   acoustic_Ipa  - (:,:,:) intensity map [W/cm²] in simulation grid
%   peak_Isppa    - scalar peak intensity within region_mask used as reference
%                   (e.g. results.Isppa_brain or results.Isppa_after_exit_plane)
%   region_mask   - logical (:,:,:), valid voxels for focal zone computation
%                   (brain mask for layered; after-exit-plane mask for free-field)
%   trans_pos     - [1×3] transducer position in grid voxels (reference point)
%   res_mm        - scalar grid resolution [mm/voxel]
%
% Output:
%   fm  struct with fields (all linear dimensions in mm, volumes in mm³):
%     peak_pos_vox          [1×3] voxel coords of spatial-peak intensity
%     peak_dist_mm          scalar distance of peak from transducer [mm]
%     mask_3dB              logical (:,:,:) — I ≥ Isppa/2 within region_mask
%     mask_6dB              logical (:,:,:) — I ≥ Isppa/4 within region_mask
%     vol_3dB_mm3           volume of −3 dB region [mm³]
%     vol_6dB_mm3           volume of −6 dB region [mm³]
%     centroid_3dB_vox      [1×3] voxel coords of centroid of −3 dB region
%     centroid_3dB_dist_mm  distance from transducer to −3 dB centroid [mm]
%
% See also: ACOUSTIC_ANALYSIS, MASKED_MAX_3D

arguments
    acoustic_Ipa  {mustBeNumeric}
    peak_Isppa    (1,1) double
    region_mask   {mustBeNumericOrLogical}
    trans_pos     (1,:) double
    res_mm        (1,1) double
end

    fm = struct();

    % Find peak intensity voxel within region
    [~, Ix, Iy, Iz] = masked_max_3d(acoustic_Ipa, region_mask);
    fm.peak_pos_vox   = [Ix, Iy, Iz];
    fm.peak_dist_mm   = norm(([Ix, Iy, Iz] - trans_pos) .* res_mm);

    % Build −3 dB and −6 dB masks
    fm.mask_3dB = acoustic_Ipa >= peak_Isppa / 2 & region_mask;  % I ≥ Imax/2
    fm.mask_6dB = acoustic_Ipa >= peak_Isppa / 4 & region_mask;  % I ≥ Imax/4 (FWHM)

    % Volumes
    fm.vol_3dB_mm3 = sum(fm.mask_3dB(:)) * res_mm^3;
    fm.vol_6dB_mm3 = sum(fm.mask_6dB(:)) * res_mm^3;

    % Centroid of −3 dB region (geometric centre of binary mask)
    idx = find(fm.mask_3dB);
    if ~isempty(idx)
        [cx, cy, cz] = ind2sub(size(fm.mask_3dB), idx);
        fm.centroid_3dB_vox     = [mean(cx), mean(cy), mean(cz)];
    else
        fm.centroid_3dB_vox     = fm.peak_pos_vox;
    end
    fm.centroid_3dB_dist_mm = norm((fm.centroid_3dB_vox - trans_pos) .* res_mm);
end
