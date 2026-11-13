function fm = compute_focal_metrics(acoustic_Ipa, peak_Isppa, region_mask, trans_pos, res_mm, focal_reference)
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
% The −3 dB/−6 dB masks and volumes above are always reported at both fixed
% ITRUSST thresholds regardless of focal_reference. Only the "FLHM centre"
% (flhm_center_*, flhm_span_mm) -- the quantity comparable to a manufacturer's
% quoted "focal distance" -- uses the configurable threshold, since different
% manufacturers define that quantity differently (see FOCAL_REFERENCE_TO_FRAC).
%
% Use as:
%   fm = compute_focal_metrics(acoustic_Ipa, peak_Isppa, region_mask, ...
%                              trans_pos, res_mm)
%   fm = compute_focal_metrics(acoustic_Ipa, peak_Isppa, region_mask, ...
%                              trans_pos, res_mm, focal_reference)
%
% Input:
%   acoustic_Ipa  - (:,:,:) intensity map [W/cm²] in simulation grid
%   peak_Isppa    - scalar peak intensity within region_mask used as reference
%                   (e.g. results.Isppa_brain or results.Isppa_after_exit_plane)
%   region_mask   - logical (:,:,:), valid voxels for focal zone computation
%                   (brain mask for layered; after-exit-plane mask for free-field)
%   trans_pos     - [1×3] reference point in grid voxels for the *_dist_mm
%                   outputs.  acoustic_analysis passes the bowl apex, so the
%                   *_dist_mm fields are apex-referenced (legacy convention).
%                   Exit-plane distances are derived by the caller from the
%                   returned *_vox positions (which are frame-agnostic).
%   res_mm        - scalar grid resolution [mm/voxel]
%   focal_reference - convention for the FLHM-centre threshold: 'flhm_3db'
%                   (default, I >= peak/2) | 'flhm_6db' (I >= peak/4) |
%                   'peak' (centre = peak position; span at -3 dB).
%                   Should match calibration.focal_reference for the
%                   equipment being simulated. See FOCAL_REFERENCE_TO_FRAC.
%
% Output:
%   fm  struct with fields (all linear dimensions in mm, volumes in mm³):
%     peak_pos_vox          [1×3] voxel coords of spatial-peak intensity
%     peak_dist_mm          scalar distance of peak from reference point [mm]
%     mask_3dB              logical (:,:,:) — I ≥ Isppa/2 within region_mask
%     mask_6dB              logical (:,:,:) — I ≥ Isppa/4 within region_mask
%     vol_3dB_mm3           volume of −3 dB region [mm³]
%     vol_6dB_mm3           volume of −6 dB region [mm³]
%     centroid_3dB_vox      [1×3] voxel coords of centroid of −3 dB region
%     centroid_3dB_dist_mm  distance from reference point to −3 dB centroid [mm]
%     flhm_center_vox       [1×3] voxel coords of the FLHM centre (on-axis)
%     flhm_center_dist_mm   distance from reference point to the FLHM centre [mm]
%     flhm_span_mm          axial focal length at half maximum, at the
%                           configured focal_reference threshold [mm]
%     flhm_truncated        true if the FLHM span hit the grid boundary
%                           (FLHM centre then unreliable — enlarge the FOV)
%
% FLHM (focal length at half maximum) centre convention:
%   The through-focus axial intensity profile is thresholded at the
%   configured focal_reference fraction of peak; the midpoint of the
%   contiguous span bracketing the peak is the FLHM centre. This is what
%   manufacturers usually quote as the "focal distance" (from the exit
%   plane), at a manufacturer-specific threshold: Imasonic/IGT steering
%   tables use −6 dB (I >= peak/4) rather than the −3 dB ITRUSST default.
%   The FLHM centre lies distal to the intensity peak, because the focal
%   lobe is skewed distally. With focal_reference = 'peak' the centre is
%   the peak itself.
%
% See also: ACOUSTIC_ANALYSIS, MASKED_MAX_3D, FOCAL_REFERENCE_TO_FRAC

arguments
    acoustic_Ipa    {mustBeNumeric}
    peak_Isppa      (1,1) double
    region_mask     {mustBeNumericOrLogical}
    trans_pos       (1,:) double
    res_mm          (1,1) double
    focal_reference (1,:) char = 'flhm_3db'
end

    % 'peak' has no FLHM threshold of its own: the centre is the peak itself,
    % and the reported span falls back to the ITRUSST -3 dB threshold.
    use_peak = strcmp(focal_reference, 'peak');
    if use_peak
        flhm_thr_frac = focal_reference_to_frac('flhm_3db');
    else
        flhm_thr_frac = focal_reference_to_frac(focal_reference);
    end

    fm = struct();

    % Find peak intensity voxel within region
    [~, Ix, Iy, Iz] = masked_max_3d(acoustic_Ipa, region_mask);
    % Voxel positions have the grid's dimensionality (2 on 2D/axisymmetric
    % grids, where masked_max_3d returns a trailing index of 1)
    nd = ndims(acoustic_Ipa);
    fm.peak_pos_vox   = [Ix, Iy, Iz];
    fm.peak_pos_vox   = fm.peak_pos_vox(1:nd);
    fm.peak_dist_mm   = norm((fm.peak_pos_vox - trans_pos) .* res_mm);

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
        fm.centroid_3dB_vox     = fm.centroid_3dB_vox(1:nd);
    else
        fm.centroid_3dB_vox     = fm.peak_pos_vox;
    end
    fm.centroid_3dB_dist_mm = norm((fm.centroid_3dB_vox - trans_pos) .* res_mm);

    % FLHM centre along the beam axis (last grid dimension).  Extract the
    % on-axis line through the peak, walk outward from the peak while the
    % intensity stays above the focal_reference threshold, and take the
    % midpoint of that contiguous span.
    beam_dim = ndims(acoustic_Ipa);
    if beam_dim == 3
        axial_line = squeeze(acoustic_Ipa(Ix, Iy, :));
        pk_idx     = Iz;
    else   % 2D grid: beam along the 2nd dimension
        axial_line = squeeze(acoustic_Ipa(Ix, :))';
        pk_idx     = Iy;
    end
    n_ax = numel(axial_line);
    thr  = peak_Isppa * flhm_thr_frac;       % per configured focal_reference
    lo = pk_idx;
    while lo > 1 && axial_line(lo - 1) >= thr
        lo = lo - 1;
    end
    hi = pk_idx;
    while hi < n_ax && axial_line(hi + 1) >= thr
        hi = hi + 1;
    end
    if use_peak
        mid_idx = pk_idx;
    else
        mid_idx = (lo + hi) / 2;
    end

    fm.flhm_center_vox = fm.peak_pos_vox;
    fm.flhm_center_vox(beam_dim) = mid_idx;
    fm.flhm_center_dist_mm = norm((fm.flhm_center_vox - trans_pos) .* res_mm);
    fm.flhm_span_mm        = (hi - lo) * res_mm;
    fm.flhm_truncated      = (lo == 1) || (hi == n_ax);
end
