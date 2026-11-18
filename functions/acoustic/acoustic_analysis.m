function [results, acoustic_Itp, acoustic_Ipa, acoustic_Ita, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, ...
                                segmentation, source_labels, reserved, provenance)
% ACOUSTIC_ANALYSIS  Compute acoustic intensity metrics and generate output plots and tables
%
% Derives spatial pressure, intensity, and mechanical index (MI) maps from
% k-Wave sensor output.  Three intensity quantities are computed:
%
%   Itp  — intensity at temporal peak: p²/(2ρc), derived from peak pressure [W/cm²]
%   Ipa  — intensity of pulse average: Itp × ramp_factor (= Itp for square pulses)
%            requires parameters.timing.pd; otherwise Ipa = Itp and saving is suppressed
%   Ita  — intensity of temporal average: Ipa × duty_cycle  [W/cm²]
%            requires parameters.timing.pd and .pri; otherwise Ita = []
%
% For layered and phantom media, computes per-tissue peak intensities (Ipa),
% pressures, and MI (skin, skull, brain) and the transcranial MI.  Writes a summary CSV and saves
% overlay plots of the intensity on the segmentation.
%
% Use as:
%   [results, acoustic_Itp, acoustic_Ipa, acoustic_Ita, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
%       acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, ...
%                         segmentation, source_labels)
%   [results, ...] = acoustic_analysis(..., [], provenance)
%
% Input:
%   parameters   - PRESTUS config with transducer, grid, analysis, io fields
%   kwave_medium - spatial medium property maps (sound_speed, density, etc.)
%   medium_masks - voxel-wise tissue layer label map
%   sensor_data  - k-Wave output; requires p_max_all [Pa] and p_min_all [Pa] (for MI)
%   segmentation - CHARM tissue label map
%   source_labels- transducer element label map (for overlay plots)
%   reserved     - unused placeholder (pass [] or omit)
%   provenance   - acoustic provenance struct from prestus_pipeline (optional)
%
% Output:
%   results          - struct with fields: Isptp, Isppa, [Ispta], Isppa_after_exit_plane,
%                      Psptp, Ipa_target, Ipa_target_radius, real_focal_distance_bowl, and
%                      per-tissue Isppa/Psptp/MI fields for layered/phantom media.
%                      ITRUSST focal metrics (Martin et al., 2025, Brain Stimulation)
%                      are computed for both layered (brain-masked) and free-field
%                      (after-exit-plane-masked) simulations:
%                        focal_mask_3dB / focal_mask_6dB  — logical focal zone masks
%                        focal_vol_3dB/6dB_mm3            — focal zone volumes [mm³]
%                        focal_centroid_3dB_bowl_mm       — distance of −3 dB centroid
%                                                           from transducer [mm]
%                        focal_fwhm_axial_mm, focal_fwhm_lateral_min/max_mm
%                                                         — −3 dB widths through the
%                                                           focal peak (compute_focal_fwhm)
%   acoustic_Itp     - intensity at temporal peak map [W/cm²]
%   acoustic_Ipa     - intensity of pulse average map [W/cm²]  (empty when timing absent)
%   acoustic_Ita     - intensity of temporal average map [W/cm²] (empty when timing absent)
%   acoustic_MI      - spatial mechanical index map [-]
%   acoustic_pressure- temporal peak pressure map [Pa]
%   highlighted_pos  - [1x3] grid coordinates of peak intensity for plots [voxels]
%
% See also: ACOUSTIC_WRAPPER, ACOUSTIC_SIMULATION, NIFTI_ACOUSTIC, MASKED_MAX_3D

arguments
    parameters   (1,1) struct
    kwave_medium (1,1) struct
    medium_masks {mustBeNumericOrLogical}
    sensor_data  (1,1) struct
    segmentation {mustBeNumericOrLogical}
    source_labels{mustBeNumericOrLogical}
    reserved     = []    %#ok<INUSA>
    provenance   struct  = struct()
end

    fprintf('[acoustic_analysis] Processing the results of acoustic simulations...\n')

    % select transducer info
    tr = parameters.transducer(1);
    trans_pos = tr.trans_pos;
    focus_pos = tr.focus_pos;

    % initialize output structure
    results = struct();

    % For async-combined caches, p_max_async holds the per-voxel pressure
    % maximum across transducers (instantaneous peak when one fires at a time).
    % p_max_all holds the thermal-equivalent field (intensity sum) and must
    % not be used for acoustic safety metrics.
    if isfield(sensor_data, 'p_max_async')
        acoustic_pressure = gather(sensor_data.p_max_async);
    else
        acoustic_pressure = gather(sensor_data.p_max_all);
    end
    results.Psptp = max(acoustic_pressure(:)); % spatial peak temporal peak pressure

    % Itp: intensity at temporal peak — p²/(2ρc) at each voxel [W/cm²].
    % Derived directly from the peak pressure field; no pulse-envelope correction.
    acoustic_Itp = acoustic_pressure.^2 ./ ...
        (2 * (kwave_medium.sound_speed .* kwave_medium.density)) .* 1e-4;

    results.Isptp = max(acoustic_Itp(:)); % spatial-peak temporal-peak intensity [W/cm²]

    % Ipa: intensity of pulse average — Itp scaled by the ramp factor (Ipa ≤ Itp).
    % For square pulses (no ramp) Ipa = Itp; ramp_f = 1.
    % Requires timing.pd to be meaningful; when timing is absent acoustic_Ipa is
    % returned as empty so that nifti_acoustic can suppress the redundant _Ipa map.
    timing_has_pd = isfield(parameters, 'timing') && isfield(parameters.timing, 'pd');
    if timing_has_pd
        ramp_f = 1;
        if isfield(parameters.timing, 'pulse_ramp_dur') && parameters.timing.pulse_ramp_dur > 0
            ramp_shape = "linear";
            if isfield(parameters.timing, 'pulse_ramp_shape')
                ramp_shape = string(parameters.timing.pulse_ramp_shape);
            end
            ramp_f = compute_ramp_factor(parameters.timing.pd, ...
                                         parameters.timing.pulse_ramp_dur, ramp_shape);
            fprintf('[acoustic_analysis] Ramp correction: Ipa = %.4f * Itp (shape: %s, ramp_dur = %.4g s)\n', ...
                    ramp_f, ramp_shape, parameters.timing.pulse_ramp_dur);
        end
        acoustic_Ipa = acoustic_Itp * ramp_f;
    else
        acoustic_Ipa = acoustic_Itp;  % fallback for analysis; not saved as a separate map
    end

    results.Isppa = max(acoustic_Ipa(:)); % spatial-peak pulse-average intensity [W/cm²]

    % Ita: intensity of temporal average — Ipa × duty_cycle [W/cm²].
    % Empty when timing.pd and .pri are not both available.
    acoustic_Ita = [];
    if timing_has_pd && isfield(parameters.timing, 'pri') && parameters.timing.pri > 0
        duty_cycle = parameters.timing.pd / parameters.timing.pri;
        acoustic_Ita = acoustic_Ipa * duty_cycle;
        results.Ispta = max(acoustic_Ita(:)); % spatial-peak temporal-average intensity [W/cm²]
        fprintf('[acoustic_analysis] Ita = %.6g * Ipa (duty_cycle = pd/pri = %.4g/%.4g s)\n', ...
                duty_cycle, parameters.timing.pd, parameters.timing.pri);
    end

    % Calculates the Mechanical Index for every gridpoint.
    % MI = |p_neg| / sqrt(f) per IEC 62359 / FDA: peak rarefaction pressure.
    % In a multi-frequency setup the lowest frequency is used (conservative).
    freqs_all  = [parameters.transducer.freq_hz];
    freq_Hz    = min(freqs_all);
    freq_MHz   = freq_Hz / 1e6;
    if numel(unique(freqs_all)) > 1
        fprintf(['[acoustic_analysis] Multi-frequency setup: MI computed using lowest frequency %g Hz ' ...
                 '(conservative; transducer frequencies: %s Hz).\n'], freq_Hz, num2str(freqs_all));
    end

    % Select peak negative pressure map (p_min_all / p_min_async are ≤ 0).
    % Guard against empty fields (e.g. jwave sets p_min_all=[] when the HDF5
    % read fails) — treat those the same as a missing field.
    if isfield(sensor_data, 'p_min_async') && ~isempty(sensor_data.p_min_async)
        acoustic_p_neg = abs(gather(sensor_data.p_min_async));
    elseif isfield(sensor_data, 'p_min_all') && ~isempty(sensor_data.p_min_all)
        acoustic_p_neg = abs(gather(sensor_data.p_min_all));
    else
        warning('acoustic_analysis:noPeakNeg', ...
            'p_min_all not found in sensor_data; falling back to p_max for MI (may underestimate).');
        acoustic_p_neg = acoustic_pressure;
    end

    acoustic_MI = (acoustic_p_neg/10^6)/sqrt(freq_MHz);
    results.MI = max(acoustic_MI(:)); % global MI, mirrors the global Isppa/Psptp convention

    fprintf('[acoustic_analysis] ---- acoustic analysis summary ----\n');
    fprintf('[acoustic_analysis] Grid dims:            %s  |  resolution: %.4g mm\n', ...
        num2str(parameters.grid.dims), parameters.grid.resolution_mm);
    fprintf('[acoustic_analysis] p_max range [Pa]:     %.4g  –  %.4g\n', ...
        min(acoustic_pressure(:)), max(acoustic_pressure(:)));
    fprintf('[acoustic_analysis] Frequency:            %.4g kHz\n', freq_Hz / 1e3);
    fprintf('[acoustic_analysis] Medium type:          %s\n', parameters.simulation.medium);
    fprintf('[acoustic_analysis] Active sensor voxels: %d\n', numel(acoustic_pressure));
    fprintf('[acoustic_analysis] ------------------------------------\n');

    % Creates a mask for voxels beyond the exit plane of any transducer.
    % For multi-transducer setups the mask is the union across transducers:
    % a voxel is "after exit" if it clears at least one transducer's exit
    % plane (the logical OR of per-transducer after-exit masks). AND would
    % yield an empty mask for bilateral / opposing-beam setups.
    comp_grid_size = size(sensor_data.p_max_all);
    after_exit_plane_mask = false(comp_grid_size);  % OR accumulation: start at 0

    for it = 1:numel(parameters.transducer)
        tr_it        = parameters.transducer(it);
        trans_pos_it = tr_it.trans_pos;
        bowl_depth_grid_it = round(annular_bowl_depth_mm(tr_it.(tr_it.type)) / parameters.grid.resolution_mm);

        % Assumption: the beam travels along the last grid dimension (z / axial),
        % and whether the transducer is on the "high-z" or "low-z" face is inferred
        % from whether trans_pos(beam_dim) > grid_size(beam_dim)/2.
        % This holds for all standard PRESTUS grid modes ('focal_axis',
        % 'transducer_axis', 'ras_plus') because grid_tissue_setup always places
        % the transducer along z.  It would NOT hold for a fully arbitrary rotation
        % applied outside of PRESTUS.
        mask_it = ones(comp_grid_size);
        if numel(parameters.grid.dims) == 3
            beam_dim = numel(parameters.grid.dims);
            if trans_pos_it(beam_dim) > comp_grid_size(beam_dim) / 2
                % HIGH-z transducer: beam propagates toward low-z.
                % Zero the near-transducer (high-z) region; keep the focus side.
                mask_it(:,:,(trans_pos_it(beam_dim) - bowl_depth_grid_it):end) = 0;
            else
                % LOW-z transducer: beam propagates toward high-z.
                % Zero the near-transducer (low-z) region; keep the focus side.
                mask_it(:,:,1:(trans_pos_it(beam_dim) + bowl_depth_grid_it)) = 0;
            end
        elseif numel(parameters.grid.dims) == 2
            beam_dim = numel(parameters.grid.dims);
            if trans_pos_it(beam_dim) > comp_grid_size(beam_dim) / 2
                mask_it(:,(trans_pos_it(beam_dim) - bowl_depth_grid_it):end) = 0;
            else
                mask_it(:,1:(trans_pos_it(beam_dim) + bowl_depth_grid_it)) = 0;
            end
        end
        % Union: voxel is "after exit" if it clears at least one transducer's exit plane.
        % AND would produce an empty mask for bilateral / opposing-beam setups.
        after_exit_plane_mask = after_exit_plane_mask | mask_it;
    end

    % [1] average Isppa within a cube around the target (half-width focus_area_radius)
    % convert the radius from mm to voxels
    avg_radius = round(parameters.analysis.focus_area_radius/parameters.grid.resolution_mm); % [voxel]
    idx = arrayfun(@(d) max(1,focus_pos(d)-avg_radius):...
        min(size(acoustic_Ipa,d),focus_pos(d)+avg_radius), ...
            1:ndims(acoustic_Ipa), 'UniformOutput', false);
    results.Ipa_target_radius = acoustic_Ipa(idx{:});
    results.Ipa_target_radius = mean(results.Ipa_target_radius(:));
    
    % [2] IPA and temporal peak pressure at the target
    idx = num2cell(focus_pos);
    results.Ipa_target  = acoustic_Ipa(idx{:});
    results.Ptp_target  = acoustic_pressure(idx{:}); % temporal peak pressure at target

    % Get tissue-specific masks. For BabelBrain runs, prefer BabelBrain's own
    % material map when it was carried through (read_babelbrain_output): it
    % matches medium_masks except at tissue boundaries, where it attributes
    % voxels (e.g. the skull-brain reflection hotspot) the way the solver did.
    if isstruct(sensor_data) && isfield(sensor_data, 'babelbrain_material_map') ...
            && ~isempty(sensor_data.babelbrain_material_map) ...
            && isequal(size(sensor_data.babelbrain_material_map), size(acoustic_Ipa))
        mask = bb_tissuemask_binary(sensor_data.babelbrain_material_map);
        fprintf(['[acoustic_analysis] Tissue masks from BabelBrain material map ' ...
            '(the exact mask the field was solved on; matches medium_masks to ' ...
            'within the label-resample round-trip).\n']);
    else
        mask = tissuemask_binary(parameters, medium_masks);
    end

    % calculate X, Y and Z coordinates of max. intensity 
    % if medium includes brain, consider only intensities within the brain
    % otherwise: consider all media beyond the exit plane
    if any(mask.brain(:))
        [results.Isppa_after_exit_plane, Ix_eplane, Iy_eplane,Iz_eplane] = ...
            masked_max_3d(acoustic_Ipa, mask.brain);
    else
        [results.Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = ...
            masked_max_3d(acoustic_Ipa, after_exit_plane_mask);
    end
    % combine coordinates into one point of max. intensity in the grid
    if numel(parameters.grid.dims)==3
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
    elseif numel(parameters.grid.dims)==2
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
    end
    fprintf('[acoustic_analysis] Final transducer, expected focus, and max ISPPA positions\n')
    
    % [3] calculate the realized focal distance in BOTH reference frames.
    % Two conventions coexist in the TUS ecosystem and both are tracked to
    % avoid ambiguity:
    %   * bowl APEX (trans_pos): BabelBrain FocalLength/ZSteering, curv_radius_mm,
    %     focal_distance_bowl, and the *_bowl_mm columns.
    %   * EXIT PLANE: manufacturer specs (FLHM centre), focal_distance_ep, and
    %     the O'Neil/Rayleigh analytical convention.
    % The exit plane sits curv_radius_mm - dist_geom_ep_mm (focal_distance_offset;
    % not the rim sagitta for library devices) in front of the apex, along the
    % beam axis toward the focus.
    if isfield(tr, tr.type)
        bowl_depth_vox = annular_bowl_depth_mm(tr.(tr.type)) / parameters.grid.resolution_mm;
    else
        bowl_depth_vox = 0;   % non-bowl source: exit plane coincides with trans_pos
    end
    beam_vec = focus_pos - trans_pos;
    if norm(beam_vec) > 0
        exit_plane_pos = trans_pos + bowl_depth_vox * beam_vec / norm(beam_vec);
    else
        exit_plane_pos = trans_pos;
    end
    % Distance helper: voxel position -> mm from a reference voxel position.
    dist_mm = @(pos, ref) norm((pos - ref) .* parameters.grid.resolution_mm);
    real_focal_distance_bowl = dist_mm(results.max_isppa_eplane_pos, trans_pos); % bowl apex
    real_focal_distance_ep   = dist_mm(results.max_isppa_eplane_pos, exit_plane_pos); % exit plane

    % [3b] cranial Thermal Index (IEC 62359), informational. The emitted-power
    % derivation is an explicit assumption isolated in compute_tic.m.
    results.TIC = compute_tic(parameters, tr, results.Isppa);

    % Layer-specific outcomes (in case a layered simulation)
    if contains(parameters.simulation.medium, {'layered'; 'phantom'})

        % extract indices in brain medium
        [results.Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = ...
            masked_max_3d(acoustic_Ipa, mask.brain);
        [results.min_Isppa_brain] = min(acoustic_Ipa(mask.brain));

        % Position highlighted in plots (blue marker): the BRAIN-masked peak
        % when brain tissue exists, consistent with Isppa_after_exit_plane /
        % real_focal_distance_bowl above (both already brain-masked) and with the
        % Isppa_brain value reported alongside the plots. The unmasked global
        % peak is NOT used: for elastic-skull solvers (BabelBrain) it lands in
        % the skull-brain interface reflection hotspot, not at the focal spot
        % of interest, and would drag the marker (and, below, the colour
        % scale) off the brain target.
        if any(mask.brain(:))
            highlighted_pos = [Ix_brain, Iy_brain, Iz_brain];
        else
            [~, Ix, Iy, Iz] = masked_max_3d(acoustic_Ipa, ones(size(medium_masks)));
            highlighted_pos = [Ix, Iy, Iz];
        end
        highlighted_pos = highlighted_pos(1:numel(parameters.grid.dims));

        % ITRUSST focal metrics (Martin et al., 2025, Brain Stimulation):
        %   −3 dB: I ≥ Isppa/2   (half-intensity; 10·log10(0.5) = −3 dB)
        %   −6 dB: I ≥ Isppa/4   (half-pressure; P ≥ Pmax/2; 10·log10(0.25) = −6 dB)
        fm = compute_focal_metrics(acoustic_Ipa, results.Isppa_brain, mask.brain, ...
                                   trans_pos, parameters.grid.resolution_mm, ...
                                   focal_reference_from_parameters(parameters));
        results.focal_mask_3dB        = fm.mask_3dB;
        results.focal_mask_6dB        = fm.mask_6dB;
        results.focal_vol_3dB_mm3          = fm.vol_3dB_mm3;
        results.focal_vol_6dB_mm3          = fm.vol_6dB_mm3;
        % Distances from both references (bowl apex, exit plane) from the frame-agnostic
        % *_vox positions returned by compute_focal_metrics.
        results.focal_centroid_3dB_bowl_mm = fm.centroid_3dB_bowl_mm;                    % bowl apex
        results.focal_centroid_3dB_ep_mm   = dist_mm(fm.centroid_3dB_vox, exit_plane_pos);
        results.focal_flhm_center_bowl_mm  = dist_mm(fm.flhm_center_vox, trans_pos);
        results.focal_flhm_center_ep_mm    = dist_mm(fm.flhm_center_vox, exit_plane_pos);
        if fm.flhm_truncated
            warning('acoustic_analysis:flhmTruncated', ...
                ['−3 dB focal lobe reaches the grid boundary; ' ...
                 'focal_flhm_center_*_mm is unreliable (enlarge the FOV).']);
        end

        % −3 dB focal widths through the brain peak, along and across the
        % transducer axis (ITRUSST reporting item)
        fw = compute_focal_fwhm(acoustic_Ipa, fm.peak_pos_vox, beam_vec, ...
                                parameters.grid.resolution_mm);
        results.focal_fwhm_axial_mm       = fw.axial_mm;
        results.focal_fwhm_lateral_min_mm = fw.lateral_min_mm;
        results.focal_fwhm_lateral_max_mm = fw.lateral_max_mm;

        % Keep backward-compatible alias (−3 dB volume)
        half_max_ISPPA_volume_brain = fm.vol_3dB_mm3;

        [results.Psptp_brain] = masked_max_3d(acoustic_pressure, mask.brain);
        [results.MI_brain] = masked_max_3d(acoustic_MI, mask.brain);

        % extract indices in skull medium
        [results.Isppa_skull] = masked_max_3d(acoustic_Ipa, mask.skull);
        [results.Psptp_skull] = masked_max_3d(acoustic_pressure, mask.skull);
        [results.MI_skull] = masked_max_3d(acoustic_MI, mask.skull);

        % extract indices in skin medium
        [results.Isppa_skin] = masked_max_3d(acoustic_Ipa, mask.skin);
        [results.Psptp_skin] = masked_max_3d(acoustic_pressure, mask.skin);
        [results.MI_skin] = masked_max_3d(acoustic_MI, mask.skin);

        % MItc: transcranial MI — max MI across intracranial voxels (WM, GM, CSF, blood)
        intracranial_mask = ismember(segmentation, charm_seg_labels().intracranial);
        if any(intracranial_mask(:))
            [results.MI_tc] = masked_max_3d(acoustic_MI, intracranial_mask);
        else
            results.MI_tc = NaN; % intracranial tissues not modelled
        end

        % Free-field ISPPA from provenance (BabelBrain water run, or NaN otherwise).
        % Enables normalized comparison across solvers: Isppa_brain / freefield_isppa_wcm2.
        freefield_isppa_wcm2 = NaN;
        if isfield(provenance, 'freefield_isppa_wcm2') && isfinite(provenance.freefield_isppa_wcm2)
            freefield_isppa_wcm2 = provenance.freefield_isppa_wcm2;
            results.freefield_isppa_wcm2 = freefield_isppa_wcm2;
        end

        Isptp_val = results.Isptp;
        Ispta_val = NaN; if isfield(results, 'Ispta'), Ispta_val = results.Ispta; end

        writetable(table(parameters.subject_id, freq_Hz, ...
            Isptp_val, results.Isppa, Ispta_val, results.Isppa_after_exit_plane, ...
            real_focal_distance_bowl, results.Isppa_skin, results.Isppa_skull, ...
            results.Isppa_brain, results.Psptp_skin, results.Psptp_skull, ...
            results.Psptp_brain, results.Psptp, results.Ptp_target, results.MI_skin, results.MI_skull, ...
            results.MI_brain, results.MI_tc, Ix_brain, Iy_brain, Iz_brain, trans_pos, focus_pos, ...
            results.Ipa_target, results.Ipa_target_radius, half_max_ISPPA_volume_brain, ...
            results.focal_vol_3dB_mm3, results.focal_vol_6dB_mm3, ...
            results.focal_centroid_3dB_bowl_mm, results.focal_centroid_3dB_ep_mm, ...
            real_focal_distance_ep, ...
            results.focal_flhm_center_bowl_mm, results.focal_flhm_center_ep_mm, ...
            freefield_isppa_wcm2, ...
            results.TIC, ...
            results.focal_fwhm_axial_mm, results.focal_fwhm_lateral_min_mm, ...
            results.focal_fwhm_lateral_max_mm, ...
            'VariableNames', { 'subject_id', 'freq_Hz', ...
            'Isptp', 'Isppa', 'Ispta', 'Isppa_after_exitplane', 'real_focal_distance_bowl_mm', ...
            'Isppa_skin', 'Isppa_skull', 'Isppa_brain', ...
            'Psptp_skin', 'Psptp_skull', 'Psptp_brain', 'Psptp', 'Ptp_target', ...
            'MI_skin', 'MI_skull', 'MI_brain', 'MI_tc', ...
            'Ix_brain_vox', 'Iy_brain_vox', 'Iz_brain_vox', ...
            'trans_pos_vox', 'focus_pos_vox', ...
            'Ipa_target', 'Ipa_target_radius', ...
            'halfmax_ISPPA_volume_brain_mm3', ...
            'focal_vol_3dB_mm3', 'focal_vol_6dB_mm3', ...
            'focal_centroid_3dB_bowl_mm', 'focal_centroid_3dB_ep_mm', ...
            'real_focal_distance_ep_mm', ...
            'focal_flhm_center_bowl_mm', 'focal_flhm_center_ep_mm', ...
            'freefield_isppa_wcm2', 'TIC', ...
            'focal_fwhm_axial_mm', 'focal_fwhm_lateral_min_mm', 'focal_fwhm_lateral_max_mm'}), ...
            parameters.io.filename_table);
    else
        % Free-field simulation: compute ITRUSST focal metrics using after-exit-plane mask.
        % Region mask restricts the focal zone to beyond the transducer exit plane.
        fm_ff = compute_focal_metrics(acoustic_Ipa, results.Isppa_after_exit_plane, ...
                                      after_exit_plane_mask, trans_pos, parameters.grid.resolution_mm, ...
                                      focal_reference_from_parameters(parameters));
        results.focal_mask_3dB        = fm_ff.mask_3dB;
        results.focal_mask_6dB        = fm_ff.mask_6dB;
        results.focal_vol_3dB_mm3          = fm_ff.vol_3dB_mm3;
        results.focal_vol_6dB_mm3          = fm_ff.vol_6dB_mm3;
        % Distances from both references (see layered branch).
        results.focal_centroid_3dB_bowl_mm = fm_ff.centroid_3dB_bowl_mm;                    % bowl apex
        results.focal_centroid_3dB_ep_mm   = dist_mm(fm_ff.centroid_3dB_vox, exit_plane_pos);
        results.focal_flhm_center_bowl_mm  = dist_mm(fm_ff.flhm_center_vox, trans_pos);
        results.focal_flhm_center_ep_mm    = dist_mm(fm_ff.flhm_center_vox, exit_plane_pos);
        if fm_ff.flhm_truncated
            warning('acoustic_analysis:flhmTruncated', ...
                ['−3 dB focal lobe reaches the grid boundary; ' ...
                 'focal_flhm_center_*_mm is unreliable (enlarge the FOV).']);
        end
        fw = compute_focal_fwhm(acoustic_Ipa, fm_ff.peak_pos_vox, beam_vec, ...
                                parameters.grid.resolution_mm);
        results.focal_fwhm_axial_mm       = fw.axial_mm;
        results.focal_fwhm_lateral_min_mm = fw.lateral_min_mm;
        results.focal_fwhm_lateral_max_mm = fw.lateral_max_mm;

        highlighted_pos = results.max_isppa_eplane_pos;
        Isptp_val = results.Isptp;
        Ispta_val = NaN; if isfield(results, 'Ispta'), Ispta_val = results.Ispta; end
        writetable(table(parameters.subject_id, freq_Hz, ...
            Isptp_val, results.Isppa, Ispta_val, results.Isppa_after_exit_plane, ...
            results.Psptp, results.MI, results.Ptp_target, real_focal_distance_bowl, trans_pos, focus_pos, ...
            results.Ipa_target, results.Ipa_target_radius, ...
            results.focal_vol_3dB_mm3, results.focal_vol_6dB_mm3, ...
            results.focal_centroid_3dB_bowl_mm, results.focal_centroid_3dB_ep_mm, ...
            real_focal_distance_ep, ...
            results.focal_flhm_center_bowl_mm, results.focal_flhm_center_ep_mm, ...
            results.TIC, ...
            results.focal_fwhm_axial_mm, results.focal_fwhm_lateral_min_mm, ...
            results.focal_fwhm_lateral_max_mm, ...
            'VariableNames', {'subject_id', 'freq_Hz', ...
            'Isptp', 'Isppa', 'Ispta', 'Isppa_after_exitplane', ...
            'Psptp', 'MI', 'Ptp_target', 'real_focal_distance_bowl_mm', 'trans_pos_vox', 'focus_pos_vox', ...
            'Ipa_target', 'Ipa_target_radius', ...
            'focal_vol_3dB_mm3', 'focal_vol_6dB_mm3', ...
            'focal_centroid_3dB_bowl_mm', 'focal_centroid_3dB_ep_mm', ...
            'real_focal_distance_ep_mm', ...
            'focal_flhm_center_bowl_mm', 'focal_flhm_center_ep_mm', 'TIC', ...
            'focal_fwhm_axial_mm', 'focal_fwhm_lateral_min_mm', 'focal_fwhm_lateral_max_mm'}), ...
        parameters.io.filename_table);
    end

    % Plot intensity on the segmented image (up to 2 transducers).
    % For multi-transducer setups all transducer geometries are drawn on
    % every slice; slices are centred on each transducer's focus in turn.
    n_plots = min(2, numel(parameters.transducer));
    if numel(parameters.transducer) > n_plots
        warning('More than two transducers: intensity plots on segmentation will be created only for the first 2 transducers');
    end

    % Build N×3 position matrices so plot_overlay can draw every transducer.
    all_tpos = zeros(n_plots, numel(parameters.transducer(1).trans_pos));   % 2 on 2D grids
    all_fpos = zeros(n_plots, numel(parameters.transducer(1).focus_pos));
    for ti = 1:n_plots
        all_tpos(ti, :) = parameters.transducer(ti).trans_pos;
        all_fpos(ti, :) = parameters.transducer(ti).focus_pos;
    end

    % Sync mode: all transducers fire together — one combined plot, no trans-XX suffix.
    sync_mode = n_plots > 1 && ~is_async_mode(parameters);
    n_iters   = 1 + (~sync_mode) * (n_plots - 1);  % 1 for sync, n_plots for async/single

    % Colour range for the intensity overlays: cap at the brain-masked peak
    % (Isppa_brain) rather than the unmasked global max. For elastic-skull
    % solvers (BabelBrain) the global max sits in the skull-brain interface
    % reflection hotspot — several times the brain focal value — which would
    % compress the entire brain field into the bottom of the colour scale.
    % Matches the convention already used by plot_intensity_t1_overlay.
    overlay_cmax = [];
    if isfield(results, 'Isppa_brain') && ~isempty(results.Isppa_brain) ...
            && isfinite(results.Isppa_brain) && results.Isppa_brain > 0
        overlay_cmax = results.Isppa_brain;
    elseif isfield(results, 'Isppa') && isfinite(results.Isppa) && results.Isppa > 0
        overlay_cmax = results.Isppa;
    end

    for ti = 1:n_iters

        fpos_sim = parameters.transducer(ti).focus_pos;

        t_suffix = trans_suffix(ti, n_iters);

        if numel(parameters.grid.dims)==3
            slices = struct( ...
                'dim', {'x', 'y', 'z'}, ...
                'pos', {fpos_sim(1), fpos_sim(2), fpos_sim(3)} ...
                );

            for i_slice = 1:numel(slices)
                slice = slices(i_slice);

                overlay_opts = {};
                if ~isempty(overlay_cmax)
                    overlay_opts = {'overlay_color_range', [0 overlay_cmax], ...
                                    'overlay_threshold_low', 0, ...
                                    'overlay_threshold_high', overlay_cmax};
                end

                [~,~,~,~,~,~,~,h]=plot_overlay(...
                    acoustic_Ipa, ...
                    segmentation, ...
                    source_labels, ...
                    parameters, ...
                    {slice.dim, slice.pos}, ...
                    all_tpos, ...
                    all_fpos, ...
                    highlighted_pos, ...
                    overlay_opts{:});

                % Construct output filename
                output_plot = fullfile( ...
                    parameters.io.dir_img, ...
                    sprintf('sub-%03d_%s_intensity_%s%s%s.png', ...
                    parameters.subject_id, ...
                    parameters.simulation.medium, ...
                    slice.dim, ...
                    t_suffix, ...
                    parameters.io.output_affix) ...
                    );

                % Keep original colors and save
                set(h, 'InvertHardcopy', 'off');
                saveas(h, output_plot, 'png');
                close(h);
            end
        else
            h = plot_overlay_2d(...
                acoustic_Ipa, ...
                segmentation, ...
                source_labels, ...
                after_exit_plane_mask, ...
                all_tpos(ti,:), ...
                all_fpos(ti,:), ...
                highlighted_pos);

            output_plot = fullfile(parameters.io.dir_img, ...
                sprintf('sub-%03d_%s_intensity%s%s.png', ...
                parameters.subject_id, ...
                parameters.simulation.medium, ...
                t_suffix, ...
                parameters.io.output_affix));
            set(h, 'InvertHardcopy', 'off'); % keep original colours
            saveas(h, output_plot, 'png')
            close(h);
        end
    end

    % Axial profile (3-D only): pressure and intensity along beam axis through focus,
    % with tissue-layer shading. Saves PNG + CSV alongside the overlay plots.
    if numel(parameters.grid.dims) == 3 && ~isempty(medium_masks)
        try
            % Free-water references at the plotted drive level; not for combined
            % (per-source scaled) fields or BabelBrain (own drive scaling).
            ref_opts = struct('IsppaScale', 1, 'FreefieldAxial', struct(), ...
                'FreeWater', ~isfield(provenance, 'combination_mode') && ...
                ~(isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'babelbrain')));
            if isfield(sensor_data, 'isppa_scale_p')
                ref_opts.IsppaScale = double(sensor_data.isppa_scale_p)^2;
            end
            if isfield(provenance, 'freefield_axial') && ~isempty(provenance.freefield_axial)
                ref_opts.FreefieldAxial = provenance.freefield_axial(1);
            end
            ref_args = namedargs2cell(ref_opts);
            plot_axial_profile(acoustic_Ipa, acoustic_pressure, medium_masks, ...
                parameters, exit_plane_pos, focus_pos, highlighted_pos, ref_args{:});
        catch ME
            warning('acoustic_analysis:axialProfile', ...
                'plot_axial_profile failed: %s', ME.message);
        end
    end
end

% -------------------------------------------------------------------------
function focal_reference = focal_reference_from_parameters(parameters)
% Defensive lookup of calibration.focal_reference, matching
% COMPUTE_FLHM_CENTER's own default so callers that construct parameters
% by hand (rather than via load_parameters, which always populates this
% from config_default.yaml) still get sensible behaviour.
    if isfield(parameters, 'calibration') && isfield(parameters.calibration, 'focal_reference')
        focal_reference = parameters.calibration.focal_reference;
    else
        focal_reference = 'flhm_3db';
    end
end

% -------------------------------------------------------------------------
function s = vec2csvstr(v)
% Format a numeric row vector as a bracketed comma-separated string, e.g.
% [1.23,4.56].  Used to store per-transducer vectors in a single CSV cell.
    s = ['[' strjoin(arrayfun(@(x) sprintf('%.6g', x), v, 'UniformOutput', false), ',') ']'];
end


