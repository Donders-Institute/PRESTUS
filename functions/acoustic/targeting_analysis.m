function tgt = targeting_analysis(parameters, options)
% TARGETING_ANALYSIS  Targeting accuracy of one run from its focal-mask NIfTIs
%
% Compares the simulated -3 dB focal volume with the intended target, in
% native (T1 RAS mm) and MNI space, and writes one row per run to
%
%   <dir_tabular>/sub-XXX_<medium><affix>_targeting.csv
%
% which the group report joins with the per-subject CSV. Measures:
%
%   Native (T1w NIfTIs, always when the -3 dB mask exists)
%     focal_centroid_ras_*   centroid of the -3 dB focal volume [RAS mm]
%     peak_ras_*             intensity peak inside that volume [RAS mm]
%     planned_focus_ras_*    focus position resolved by placement [RAS mm]
%     focus_error_mm         |centroid - planned focus|
%     focus_error_axial_mm   component along the beam (transducer -> planned
%                            focus); > 0 = simulated focus deeper than planned
%     focus_error_lateral_mm component perpendicular to the beam
%     peak_error_mm          |peak - planned focus|
%   MNI (MNI NIfTIs, when io.save_MNI produced them)
%     focal_centroid_mni_*, peak_mni_*
%     target_mni_*           analysis.target_mni_mm, else the MNI target of the
%                            placement mode (target_mni_source says which)
%     target_error_mni_mm    |MNI centroid - MNI target|
%     roi_*                  with analysis.target_roi_mni: ROI volume, overlap
%                            with the focal volume, fraction of the focal
%                            volume inside the ROI, ROI coverage, mean/max
%                            intensity in the ROI, centroid-to-ROI-centroid
%
% With several transducers, the planned focus and beam axis of
% transducer(1) are used. Missing inputs give NaN, never an error.
%
% Use as:
%   tgt = targeting_analysis(parameters)
%   tgt = targeting_analysis(parameters, 'PlannedFocus', f, 'PlannedTrans', t)
%
% Input:
%   parameters - PRESTUS config of the run (after path_log_setup). The planned
%                T1 positions are read from parameters.state.planned_focus_t1_vox
%                / planned_trans_t1_vox, recorded by PRESTUS_PIPELINE before
%                the grid stage (in the RAS+ T1 grid of the T1w outputs; see
%                GRID_TISSUE_SETUP).
%   options    - name-value pairs:
%     'PlannedFocus' [N×3] planned focus per transducer [T1 voxels, 1-based]
%     'PlannedTrans' [N×3] planned transducer position [T1 voxels, 1-based]
%     'Write'        write the CSV (default true)
%
% Output:
%   tgt - one-row table (also written as CSV)
%
% See also: STUDY_RUN_DESCRIPTORS, SAMPLE_VOLUME_AT_WORLD, GRID_TO_RAS,
%           GENERATE_GROUP_REPORT

arguments
    parameters (1,1) struct
    options.PlannedFocus double = []
    options.PlannedTrans double = []
    options.Write (1,1) logical = true
end

    id     = parameters.subject_id;
    medium = char(parameters.simulation.medium);
    affix  = char(parameters.io.output_affix);
    [dir_nii, dir_tab] = output_dirs(parameters);

    focus_vox = options.PlannedFocus;
    trans_vox = options.PlannedTrans;
    if isempty(focus_vox) && isfield(parameters, 'state') && isfield(parameters.state, 'planned_focus_t1_vox')
        focus_vox = parameters.state.planned_focus_t1_vox;
    end
    if isempty(trans_vox) && isfield(parameters, 'state') && isfield(parameters.state, 'planned_trans_t1_vox')
        trans_vox = parameters.state.planned_trans_t1_vox;
    end

    d = study_run_descriptors(parameters);
    row = struct();
    row.subject_id        = id;
    row.condition         = d.condition;
    row.affix             = string(affix);
    row.placement_mode    = d.placement_mode;
    row.target_label      = d.target_label;
    row.transducer_name   = d.transducer_name;
    row.n_transducers     = d.n_transducers;
    row.focal_distance_ep_planned_mm = first_or_nan(d.focal_distance_ep_mm);
    row.target_isppa_wcm2 = first_or_nan(d.target_isppa_wcm2);

    nii = @(space, type) fullfile(dir_nii, sprintf('sub-%03d_%s_%s%s_%s.nii.gz', id, medium, space, affix, type));

    %% Native space
    [cen, peak, vol_mm3, iq] = focal_summary(nii('T1w', 'focal_mask_3dB'), ...
        {nii('T1w', 'Ipa'), nii('T1w', 'Itp')});
    row.intensity_quantity = string(iq);
    row.focal_vol_3dB_T1w_mm3 = vol_mm3;
    row = put_xyz(row, 'focal_centroid_ras', cen);
    row = put_xyz(row, 'peak_ras', peak);

    planned_focus = nan(1, 3); planned_trans = nan(1, 3);
    t1_hdr = header_of(nii('T1w', 'focal_mask_3dB'));
    if ~isempty(t1_hdr) && size(focus_vox, 2) == 3 && ~isempty(focus_vox)
        planned_focus = grid_to_ras(focus_vox(1, :), t1_hdr);
        if size(trans_vox, 2) == 3 && ~isempty(trans_vox)
            planned_trans = grid_to_ras(trans_vox(1, :), t1_hdr);
        end
    end
    row = put_xyz(row, 'planned_focus_ras', planned_focus);
    row = put_xyz(row, 'planned_trans_ras', planned_trans);

    err = cen - planned_focus;
    row.focus_error_mm = norm_or_nan(err);
    beam = planned_focus - planned_trans;
    if all(isfinite(beam)) && norm(beam) > 0 && all(isfinite(err))
        u = beam / norm(beam);
        row.focus_error_axial_mm   = dot(err, u);
        row.focus_error_lateral_mm = norm(err - dot(err, u) * u);
    else
        row.focus_error_axial_mm   = NaN;
        row.focus_error_lateral_mm = NaN;
    end
    row.peak_error_mm = norm_or_nan(peak - planned_focus);

    %% MNI space
    mni_mask = nii('MNI', 'focal_mask_3dB');
    [cen_mni, peak_mni] = focal_summary(mni_mask, {nii('MNI', 'Ipa'), nii('MNI', 'Itp')});
    row = put_xyz(row, 'focal_centroid_mni', cen_mni);
    row = put_xyz(row, 'peak_mni', peak_mni);

    [target_mni, source] = resolve_target_mni(parameters);
    row = put_xyz(row, 'target_mni', target_mni);
    row.target_mni_source = string(source);
    row.target_error_mni_mm      = norm_or_nan(cen_mni - target_mni);
    row.peak_target_error_mni_mm = norm_or_nan(peak_mni - target_mni);

    row = roi_metrics(row, parameters, mni_mask, {nii('MNI', 'Ipa'), nii('MNI', 'Itp')}, cen_mni);

    tgt = struct2table(row, 'AsArray', true);

    if options.Write
        csv = fullfile(dir_tab, sprintf('sub-%03d_%s%s_targeting.csv', id, medium, affix));
        writetable(tgt, csv);
        fprintf('Targeting: focus error %.1f mm (axial %+.1f, lateral %.1f); MNI target error %.1f mm -> %s\n', ...
            row.focus_error_mm, row.focus_error_axial_mm, row.focus_error_lateral_mm, ...
            row.target_error_mni_mm, csv);
    end
end

% =========================================================================
function [dir_nii, dir_tab] = output_dirs(parameters)
    out = '';
    if isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        out = parameters.io.dir_output;
    else
        out = get_output_dir(parameters);
    end
    dir_nii = fullfile(out, 'nii');
    if isfield(parameters.io, 'dir_nii_T1w') && ~isempty(parameters.io.dir_nii_T1w)
        dir_nii = parameters.io.dir_nii_T1w;
    end
    dir_tab = out;
    if isfield(parameters.io, 'dir_tabular') && ~isempty(parameters.io.dir_tabular)
        dir_tab = parameters.io.dir_tabular;
    end
end

function [cen, peak, vol_mm3, quantity] = focal_summary(mask_file, intensity_files)
% Centroid of the focal mask, peak of the first available intensity map in it.
    cen = nan(1, 3); peak = nan(1, 3); vol_mm3 = NaN; quantity = '';
    if ~isfile(mask_file); return; end
    hdr = niftiinfo(mask_file);
    m = niftiread(hdr) > 0;
    if ~any(m(:)); return; end
    [i, j, k] = ind2sub(size(m), find(m));
    ras = grid_to_ras([i j k], hdr);
    cen = mean(ras, 1);
    vol_mm3 = nnz(m) * voxel_volume(hdr);
    for f = 1:numel(intensity_files)
        if ~isfile(intensity_files{f}); continue; end
        ih = niftiinfo(intensity_files{f});
        if ~isequal(ih.ImageSize(1:3), hdr.ImageSize(1:3)); continue; end
        I = niftiread(ih);
        [~, imax] = max(double(I(m)));
        peak = ras(imax, :);
        [~, nm] = fileparts(intensity_files{f});
        quantity = regexprep(nm, '^.*_(\w+)\.nii$', '$1');
        return
    end
end

function row = roi_metrics(row, parameters, mni_mask_file, mni_int_files, cen_mni)
    names = {'roi_volume_mm3', 'roi_overlap_mm3', 'roi_focal_fraction', 'roi_coverage', ...
             'roi_mean_intensity', 'roi_max_intensity', 'roi_centroid_dist_mm'};
    for k = 1:numel(names); row.(names{k}) = NaN; end
    roi_file = '';
    if isfield(parameters, 'analysis') && isfield(parameters.analysis, 'target_roi_mni')
        roi_file = char(string(parameters.analysis.target_roi_mni));
    end
    if isempty(roi_file); return; end
    if ~isfile(roi_file)
        warning('targeting_analysis:roiMissing', 'analysis.target_roi_mni not found: %s', roi_file);
        return
    end
    rh = niftiinfo(roi_file);
    roi = niftiread(rh) > 0;
    if ~any(roi(:)); return; end
    [i, j, k] = ind2sub(size(roi), find(roi));
    roi_pts = grid_to_ras([i j k], rh);
    row.roi_volume_mm3 = nnz(roi) * voxel_volume(rh);
    row.roi_centroid_dist_mm = norm_or_nan(cen_mni - mean(roi_pts, 1));

    if ~isfile(mni_mask_file); return; end
    mh = niftiinfo(mni_mask_file);
    m = niftiread(mh) > 0;
    [i, j, k] = ind2sub(size(m), find(m));
    focal_pts = grid_to_ras([i j k], mh);
    in_roi = sample_volume_at_world(roi, rh, focal_pts) > 0;
    row.roi_overlap_mm3    = nnz(in_roi) * voxel_volume(mh);
    row.roi_focal_fraction = nnz(in_roi) / max(1, numel(in_roi));
    row.roi_coverage       = min(1, row.roi_overlap_mm3 / row.roi_volume_mm3);

    for f = 1:numel(mni_int_files)
        if ~isfile(mni_int_files{f}); continue; end
        ih = niftiinfo(mni_int_files{f});
        v = sample_volume_at_world(niftiread(ih), ih, roi_pts);
        v = v(isfinite(v));
        if ~isempty(v)
            row.roi_mean_intensity = mean(v);
            row.roi_max_intensity  = max(v);
        end
        return
    end
end

function [t, source] = resolve_target_mni(parameters)
    t = nan(1, 3); source = '';
    cands = {{'analysis', 'target_mni_mm'}, 'analysis.target_mni_mm'};
    mode = '';
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'mode')
        mode = char(string(parameters.placement.mode));
    end
    switch mode
        case 'heuristic', cands(end+1, :) = {{'placement', 'heuristic', 'mni_target_mm'}, 'placement.heuristic.mni_target_mm'};
        case 'plantus',   cands(end+1, :) = {{'placement', 'plantus', 'mni_target_mm'}, 'placement.plantus.mni_target_mm'};
        case 'mni',       cands(end+1, :) = {{'placement', 'mni', 'focus_pos_mm'}, 'placement.mni.focus_pos_mm'};
    end
    for c = 1:size(cands, 1)
        s = parameters; ok = true;
        for p = 1:numel(cands{c, 1})
            if ~isstruct(s) || ~isfield(s, cands{c, 1}{p}); ok = false; break; end
            s = s.(cands{c, 1}{p});
        end
        if ok && isnumeric(s) && numel(s) == 3 && all(isfinite(s))
            t = double(s(:)'); source = cands{c, 2};
            return
        end
    end
end

function v = voxel_volume(hdr)
    v = abs(det(hdr.Transform.T(1:3, 1:3)));
end

function hdr = header_of(file)
    hdr = [];
    if isfile(file); hdr = niftiinfo(file); end
end

function row = put_xyz(row, name, v)
    row.([name '_x']) = v(1);
    row.([name '_y']) = v(2);
    row.([name '_z']) = v(3);
end

function n = norm_or_nan(v)
    if all(isfinite(v)); n = norm(v); else, n = NaN; end
end

function v = first_or_nan(x)
    if isempty(x) || ~isnumeric(x); v = NaN; else, v = double(x(1)); end
end
