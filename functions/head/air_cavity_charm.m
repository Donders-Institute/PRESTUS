function [mask, note] = air_cavity_charm(parameters, planimg, grid_size)
% AIR_CAVITY_CHARM  Internal air from the SimNIBS charm atlas labelling, on the simulation grid
%
% SimNIBS charm labels internal air (sinuses, nasal cavity, airway, mastoid
% air cells) as 'Air-Internal' (label 501) in its atlas labelling,
% <m2m>/segmentation/labeling.nii.gz, but merges it into the background (0)
% of final_tissues.nii.gz. This function reads the labelling, selects
% headmodel.air_cavities.simnibs_labels and maps the mask onto the
% simulation grid with the transform of the head model (planimg.transf),
% nearest neighbour. The labelling shares the grid of final_tissues.nii.gz.
%
% Use as:
%   [mask, note] = air_cavity_charm(parameters, planimg, grid_size)
%
% Input:
%   parameters - PRESTUS config (path.seg, subject_id, headmodel.air_cavities)
%   planimg    - planning-image struct from GRID_STAGE (transf: T1 voxel -> grid voxel)
%   grid_size  - size of the segmentation on the simulation grid
%
% Output:
%   mask - logical volume of grid_size, or [] when the source is not
%          'simnibs' or the mask cannot be built
%   note - why the mask is absent ('' when built)
%
% Non-layered media return [] silently. A missing labelling warns
% (prestus:airCavitiesUnavailable): charm writes it with every
% segmentation, so it is absent only when removed afterwards.
%
% See also: AIR_CAVITY_MASK, PREPROC_HEAD, AFFINE_RESAMPLE_3D

    mask = [];
    note = '';
    cfg = get_field(parameters, {'headmodel', 'air_cavities'}, struct());
    if ~strcmp(char(get_field(cfg, {'source'}, 'none')), 'simnibs')
        return
    end
    if ~strcmp(get_field(parameters, {'simulation', 'medium'}, ''), 'layered')
        note = 'air-cavity check needs simulation.medium: layered';     % water, phantom: nothing to check
        return
    end
    if numel(grid_size) ~= 3 || ~isfield(planimg, 'transf') || isempty(planimg.transf)
        note = 'needs a 3D head-model grid (not axisymmetric or 2D)';
        fprintf('[air_cavity_charm] Skipped: %s.\n', note);
        return
    end
    seg_dir = fullfile(char(get_field(parameters, {'path', 'seg'}, '')), ...
        sprintf('m2m_sub-%03d', get_field(parameters, {'subject_id'}, 0)));
    f = fullfile(seg_dir, 'segmentation', 'labeling.nii.gz');
    if ~isfile(f)
        % charm writes it with every segmentation; it is missing only when removed afterwards
        note = sprintf('charm labelling not found (%s)', f);
        warning('prestus:airCavitiesUnavailable', ...
            ['Air-cavity check skipped: %s. Re-run the SimNIBS segmentation (io.overwrite_simnibs: 1), ' ...
             'or set headmodel.air_cavities.source to ''hu'' or ''none''.'], note);
        return
    end

    labels = double(get_field(cfg, {'simnibs_labels'}, 501));
    img = niftiread(f);
    hdr = niftiinfo(f);
    [img, ~] = ensure_ras_plus(img, hdr);
    air = uint8(ismember(img, labels));
    clear img
    mask = affine_resample_3d(air, planimg.transf, grid_size(:)', 'nearest', 0) > 0;
    fprintf('[air_cavity_charm] %s labels %s: %d voxels in the labelling, %d on the simulation grid.\n', ...
        f, mat2str(labels), nnz(air), nnz(mask));
end
