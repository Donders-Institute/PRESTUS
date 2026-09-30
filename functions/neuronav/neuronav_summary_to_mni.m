function T = neuronav_summary_to_mni(T, parameters, varargin)
% NEURONAV_SUMMARY_TO_MNI  Add MNI coordinates to a neuronav ingest summary.
%
% Maps the Localite pose, transducer and focus coordinates of each summary
% row (NEURONAV_INGEST_BATCH, planning-image RAS mm) to MNI RAS mm with the
% subject's SimNIBS nonlinear field, as in NEURONAV_CONVERT_NATIVE_TO_MNI:
% planning RAS -> planning voxel -> segmentation RAS (the planning T1 and
% final_tissues share the voxel grid) -> Conform2MNI_nonl. Unlike that
% function, the result is not rounded or clamped to the MNI image.
%
% Use as:
%   T = neuronav_summary_to_mni(T, parameters)
%   T = neuronav_summary_to_mni(T, parameters, 'progress', fcn)
%
% Input:
%   T          - summary table (or path of neuronav_ingest_summary.csv)
%   parameters - config struct; uses path.seg (m2m_<sub>/final_tissues.nii.gz
%                and m2m_<sub>/toMNI/Conform2MNI_nonl.nii.gz) and the
%                planning T1 settings of NEURONAV_FIND_PLANNING_T1
%
% Name-value options:
%   'progress' - @(k, n, sub, status) called after each subject
%
% Output:
%   T - input table plus pose_mni_x/y/z, trans_mni_x/y/z, focus_mni_x/y/z
%       (NaN where the input is NaN or the subject cannot be mapped) and
%       mni_status ('' = mapped, otherwise the reason)
%
% See also: NEURONAV_INGEST_BATCH, NEURONAV_CONVERT_NATIVE_TO_MNI,
%           NEURONAV_APPLY_DEFORMATION

    ip = inputParser;
    ip.addParameter('progress', []);
    ip.parse(varargin{:});
    opt = ip.Results;

    if ischar(T) || isstring(T)
        opts = detectImportOptions(char(T));
        opts = setvartype(opts, intersect({'sub', 'ses', 'target', 'rule', 'label', ...
            'label_source', 'markertype', 'source_file', 'xml', 'status'}, opts.VariableNames), 'char');
        T = readtable(char(T), opts);
    end

    kinds = {'pose', 'pose_%s'; 'trans', 'trans_ras_%s'; 'focus', 'focus_ras_%s'};
    n = height(T);
    for k = 1:size(kinds, 1)
        for ax = 'xyz'
            T.(sprintf('%s_mni_%s', kinds{k, 1}, ax)) = nan(n, 1);
        end
    end
    T.mni_status = repmat({''}, n, 1);

    seg_root = '';
    if isfield(parameters, 'path') && isfield(parameters.path, 'seg')
        seg_root = char(parameters.path.seg);
    end

    subs = unique(T.sub, 'stable');
    for s = 1:numel(subs)
        sub_id = subs{s};
        rows = strcmp(T.sub, sub_id);
        status = '';
        try
            if isempty(seg_root), error('path.seg is not set.'); end
            X = load_subject(sub_id, parameters, seg_root);
            for k = 1:size(kinds, 1)
                src = cellfun(@(a) sprintf(kinds{k, 2}, a), {'x', 'y', 'z'}, 'UniformOutput', false);
                if ~all(ismember(src, T.Properties.VariableNames)), continue; end
                xyz = [T.(src{1})(rows), T.(src{2})(rows), T.(src{3})(rows)];
                mni = nan(size(xyz));
                ok = all(isfinite(xyz), 2);
                if any(ok), mni(ok, :) = planning_to_mni(xyz(ok, :), X); end
                axes_ = 'xyz';
                for a = 1:3
                    T.(sprintf('%s_mni_%s', kinds{k, 1}, axes_(a)))(rows) = mni(:, a);
                end
            end
        catch ME
            status = ME.message;
        end
        T.mni_status(rows) = {status};
        if ~isempty(opt.progress), opt.progress(s, numel(subs), sub_id, status); end
    end
end

function X = load_subject(sub_id, parameters, seg_root)
    m2m = fullfile(seg_root, sprintf('m2m_%s', sub_id));
    seg_file  = fullfile(m2m, 'final_tissues.nii.gz');
    warp_file = fullfile(m2m, 'toMNI', 'Conform2MNI_nonl.nii.gz');
    if ~isfile(seg_file), error('No %s.', seg_file); end
    if ~isfile(warp_file), error('No %s.', warp_file); end
    X.plan = niftiinfo(neuronav_find_planning_t1(sub_id, parameters));
    X.seg  = niftiinfo(seg_file);
    if ~isequal(X.plan.ImageSize(1:3), X.seg.ImageSize(1:3))
        % The planning -> segmentation step reinterprets planning voxels
        % under the segmentation header, so the grids must be identical.
        error('Planning T1 grid [%s] differs from final_tissues grid [%s].', ...
            num2str(X.plan.ImageSize(1:3)), num2str(X.seg.ImageSize(1:3)));
    end
    X.warp   = niftiread(warp_file);
    warp_info = niftiinfo(warp_file);
    X.warp_T = warp_info.Transform.T;
    X.nearest_idx = neuronav_build_deformation_cache(X.warp);
end

function mni = planning_to_mni(xyz, X)
    % Planning voxel kept continuous (see NEURONAV_CONVERT_NATIVE_TO_MNI),
    % then the forward segmentation affine as in TRANSFORM_COORDINATES
    % ('grid' -> 'ras_plus'), vectorised over rows.
    vox = [xyz, ones(size(xyz, 1), 1)] / X.plan.Transform.T;
    ras_seg = transformPointsForward(X.seg.Transform, vox(:, 1:3));
    mni = neuronav_apply_deformation(ras_seg, X.warp, X.warp_T, X.nearest_idx);
end
