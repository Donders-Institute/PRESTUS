function out = seg_qc(action, varargin)
% SEG_QC  Record head-model postprocessing steps for the segmentation report
%
% Each step of turning the SimNIBS segmentation into the simulation medium
% (tissue mapping and smoothing, skull hole filling, eyes to water, CSF-based
% crop, grid crop) adds one record while a session is open: tissue voxel
% counts before and after, the voxels moved between tissues, and three
% orthogonal slices of both states through a reference point (the focus).
% Whole volumes are not kept. GENERATE_SEGMENTATION_REPORT renders the
% saved records.
%
% A session is opened by PREPROC_HEAD when report.segmentation is on
% (default) and the medium is layered; without an open session 'add' and
% 'add_mapping' do nothing, so the step functions can be called on their own.
%
% Use as:
%   seg_qc('begin', parameters)               % open a session
%   tf = seg_qc('active')
%   seg_qc('add', id, before, after, ref_pos, meta)
%   seg_qc('add_mapping', parameters, segmentation, medium_masks, ref_pos, grid_note)
%   steps = seg_qc('steps')                   % records so far
%   seg_qc('import', steps)                   % prepend earlier records
%   seg_qc('save', parameters)                % write SEG_QC('file', parameters)
%   seg_qc('end')                             % close the session
%   file = seg_qc('file', parameters)
%   note = seg_qc('grid_note', parameters, where)  % 'placement' | 'headmodel' | 'cropped'
%
% Input ('add'):
%   id       - step id ('segmentation', 'smoothing', 'skull_fill', 'eyes',
%              'csf_crop', 'grid_crop')
%   before   - label volume before the step ([] for a state without a before)
%   after    - label volume after the step
%   ref_pos  - [1x3] voxel through which the slices are taken
%   meta     - struct: title, description, grid, label_space ('medium' |
%              'charm'), and optional notes
%
% The record file is <io.dir_cache>/sub-XXX_<medium>_segqc<affix>.mat next
% to the cropped head model it describes (affix: io.preproc_affix when set,
% else io.output_affix), so a run that reuses the head model reuses it too.
%
% See also: GENERATE_SEGMENTATION_REPORT, PREPROC_HEAD, HEAD_SMOOTH_AND_CROP

persistent S is_active info

    out = [];
    switch action
        case 'begin'
            p = varargin{1};
            S = empty_steps();
            info = session_info(p);
            is_active = report_enabled(p);

        case 'end'
            is_active = false;
            S = empty_steps();

        case 'active'
            out = ~isempty(is_active) && is_active;

        case 'add'
            if isempty(is_active) || ~is_active, return; end
            [id, before, after, ref_pos, meta] = varargin{1:5};
            S(end+1) = make_step(id, before, after, ref_pos, meta);

        case 'add_mapping'
            if isempty(is_active) || ~is_active, return; end
            [p, segmentation, medium_masks, ref_pos, grid_note] = varargin{1:5};
            S = [S, mapping_steps(p, segmentation, medium_masks, ref_pos, grid_note)];

        case 'steps'
            if isempty(S), S = empty_steps(); end
            out = S;

        case 'import'
            if isempty(is_active) || ~is_active, return; end
            prev = varargin{1};
            if isempty(prev), return; end
            S = [prev(:)', S];

        case 'save'
            if isempty(is_active) || ~is_active, return; end
            p = varargin{1};
            file = seg_qc('file', p);
            steps = S; %#ok<NASGU>
            qc_info = info; %#ok<NASGU>
            try
                save(file, 'steps', 'qc_info', '-v7');
                fprintf('[seg_qc] head-model QC records saved: %s\n', file);
            catch ME
                warning('PRESTUS:segqc:save', 'Could not save the head-model QC records (%s).', ME.message);
            end

        case 'grid_note'
            % seg_qc('grid_note', parameters, where): where = 'placement' | 'headmodel' | 'cropped'
            p = varargin{1};
            mode = 'focal_axis';
            if isfield(p, 'grid') && isfield(p.grid, 'mode'); mode = p.grid.mode; end
            switch varargin{2}
                case 'headmodel'
                    out = 'shared head model (RAS+, unrotated; slices through the volume centre)';
                case 'cropped'
                    out = 'simulation grid (cropped)';
                otherwise
                    if strcmp(mode, 'ras_plus')
                        out = 'RAS+ grid at the simulation resolution';
                    else
                        out = 'grid rotated to the beam axis (z = beam axis)';
                    end
            end

        case 'file'
            p = varargin{1};
            affix = '';
            if isfield(p.io, 'preproc_affix')
                affix = p.io.preproc_affix;
            elseif isfield(p.io, 'output_affix')
                affix = p.io.output_affix;
            end
            out = fullfile(p.io.dir_cache, sprintf('sub-%03d_%s_segqc%s.mat', ...
                p.subject_id, p.simulation.medium, affix));

        otherwise
            error('PRESTUS:segqc:action', 'Unknown seg_qc action ''%s''.', action);
    end

end

function S = mapping_steps(p, segmentation, medium_masks, ref_pos, grid_note)
% Records of the SimNIBS segmentation and of the tissue mapping with smoothing
    S = empty_steps();
    hm = p.headmodel;
    S(end+1) = make_step('segmentation', [], segmentation, ref_pos, struct( ...
        'title', 'SimNIBS segmentation', ...
        'description', 'charm tissue labels resampled to the simulation resolution.', ...
        'grid', grid_note, 'label_space', 'charm'));

    % Same mapping without smoothing: the "before" of the smoothing step
    p_raw = p;
    p_raw.headmodel.smooth_method = 'off';
    raw = preproc_medium_mask(segmentation, p_raw);
    if strcmp(hm.smooth_method, 'off')
        desc = 'Tissue labels mapped to medium layers (layers.*) without smoothing (headmodel.smooth_method: off).';
    else
        desc = sprintf(['Tissue labels mapped to medium layers (layers.*); each layer mask smoothed ' ...
            '(%s, FWHM %g mm) and thresholded (skull %g, other tissues %g). Later layers overwrite ' ...
            'earlier ones. Before: the same mapping without smoothing.'], hm.smooth_method, ...
            hm.smooth_fwhm_mm, hm.smooth_threshold_skull, hm.smooth_threshold_other);
    end
    S(end+1) = make_step('smoothing', raw, medium_masks, ref_pos, struct( ...
        'title', 'Tissue mapping and smoothing', 'description', desc, ...
        'grid', grid_note, 'label_space', 'medium'));
end

% ── Records ─────────────────────────────────────────────────────────────
function S = empty_steps()
    S = struct('id', {}, 'title', {}, 'description', {}, 'grid', {}, 'label_space', {}, ...
        'notes', {}, 'size', {}, 'ref_pos', {}, 'counts_before', {}, 'counts_after', {}, ...
        'transitions', {}, 'n_changed', {}, 'slices_before', {}, 'slices_after', {});
end

function st = make_step(id, before, after, ref_pos, meta)
    sz = size(after, 1:3);
    ref = min(max(round(double(ref_pos(:)')), 1), sz);
    has_before = ~isempty(before) && isequal(size(before, 1:3), sz);

    K = double(max(after(:)));
    if has_before, K = max(K, double(max(before(:)))); end
    K = max(0, K);

    counts_after = label_counts(after, K);
    counts_before = [];
    transitions = [];
    n_changed = 0;
    slices_before = {};
    if has_before
        counts_before = label_counts(before, K);
        ch = before ~= after;
        n_changed = nnz(ch);
        transitions = zeros(K+1, K+1);
        if n_changed > 0
            transitions = accumarray([double(before(ch)) + 1, double(after(ch)) + 1], 1, [K+1, K+1]);
        end
        clear ch
        slices_before = take_slices(before, ref);
    end

    notes = '';
    if isfield(meta, 'notes'); notes = meta.notes; end
    st = struct('id', id, 'title', meta.title, 'description', meta.description, ...
        'grid', meta.grid, 'label_space', meta.label_space, 'notes', notes, ...
        'size', sz, 'ref_pos', ref, 'counts_before', counts_before, 'counts_after', counts_after, ...
        'transitions', transitions, 'n_changed', n_changed, ...
        'slices_before', {slices_before}, 'slices_after', {take_slices(after, ref)});
end

function c = label_counts(vol, K)
    c = zeros(1, K+1);
    for k = 0:K
        c(k+1) = nnz(vol == k);
    end
end

function sl = take_slices(vol, ref)
% Planes through ref: fixed dim 1 (j,k), fixed dim 2 (i,k), fixed dim 3 (i,j)
    v = uint8(min(max(vol, 0), 255));
    sl = { squeeze(v(ref(1), :, :)), squeeze(v(:, ref(2), :)), v(:, :, ref(3)) };
end

% ── Session settings ───────────────────────────────────────────────────
function tf = report_enabled(p)
    tf = isfield(p, 'simulation') && isfield(p.simulation, 'medium') && ...
        strcmp(p.simulation.medium, 'layered');
    if tf && isfield(p, 'report') && isfield(p.report, 'segmentation') && ~isempty(p.report.segmentation)
        tf = isequal(double(p.report.segmentation), 1) || isequal(p.report.segmentation, true);
    end
end

function info = session_info(p)
    info = struct();
    info.labels_medium = {};
    if isfield(p, 'medium_properties'); info.labels_medium = fieldnames(p.medium_properties)'; end
    info.resolution_mm = p.grid.resolution_mm;
    info.headmodel = struct();
    if isfield(p, 'headmodel'); info.headmodel = p.headmodel; end
    info.layers = struct();
    if isfield(p, 'layers'); info.layers = p.layers; end
    info.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end
