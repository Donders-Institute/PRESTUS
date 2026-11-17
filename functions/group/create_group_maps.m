function maps = create_group_maps(runs, out_dir, options)
% CREATE_GROUP_MAPS  Group maps in MNI space from the per-run MNI NIfTIs, per condition
%
% For each condition (all runs when there are none) combines the MNI-space
% outputs that the pipeline writes with io.save_MNI:
%
%   group_<cond>_focal_overlap_3dB.nii.gz  fraction of runs whose -3 dB focal
%                                          volume covers each voxel (0..1)
%   group_<cond>_<Q>_mean.nii.gz / _sd     mean / SD intensity, Q = Ipa (Itp
%                                          when a run has no Ipa map)
%   group_<cond>_heatrise_mean.nii.gz      mean temperature rise (when all
%                                          runs have heatrise maps)
%
% in <out_dir>/maps/, plus one figure per condition in <out_dir>/plots/:
% three orthogonal MNI slices through the voxel of highest overlap over the
% mean MNI-space T1 of the subjects (GET_T1_MNI_FILE, when present), with
% the ROI outline, the per-run focal centroids, and the targets from the
% targeting sidecars. Runs whose MNI grid differs from the first run's are
% skipped with a warning.
%
% Use as:
%   maps = create_group_maps(runs, out_dir)
%   maps = create_group_maps(runs, out_dir, 'SegDir', p.path.seg, 'RoiMni', roi)
%
% Input:
%   runs    - table with subject_id, condition, affix, medium, sim_dir
%             (STUDY_RUN_TABLE, included rows only)
%   out_dir - group output folder
%   options - name-value pairs:
%     'SegDir' path.seg, locates the background T1 (GET_T1_MNI_FILE); '' = none
%     'RoiMni' MNI ROI mask (NIfTI) to outline; '' = none
%     'Plot'   write the figures (default true)
%
% Output:
%   maps - struct array, one per condition: condition, tag, n (runs used),
%          n_missing, subjects, quantity, files (struct of NIfTI paths),
%          figure (PNG path or '')
%
% See also: STUDY_RUN_TABLE, STUDY_RUN_FILES, TARGETING_ANALYSIS,
%           GENERATE_GROUP_REPORT, PRESTUS_GROUP_REPORT_START

arguments
    runs table
    out_dir {mustBeTextScalar}
    options.SegDir {mustBeTextScalar} = ''
    options.RoiMni {mustBeTextScalar} = ''
    options.Plot (1,1) logical = true
end

    maps = struct('condition', {}, 'tag', {}, 'n', {}, 'n_missing', {}, 'subjects', {}, ...
                  'quantity', {}, 'files', {}, 'figure', {});
    if isempty(runs); return; end
    map_dir  = fullfile(out_dir, 'maps');
    plot_dir = fullfile(out_dir, 'plots');

    conds = unique(string(runs.condition), 'stable');
    for c = 1:numel(conds)
        sel = runs(string(runs.condition) == conds(c), :);
        m = one_condition(sel, conds(c), map_dir, plot_dir, options);
        if ~isempty(m); maps(end+1) = m; end %#ok<AGROW>
    end
end

% =========================================================================
function m = one_condition(sel, cond, map_dir, plot_dir, options)
    m = [];
    tag = char(cond);
    if isempty(tag); tag = 'all'; end
    tag = regexprep(tag, '[^a-zA-Z0-9_-]', '_');

    files = cell(height(sel), 1);
    for r = 1:height(sel)
        files{r} = study_run_files(sel.subject_id(r), sel.medium(r), sel.affix(r), sel.sim_dir(r));
    end
    has_mask = cellfun(@(f) isfile(f.mni('focal_mask_3dB')), files);
    if ~any(has_mask)
        fprintf('create_group_maps: no MNI focal masks for condition "%s" — skipped.\n', tag);
        return
    end
    files = files(has_mask);
    subjects = sel.subject_id(has_mask);
    n_missing = nnz(~has_mask);

    quantity = '';
    if all(cellfun(@(f) isfile(f.mni('Ipa')), files))
        quantity = 'Ipa';
    elseif all(cellfun(@(f) isfile(f.mni('Itp')), files))
        quantity = 'Itp';
    end
    has_heat = all(cellfun(@(f) isfile(f.mni('heatrise')), files));

    ref = niftiinfo(files{1}.mni('focal_mask_3dB'));
    sz = ref.ImageSize(1:3);
    count = zeros(sz, 'single');
    isum = zeros(sz); isq = zeros(sz); hsum = zeros(sz);
    used = false(numel(files), 1);
    for r = 1:numel(files)
        h = niftiinfo(files{r}.mni('focal_mask_3dB'));
        if ~same_grid(h, ref)
            warning('create_group_maps:gridMismatch', ...
                'sub-%03d: MNI grid differs from sub-%03d; skipped for condition "%s".', ...
                subjects(r), subjects(1), tag);
            continue
        end
        count = count + single(niftiread(h) > 0);
        if ~isempty(quantity)
            v = double(niftiread(files{r}.mni(quantity)));
            v(~isfinite(v)) = 0;
            isum = isum + v; isq = isq + v.^2;
        end
        if has_heat
            v = double(niftiread(files{r}.mni('heatrise')));
            v(~isfinite(v)) = 0;
            hsum = hsum + v;
        end
        used(r) = true;
    end
    n = nnz(used);
    if n == 0; return; end
    files = files(used); subjects = subjects(used);

    if ~isfolder(map_dir); mkdir(map_dir); end
    out = struct();
    out.overlap = write_map(count / n, ref, fullfile(map_dir, sprintf('group_%s_focal_overlap_3dB', tag)));
    if ~isempty(quantity)
        mu = isum / n;
        sd = sqrt(max(0, isq / n - mu.^2) * n / max(1, n - 1));
        out.intensity_mean = write_map(mu, ref, fullfile(map_dir, sprintf('group_%s_%s_mean', tag, quantity)));
        out.intensity_sd   = write_map(sd, ref, fullfile(map_dir, sprintf('group_%s_%s_sd', tag, quantity)));
    end
    if has_heat
        out.heatrise_mean = write_map(hsum / n, ref, fullfile(map_dir, sprintf('group_%s_heatrise_mean', tag)));
    end

    fig = '';
    if options.Plot
        try
            bg = mean_background(subjects, options.SegDir, ref);
            [cen, tgt] = sidecar_points(files);
            roi = [];
            if ~isempty(options.RoiMni) && isfile(options.RoiMni)
                roi = struct('hdr', niftiinfo(options.RoiMni));
                roi.vol = niftiread(roi.hdr) > 0;
            end
            if ~isfolder(plot_dir); mkdir(plot_dir); end
            fig = fullfile(plot_dir, sprintf('group_%s_focal_overlap.png', tag));
            title_str = sprintf('Focal overlap (-3 dB), %s: N = %d', tag, n);
            render_overlap(count / n, ref, bg, roi, cen, tgt, title_str, fig);
        catch ME
            warning('create_group_maps:plotFailed', 'Overlap figure failed (%s): %s', tag, ME.message);
            fig = '';
        end
    end

    m = struct('condition', char(cond), 'tag', tag, 'n', n, 'n_missing', n_missing, ...
        'subjects', subjects(:)', 'quantity', quantity, 'files', out, 'figure', fig);
    fprintf('create_group_maps: %s — %d run(s) -> %s\n', tag, n, map_dir);
end

function tf = same_grid(a, b)
    tf = isequal(a.ImageSize(1:3), b.ImageSize(1:3)) && ...
         max(abs(a.Transform.T(:) - b.Transform.T(:))) < 1e-3;
end

function file = write_map(vol, ref, stem)
    hdr = ref;
    hdr.Datatype = 'single';
    hdr.BitsPerPixel = 32;
    if isfield(hdr, 'MultiplicativeScaling'); hdr.MultiplicativeScaling = 1; end
    if isfield(hdr, 'AdditiveOffset');        hdr.AdditiveOffset = 0; end
    niftiwrite(single(vol), stem, hdr, 'Compressed', true);
    file = [stem '.nii.gz'];
end

function bg = mean_background(subjects, seg_dir, ref)
% Mean MNI-space T1 of the subjects on the reference grid ([] when none).
    bg = [];
    if isempty(seg_dir); return; end
    acc = []; k = 0;
    for s = subjects(:)'
        [f, f_m2m] = get_t1_mni_file(seg_dir, s);
        if ~isfile(f); f = f_m2m; end                   % earlier outputs (read only)
        if ~isfile(f); continue; end
        h = niftiinfo(f);
        if ~same_grid(h, ref); continue; end
        v = double(niftiread(h));
        v(~isfinite(v)) = 0;
        vs = sort(v(v > 0));
        if isempty(vs); continue; end
        hi = vs(max(1, round(0.99 * numel(vs))));      % 99th percentile, no toolbox
        if hi <= 0; continue; end
        v = min(v / hi, 1);
        if isempty(acc); acc = v; else, acc = acc + v; end
        k = k + 1;
    end
    if k > 0; bg = acc / k; end
end

function [cen, tgt] = sidecar_points(files)
% Per-run MNI focal centroids and MNI targets from the targeting sidecars.
    cen = zeros(0, 3); tgt = zeros(0, 3);
    for r = 1:numel(files)
        if ~isfile(files{r}.targeting_csv); continue; end
        t = readtable(files{r}.targeting_csv, 'VariableNamingRule', 'preserve');
        c = xyz(t, 'focal_centroid_mni');
        g = xyz(t, 'target_mni');
        if all(isfinite(c)); cen(end+1, :) = c; end %#ok<AGROW>
        if all(isfinite(g)); tgt(end+1, :) = g; end %#ok<AGROW>
    end
    tgt = unique(round(tgt, 2), 'rows');
end

function v = xyz(t, name)
    v = nan(1, 3); ax = 'xyz';
    for a = 1:3
        col = sprintf('%s_%s', name, ax(a));
        if ismember(col, t.Properties.VariableNames) && isnumeric(t.(col))
            v(a) = t.(col)(end);
        end
    end
end

% -------------------------------------------------------------------------
function render_overlap(ov, hdr, bg, roi, cen, tgt, title_str, png)
% Three orthogonal slices through the overlap peak, in MNI mm.
    T = hdr.Transform.T;
    sz = size(ov);
    if max(abs(T(1:3, 1:3) - diag(diag(T(1:3, 1:3)))), [], 'all') > 1e-6
        error('create_group_maps:obliqueGrid', 'MNI grid is not axis-aligned; figure skipped.');
    end
    w = cell(1, 3);
    for a = 1:3
        w{a} = T(4, a) + T(a, a) * (0:sz(a) - 1);
    end
    peak = max(ov(:));
    [i, j, k] = ind2sub(sz, find(ov >= peak - 1e-6));
    c = round(mean([i j k], 1));

    cmap = hot(256);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 430]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    planes = {'Sagittal (x = %.0f mm)', 'Coronal (y = %.0f mm)', 'Axial (z = %.0f mm)'};
    fixed = [1 2 3];
    for p = 1:3
        ax = nexttile(tl);
        dims = setdiff(1:3, fixed(p));             % in-plane axes (horizontal, vertical)
        sl_ov = take_slice(ov, fixed(p), c(fixed(p)));
        if isempty(bg); sl_bg = zeros(size(sl_ov)); else, sl_bg = take_slice(bg, fixed(p), c(fixed(p))); end
        rgb_bg = repmat(sl_bg', 1, 1, 3);
        image(ax, w{dims(1)}, w{dims(2)}, rgb_bg);
        hold(ax, 'on');
        idx = 1 + round(255 * min(1, sl_ov' / max(peak, eps)));
        rgb = reshape(cmap(idx, :), [size(idx) 3]);
        image(ax, w{dims(1)}, w{dims(2)}, rgb, 'AlphaData', 0.85 * (sl_ov' > 0));
        if ~isempty(roi)
            [g1, g2] = ndgrid(1:sz(dims(1)), 1:sz(dims(2)));
            vox = zeros(numel(g1), 3);
            vox(:, dims(1)) = g1(:); vox(:, dims(2)) = g2(:); vox(:, fixed(p)) = c(fixed(p));
            rv = sample_volume_at_world(roi.vol, roi.hdr, grid_to_ras(vox, hdr));
            rv = reshape(rv > 0, size(g1));
            if any(rv(:))
                contour(ax, w{dims(1)}, w{dims(2)}, double(rv'), [0.5 0.5], 'c', 'LineWidth', 1.2);
            end
        end
        if ~isempty(cen)
            plot(ax, cen(:, dims(1)), cen(:, dims(2)), 'o', 'MarkerSize', 4, ...
                'MarkerFaceColor', [0.2 0.6 1], 'MarkerEdgeColor', 'w');
        end
        if ~isempty(tgt)
            plot(ax, tgt(:, dims(1)), tgt(:, dims(2)), '+', 'MarkerSize', 12, ...
                'Color', [0 1 0.4], 'LineWidth', 2);
        end
        set(ax, 'YDir', 'normal', 'Color', 'k', 'FontSize', 8);
        axis(ax, 'image');
        % zoom to +-50 mm around the overlap peak (within the volume)
        for q = 1:2
            ctr = w{dims(q)}(c(dims(q)));
            lim = sort([max(min(w{dims(q)}), ctr - 50), min(max(w{dims(q)}), ctr + 50)]);
            if q == 1; xlim(ax, lim); else, ylim(ax, lim); end
        end
        labels = 'xyz';
        xlabel(ax, sprintf('%s (mm)', labels(dims(1)))); ylabel(ax, sprintf('%s (mm)', labels(dims(2))));
        title(ax, sprintf(planes{p}, w{fixed(p)}(c(fixed(p)))), 'FontWeight', 'normal');
        colormap(ax, cmap); clim(ax, [0 peak]);
    end
    cb = colorbar(ax); cb.Label.String = 'fraction of runs';
    title(tl, sprintf('%s   (dots: focal centroids, +: targets)', title_str), 'Interpreter', 'none');
    if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
        exportgraphics(fig, png, 'Resolution', 150);
    else
        print(fig, png, '-dpng', '-r150');
    end
    close(fig);
end

function s = take_slice(v, dim, idx)
    switch dim
        case 1, s = squeeze(v(idx, :, :));
        case 2, s = squeeze(v(:, idx, :));
        otherwise, s = squeeze(v(:, :, idx));
    end
end
