function plot_intensity_t1_overlay(intensity_t1, planimg, parameters, results_acoustic, highlighted_pos)
% PLOT_INTENSITY_T1_OVERLAY  Overlay back-transformed intensity on T1 image
%
% Plots the intensity map (already in T1 space) over the T1 anatomical image
% at the coronal slice through the focus, with correct orientation. Saves one
% PNG per transducer into parameters.io.dir_img.
%
% Use as:
%   plot_intensity_t1_overlay(intensity_t1, planimg, parameters, results_acoustic, highlighted_pos)
%
% Input:
%   intensity_t1     - (:,:,:) single, Isppa map back-transformed to T1 space
%   planimg          - struct with fields: t1_image_orig, inv_transf, t1_header
%   parameters       - PRESTUS config struct
%   results_acoustic - struct with fields: Isppa (and optionally Isppa_brain)
%   highlighted_pos  - [1x3] grid-space peak intensity position [voxels]
%
% See also: NIFTI_ACOUSTIC, PLOT_OVERLAY

n_trans = numel(parameters.transducer);
sync_mode = n_trans > 1 && ~is_async_mode(parameters);

if sync_mode
    max_plots = 1;
elseif n_trans > 2
    warn('More than two transducers: intensity-over-T1 plots will be created only for the first 2 transducers');
    max_plots = 2;
else
    max_plots = n_trans;
end

if isfield(results_acoustic, 'Isppa_brain') && ~isempty(results_acoustic.Isppa_brain) && ~isnan(results_acoustic.Isppa_brain)
    max_val = results_acoustic.Isppa_brain;
else
    max_val = results_acoustic.Isppa;
end

% For sync mode, collect all transducer/focus positions as stacked rows so
% plot_overlay draws every bowl outline, and use trans-01 focus as the
% reference coronal slice position.
if sync_mode
    all_tpos = zeros(n_trans, 3);
    all_fpos = zeros(n_trans, 3);
    all_tpos_t1 = zeros(n_trans, 3);
    all_fpos_t1 = zeros(n_trans, 3);
    source_labels = zeros(size(planimg.t1_image_orig));
    for ti = 1:n_trans
        all_tpos(ti,:) = parameters.transducer(ti).trans_pos;
        all_fpos(ti,:) = parameters.transducer(ti).focus_pos;
        btc_i = round(affine_apply_pts( ...
            [all_tpos(ti,:); all_fpos(ti,:); highlighted_pos], planimg.inv_transf));
        all_tpos_t1(ti,:) = btc_i(1,:);
        all_fpos_t1(ti,:) = btc_i(2,:);
        try
            [~, lbl_i] = transducer_setup( ...
                parameters.transducer(ti), btc_i(1,:), btc_i(2,:), ...
                size(planimg.t1_image_orig), planimg.t1_header.PixelDimensions(1), parameters);
            source_labels = source_labels | lbl_i;
        catch ME
            warning('plot_intensity_t1_overlay: transducer_setup failed for transducer %d (%s) — bowl outline will be omitted.', ti, ME.message);
        end
    end
    % Use trans-01 focus as the reference coronal slice
    ref_btc = round(affine_apply_pts( ...
        [all_tpos(1,:); all_fpos(1,:); highlighted_pos], planimg.inv_transf));

    if ~isempty(planimg.t1_header) && isfield(planimg.t1_header, 'Transform')
        [do_t, flip_r, flip_c] = nifti_slice_orientation(planimg.t1_header, 'y');
    else
        [do_t, flip_r, flip_c] = deal(false, false, false);
    end

    try
        [~,~,~,~,ax_ov,~,~,h] = plot_overlay( ...
            intensity_t1, planimg.t1_image_orig, source_labels, parameters, ...
            {'y', ref_btc(2,2)}, ...
            all_tpos_t1, all_fpos_t1, ref_btc(3,:), ...
            'show_rectangles', 0, ...
            'grid_step', planimg.t1_header.PixelDimensions(1), ...
            'overlay_color_range', [0 max_val], ...
            'overlay_threshold_low', 0, ...
            'overlay_threshold_high', max_val, ...
            'rotation', 0, ...
            'do_transpose', do_t, 'flip_rows', flip_r, 'flip_cols', flip_c);
        draw_fov_outline(ax_ov, parameters, planimg, ref_btc(2,2), do_t, flip_r, flip_c);

        output_plot_filename = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s_intensity_t1w%s.png', ...
            parameters.subject_id, parameters.simulation.medium, ...
            parameters.io.output_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    catch ME_plot
        warning('plot_intensity_t1_overlay: plot_overlay failed for sync combined plot (%s) — intensity overlay plot will be skipped.', ME_plot.message);
    end

else
    for ti = 1:max_plots
        tpos_sim = parameters.transducer(ti).trans_pos;
        fpos_sim = parameters.transducer(ti).focus_pos;

        backtransf_coordinates = round(affine_apply_pts( ...
            [tpos_sim; fpos_sim; highlighted_pos], planimg.inv_transf));

        try
            [~, source_labels] = transducer_setup( ...
                parameters.transducer(ti), ...
                backtransf_coordinates(1,:), ...
                backtransf_coordinates(2,:), ...
                size(planimg.t1_image_orig), ...
                planimg.t1_header.PixelDimensions(1), parameters);
        catch ME
            warning('plot_intensity_t1_overlay: transducer_setup failed for transducer %d (%s) — bowl outline will be omitted from plot.', ti, ME.message);
            source_labels = zeros(size(planimg.t1_image_orig));
        end

        if ~isempty(planimg.t1_header) && isfield(planimg.t1_header, 'Transform')
            [do_t, flip_r, flip_c] = nifti_slice_orientation(planimg.t1_header, 'y');
        else
            [do_t, flip_r, flip_c] = deal(false, false, false);
        end

        try
            [~,~,~,~,ax_ov,~,~,h] = plot_overlay( ...
                intensity_t1, planimg.t1_image_orig, source_labels, parameters, ...
                {'y', backtransf_coordinates(2,2)}, ...
                backtransf_coordinates(1,:), backtransf_coordinates(2,:), backtransf_coordinates(3,:), ...
                'show_rectangles', 0, ...
                'grid_step', planimg.t1_header.PixelDimensions(1), ...
                'overlay_color_range', [0 max_val], ...
                'overlay_threshold_low', 0, ...
                'overlay_threshold_high', max_val, ...
                'rotation', 0, ...
                'do_transpose', do_t, 'flip_rows', flip_r, 'flip_cols', flip_c);
            draw_fov_outline(ax_ov, parameters, planimg, backtransf_coordinates(2,2), do_t, flip_r, flip_c);

            t_suffix = trans_suffix(ti, max_plots);
            output_plot_filename = fullfile(parameters.io.dir_img, ...
                sprintf('sub-%03d_%s_intensity_t1w%s%s.png', ...
                parameters.subject_id, parameters.simulation.medium, ...
                t_suffix, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png');
            close(h);
        catch ME_plot
            warning('plot_intensity_t1_overlay: plot_overlay failed for transducer %d (%s) — intensity overlay plot will be skipped.', ti, ME_plot.message);
        end
    end
end
end

function draw_fov_outline(ax, parameters, planimg, y_slice, do_t, flip_r, flip_c)
% Limited-FOV mode: outline the acoustic FOV box on the displayed coronal slice.
if ~isfield(parameters.grid, 'acoustic_fov_offset_vox') || ...
        isempty(parameters.grid.acoustic_fov_offset_vox)
    return
end
fov_mask = fov_box_slice(parameters.grid.acoustic_fov_offset_vox, ...
    parameters.grid.acoustic_fov_dims_vox, planimg.inv_transf, ...
    size(planimg.t1_image_orig), y_slice);
if do_t;   fov_mask = fov_mask';        end
if flip_r; fov_mask = flipud(fov_mask); end
if flip_c; fov_mask = fliplr(fov_mask); end
if ~any(fov_mask(:)); return; end
% draw on the top (intensity) layer so the outline stays visible;
% zero-pad so the outline also closes along the image border
[nr, nc] = size(fov_mask);
padded = zeros(nr + 2, nc + 2);
padded(2:end-1, 2:end-1) = fov_mask;
hold(ax, 'on');
contour(ax, 0:nc+1, 0:nr+1, padded, [0.5 0.5], ...
    'LineColor', [1 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1.2);
hold(ax, 'off');
end

function mask = fov_box_slice(fov_offset, fov_dims, inv_transf, t1_dims, y_slice)
% Coronal (dim-2) T1 slice of the acoustic FOV box. Each T1 voxel is mapped
% into full-grid voxels and tested against the box [offset, offset+dims-1].
if isstruct(inv_transf) && isfield(inv_transf, 'tdata')
    inv_transf = inv_transf.tdata.T;
end
[xi, zi] = ndgrid(1:t1_dims(1), 1:t1_dims(3));
sim_pts  = affine_apply_pts([xi(:), repmat(y_slice, numel(xi), 1), zi(:)], inv(inv_transf));
lo = fov_offset - 0.5;
hi = fov_offset + fov_dims - 0.5;
mask = reshape(all(sim_pts >= lo & sim_pts <= hi, 2), t1_dims([1 3]));
end
