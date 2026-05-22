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

max_plots = min(2, numel(parameters.transducer));
if numel(parameters.transducer) > max_plots
    warn('More than two transducers: intensity-over-T1 plots will be created only for the first 2 transducers');
end

if isfield(results_acoustic, 'Isppa_brain') && ~isempty(results_acoustic.Isppa_brain) && ~isnan(results_acoustic.Isppa_brain)
    max_val = results_acoustic.Isppa_brain;
else
    max_val = results_acoustic.Isppa;
end

for ti = 1:max_plots
    tpos_sim = parameters.transducer(ti).trans_pos;
    fpos_sim = parameters.transducer(ti).focus_pos;

    backtransf_coordinates = round(affine_apply_pts( ...
        [tpos_sim; fpos_sim; highlighted_pos], planimg.inv_transf));

    [~, source_labels] = transducer_setup( ...
        parameters.transducer(1), ...
        backtransf_coordinates(1,:), ...
        backtransf_coordinates(2,:), ...
        size(planimg.t1_image_orig), ...
        planimg.t1_header.PixelDimensions(1), parameters);

    if ~isempty(planimg.t1_header) && isfield(planimg.t1_header, 'Transform')
        [do_t, flip_r, flip_c] = nifti_slice_orientation(planimg.t1_header, 'y');
    else
        [do_t, flip_r, flip_c] = deal(false, false, false);
    end

    [~,~,~,~,~,~,~,h] = plot_overlay( ...
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

    trans_suffix = '';
    if max_plots > 1; trans_suffix = sprintf('_T%02d', ti); end
    output_plot_filename = fullfile(parameters.io.dir_img, ...
        sprintf('sub-%03d_%s_intensity_t1%s%s.png', ...
        parameters.subject_id, parameters.simulation.medium, ...
        trans_suffix, parameters.io.output_affix));
    saveas(h, output_plot_filename, 'png');
    close(h);
end
end
