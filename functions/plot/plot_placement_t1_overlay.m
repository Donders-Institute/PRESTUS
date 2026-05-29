function plot_placement_t1_overlay(parameters, trans_pos, focus_pos, label)
% PLOT_PLACEMENT_T1_OVERLAY  Quick placement QC figure on the native T1
%
% Loads the subject T1 and produces a 2×3 grid of orthogonal slices:
%   Row 1 — slices through the transducer position (scalp view)
%   Row 2 — slices through the focus position (brain view)
% Each panel overlays the transducer bowl mask (red) and focus marker (green)
% via plot_t1_with_transducer, with a plain-marker fallback if the transducer
% geometry is not yet initialised. Images are oriented superior-at-top.
% Saves to parameters.io.dir_img as sub-NNN_placement_<label>.png.
%
% Use as:
%   plot_placement_t1_overlay(parameters, trans_pos, focus_pos, label)
%
% Input:
%   parameters - (1,1) simulation parameters struct with fields:
%                  path.anat, path.t1_pattern, subject_id,
%                  io.dir_img, io.output_affix
%   trans_pos  - [1×3] transducer position in T1 voxel space
%   focus_pos  - [1×3] focus position in T1 voxel space
%   label      - short string appended to the output filename (e.g. 'plantus')

    filename_t1 = fullfile(parameters.path.anat, ...
        sprintf(parameters.path.t1_pattern, parameters.subject_id));
    if ~isfile(filename_t1)
        warning('plot_placement_t1_overlay: T1 not found, skipping QC figure:\n  %s', filename_t1);
        return;
    end

    t1  = double(niftiread(filename_t1));
    hdr = niftiinfo(filename_t1);
    vox = hdr.PixelDimensions(1);  % assume isotropic; used for axis labels only

    tp = round(trans_pos(:)');
    fp = round(focus_pos(:)');

    % Clamp slice indices to volume bounds
    sz = size(t1);
    tp_clamped = max(1, min(tp, sz));
    fp_clamped = max(1, min(fp, sz));

    dims       = {'X (sag)','Y (cor)','Z (ax)'};
    slice_axes = {'x','y','z'};
    inplane    = {[2 3],[1 3],[1 2]};

    % Row 1: slices through transducer (scalp view); Row 2: slices through focus (brain view)
    row_defs = { ...
        tp_clamped, 'transducer'; ...
        fp_clamped, 'focus'       ...
    };

    h = figure('Visible','off','Position',[0 0 900 600]);
    for r = 1:2
        ref    = row_defs{r, 1};
        rlabel = row_defs{r, 2};
        for d = 1:3
            ax = subplot(2, 3, (r-1)*3 + d);
            ia = inplane{d};

            % --- Get image data ---
            % Render transducer bowl via plot_t1_with_transducer; fall back to
            % plain grayscale if the transducer geometry is not yet initialised.
            use_rgb = false;
            try
                rgb = plot_t1_with_transducer(t1, vox, tp, fp, parameters, ...
                    'slice_dim', d, 'slice_ind', ref(d));
                use_rgb = true;
            catch
                idx = repmat({':'}, 1, 3);
                idx{d} = ref(d);
                sl = squeeze(t1(idx{:}));
            end

            % --- Compute orientation transforms ---
            % nifti_slice_orientation returns (do_transpose, flip_rows, flip_cols)
            % so that the resulting image has superior at top and patient-left on
            % viewer's left.  Markers are adjusted to match.
            [do_tp, fl_r, fl_c] = nifti_slice_orientation(hdr, slice_axes{d});

            % Initial marker coords [row, col] before any transform
            tp_rc = [tp_clamped(ia(1)), tp_clamped(ia(2))];
            fp_rc = [fp_clamped(ia(1)), fp_clamped(ia(2))];

            if do_tp
                if use_rgb
                    rgb = permute(rgb, [2 1 3]);
                else
                    sl = sl';
                end
                tp_rc = fliplr(tp_rc);
                fp_rc = fliplr(fp_rc);
            end

            % Size after potential transpose
            if use_rgb
                nr = size(rgb, 1);  nc = size(rgb, 2);
            else
                nr = size(sl, 1);   nc = size(sl, 2);
            end

            if fl_r
                if use_rgb; rgb = rgb(end:-1:1,:,:); else; sl = sl(end:-1:1,:); end
                tp_rc(1) = nr - tp_rc(1) + 1;
                fp_rc(1) = nr - fp_rc(1) + 1;
            end
            if fl_c
                if use_rgb; rgb = rgb(:,end:-1:1,:); else; sl = sl(:,end:-1:1); end
                tp_rc(2) = nc - tp_rc(2) + 1;
                fp_rc(2) = nc - fp_rc(2) + 1;
            end

            % --- Display ---
            if use_rgb
                image(ax, rgb);
            else
                imagesc(ax, sl); colormap(ax, gray);
                hold(ax,'on');
                plot(ax, tp_rc(2), tp_rc(1), 's', 'MarkerSize', 10, ...
                     'Color', [0.2 0.4 1], 'LineWidth', 2);
                plot(ax, fp_rc(2), fp_rc(1), '+r', 'MarkerSize', 12, 'LineWidth', 2);
                hold(ax,'off');
            end
            axis(ax,'image','off');
            title(ax, sprintf('%s  %s=%d', rlabel, dims{d}, ref(d)), 'FontSize', 8);
        end
    end

    sgtitle(sprintf('sub-%03d  placement: %s  |  trans=[%d %d %d]  focus=[%d %d %d]', ...
        parameters.subject_id, label, tp(1),tp(2),tp(3), fp(1),fp(2),fp(3)), ...
        'FontSize', 9, 'Interpreter', 'none');

    % dir_img is only populated when the full pipeline's path_log_setup has run.
    % Fall back to dir_output (always present) when running placement standalone.
    if isfield(parameters.io, 'dir_img') && ~isempty(parameters.io.dir_img)
        out_dir = parameters.io.dir_img;
    elseif isfield(parameters.io, 'dir_output') && ~isempty(parameters.io.dir_output)
        out_dir = parameters.io.dir_output;
    else
        warning('plot_placement_t1_overlay: no output directory found in parameters.io, skipping QC figure.');
        return;
    end
    if ~isfolder(out_dir); mkdir(out_dir); end

    affix = '';
    if isfield(parameters.io, 'output_affix'); affix = parameters.io.output_affix; end
    out_file = fullfile(out_dir, ...
        sprintf('sub-%03d_placement_%s%s.png', parameters.subject_id, label, affix));
    saveas(h, out_file, 'png');
    close(h);
    fprintf('Placement QC figure saved:\n  %s\n', out_file);
end
