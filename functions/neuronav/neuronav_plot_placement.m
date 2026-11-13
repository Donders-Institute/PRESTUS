function fig = neuronav_plot_placement(sub_id, session, active_sides, trans_pos, target_pos, target_map, parameters, output_dir, params_plot)
% NEURONAV_PLOT_PLACEMENT  Visualise transducer placements on a segmented head.
%
% Loads the SimNIBS segmentation for the subject and renders a two-panel
% 3-D figure — both panels show every active transducer/target pair from
% target_map, viewed from two camera angles (left view, right view).
% Positions come from NEURONAV_CONVERT_TRIGGER_TO_VOXELS.
% Figures are saved as PNG to output_dir.
%
% Use as:
%   fig = neuronav_plot_placement(sub_id, session, active_sides, trans_pos, target_pos, ...
%             target_map, parameters, output_dir)
%
% Input:
%   sub_id       - subject identifier string (e.g. 'sub-001')
%   active_sides - integer array of series indices with valid data
%   trans_pos    - [Nx3] transducer voxel positions (indexed by series index)
%   target_pos   - [Nx3] target voxel positions (indexed by series index)
%   target_map   - struct array (one entry per coil).  Required fields:
%                    .series_index   integer — maps into trans_pos / target_pos rows
%                    .transducer_id  integer — which parameters.transducer() to use
%                    .target_name    char    — must start with 'L_' or 'R_'
%   parameters   - PRESTUS config struct.  Must contain parameters.path.seg.
%   output_dir   - folder where the PNG is saved
%   params_plot  - (optional) separate config struct whose transducer array
%                  is used for geometry/visualisation.  Pass when the main
%                  parameters use a legacy single-transducer format.
%                  Defaults to parameters if omitted.
%
% Output:
%   fig  - handle to the figure (invisible, already saved)
%
% See also: SHOW_3D_HEAD, NEURONAV_CONVERT_TRIGGER_TO_VOXELS

    if nargin < 9 || isempty(params_plot)
        params_plot = parameters;
    end

    % ------------------------------------------------------------------
    % Load segmentation
    % ------------------------------------------------------------------
    seg_file = fullfile(parameters.path.seg, sprintf('m2m_%s', sub_id), 'final_tissues.nii.gz');
    if ~exist(seg_file, 'file')
        warning('neuronav_plot_placement: segmentation not found for %s, skipping plot.', sub_id);
        fig = [];
        return
    end

    seg_img  = niftiread(seg_file);
    seg_info = niftiinfo(seg_file);
    pixel_size = seg_info.PixelDimensions(1);

    [gx, gy, gz] = ndgrid(1:size(seg_img,1), 1:size(seg_img,2), 1:size(seg_img,3));
    coord_mesh = [reshape(gx,[],1), reshape(gy,[],1), reshape(gz,[],1)];

    % ------------------------------------------------------------------
    % Gather every active target_map entry's position, per-pair transducer
    % geometry, and name, in target_map order. Both panels show every pair
    % (both lateralities); only the camera azimuth differs between panels.
    % ------------------------------------------------------------------
    all_trans  = zeros(0,3);
    all_target = zeros(0,3);
    all_td     = {};   % one transducer struct per pair, row-aligned with all_trans
    all_names  = {};

    for k = 1:numel(target_map)
        si = target_map(k).series_index;
        if ~ismember(si, active_sides)
            continue
        end
        ti   = target_map(k).transducer_id;
        name = target_map(k).target_name;

        % target_map(k).transducer_id is the hardware coil index, which is
        % 0-based (see NEURONAV_BUILD_TARGET_MAP), so it must be shifted by
        % one to become a valid 1-based index into params_plot.transducer.
        % Clamp to the available transducers (single-transducer configs only
        % have index 1).
        ti_clamped = max(1, min(ti + 1, numel(params_plot.transducer)));

        all_trans(end+1, :)  = trans_pos(si, :);
        all_target(end+1, :) = target_pos(si, :);
        all_td{end+1}        = coerce_transducer(params_plot.transducer(ti_clamped));
        all_names{end+1}     = display_target_name(name);
    end

    % ------------------------------------------------------------------
    % Plot — both pairs in each panel (viewed from two camera angles),
    % each pair rendered with its own transducer geometry.
    % ------------------------------------------------------------------
    apply_deface = ~isfield(parameters, 'io') || ~isfield(parameters.io, 'deface_plots') || parameters.io.deface_plots;

    n_pairs = size(all_trans, 1);
    plot_params = params_plot;
    if n_pairs > 0
        plot_params.transducer = [all_td{:}];
    end

    fig = figure('Position', [100 100 900 380], 'Visible', 'off');

    view_azimuths = [-90, 90];
    view_names    = {'left view', 'right view'};

    for i_view = 1:2
        subplot(1, 2, i_view);
        if n_pairs == 0
            axis off;
            title(sprintf('%s — no data', sub_id), 'Interpreter', 'none');
            continue
        end
        show_3d_head(seg_img, all_target, all_trans, plot_params, pixel_size, coord_mesh, ...
            [0 0 0], [view_azimuths(i_view) 10], 0, apply_deface);
        title(sprintf('%s — %s', sub_id, view_names{i_view}), 'Interpreter', 'none');

        % Legend proxies matching the color_list used inside show_3d_head,
        % one per actual target name (row-aligned with all_trans), so every
        % pair gets its own correctly colored legend entry regardless of how
        % many targets fall on each side.
        color_list = lines(max(7, n_pairs));
        hold on;
        leg_h = gobjects(1, n_pairs);
        for ip = 1:n_pairs
            leg_h(ip) = patch('XData', NaN, 'YData', NaN, ...
                'FaceColor', color_list(ip,:), 'EdgeColor', 'none', ...
                'DisplayName', all_names{ip});
        end
        legend(leg_h, all_names, 'Location', 'northeast', 'Interpreter', 'none');
    end

    saveas(fig, fullfile(output_dir, sprintf('%s_%s_placement', sub_id, session)), 'png');
    close(fig)
end

% -------------------------------------------------------------------------
function td = coerce_transducer(td)
% Ensure string fields that show_3d_head switches on are plain chars,
% not cell arrays (which can happen when YAML loads struct arrays or when
% configs lack explicit type fields).
    if isfield(td, 'type')
        if iscell(td.type),   td.type = char(td.type{1}); end
    else
        % Infer type from whichever geometry sub-struct is populated
        if isfield(td, 'annular') && ~isempty(td.annular)
            td.type = 'annular';
        elseif isfield(td, 'matrix') && ~isempty(td.matrix)
            td.type = 'matrix';
        else
            td.type = 'annular';  % safest default for this project
        end
    end
end

% -------------------------------------------------------------------------
function label = display_target_name(name)
% Convert an internal target_name (e.g. 'L_V1', 'R_M1_2') into a more
% readable legend label ('LEFT V1', 'RIGHT M1 (2)'). Cosmetic only -- this
% is never fed back into EntryTarget matching or any file-naming/lookup
% logic, which continue to use the original 'L_'/'R_'-prefixed target_name.
    name = char(name);
    suffix = '';
    tok = regexp(name, '^(.*)_(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        name = tok{1};
        suffix = sprintf(' (%s)', tok{2});
    end
    if startsWith(name, 'L_', 'IgnoreCase', true)
        label = ['LEFT ' upper(name(3:end))];
    elseif startsWith(name, 'R_', 'IgnoreCase', true)
        label = ['RIGHT ' upper(name(3:end))];
    else
        label = upper(name);
    end
    label = [label suffix];
end
