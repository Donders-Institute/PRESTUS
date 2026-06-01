function fig = neuronav_plot_placement(sub_id, session, active_sides, trans_pos, target_pos, target_map, parameters, output_dir, params_plot)
% NEURONAV_PLOT_PLACEMENT  Visualise transducer placements on a segmented head.
%
% Loads the SimNIBS segmentation for the subject and renders a two-panel
% 3-D figure (one panel per laterality: left / right) showing the transducer
% and target positions returned by NEURONAV_CONVERT_TRIGGER_TO_VOXELS.
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
%   parameters   - PRESTUS config struct.  Must contain parameters.seg_path.
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
    seg_file = fullfile(parameters.seg_path, sprintf('m2m_%s', sub_id), 'final_tissues.nii.gz');
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
    % Separate active_sides into left and right groups via target_map
    % ------------------------------------------------------------------
    sides_L = [];  trans_L = zeros(0,3);  target_L = zeros(0,3);  params_L = [];
    sides_R = [];  trans_R = zeros(0,3);  target_R = zeros(0,3);  params_R = [];

    for k = 1:numel(target_map)
        si = target_map(k).series_index;
        if ~ismember(si, active_sides)
            continue
        end
        ti  = target_map(k).transducer_id;
        name = target_map(k).target_name;

        % Build a parameters struct with only the relevant transducer at index 1.
        % Clamp ti to available transducers (single-transducer configs use index 1).
        p1 = params_plot;
        ti_clamped = min(ti, numel(params_plot.transducer));
        p1.transducer = coerce_transducer(params_plot.transducer(ti_clamped));

        if startsWith(name, 'L', 'IgnoreCase', true)
            trans_L  (end+1, :) = trans_pos(si, :);
            target_L (end+1, :) = target_pos(si, :);
            params_L = p1;
        elseif startsWith(name, 'R', 'IgnoreCase', true)
            trans_R  (end+1, :) = trans_pos(si, :);
            target_R (end+1, :) = target_pos(si, :);
            params_R = p1;
        end
    end

    % Fall back to first transducer if one side had no target_map entry
    if isempty(params_L)
        params_L = params_plot;
        params_L.transducer = coerce_transducer(params_plot.transducer(1));
    end
    if isempty(params_R)
        params_R = params_plot;
        params_R.transducer = coerce_transducer(params_plot.transducer(min(2, numel(params_plot.transducer))));
    end

    % ------------------------------------------------------------------
    % Plot
    % ------------------------------------------------------------------
    apply_deface = ~isfield(parameters, 'io') || ~isfield(parameters.io, 'deface_plots') || parameters.io.deface_plots;

    fig = figure('Position', [100 100 900 380], 'Visible', 'off');

    subplot(1, 2, 1)
    if ~isempty(trans_L)
        show_3d_head(seg_img, target_L, trans_L, params_L, pixel_size, coord_mesh, ...
            [0 0 0], [-90 10], 0, apply_deface)
        title(sprintf('%s — left', sub_id), 'Interpreter', 'none')
    else
        axis off
        title(sprintf('%s — left (no data)', sub_id), 'Interpreter', 'none')
    end

    subplot(1, 2, 2)
    if ~isempty(trans_R)
        show_3d_head(seg_img, target_R, trans_R, params_R, pixel_size, coord_mesh, ...
            [0 0 0], [90 10], 0, apply_deface)
        title(sprintf('%s — right', sub_id), 'Interpreter', 'none')
    else
        axis off
        title(sprintf('%s — right (no data)', sub_id), 'Interpreter', 'none')
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
