function plot_axial_profile(acoustic_Ipa, acoustic_pressure, medium_masks, ...
                             parameters, exitplane_pos, focus_pos, highlighted_pos)
% PLOT_AXIAL_PROFILE  Save axial pressure/intensity profile with tissue-layer shading
%
% Extracts the 1D profile along the beam axis (z-dimension) through the
% focus position, saves the numerical data as a CSV, and saves a two-panel
% figure (pressure + intensity) with coloured background bands indicating
% tissue layers.
%
% The axial-distance origin is the transducer EXIT PLANE (not the bowl apex),
% so distances match the target focal_distance_ep and the O'Neill/Rayleigh
% analytical convention.  The caller (acoustic_analysis) passes the exit-plane
% voxel position, which is the bowl apex shifted toward the focus by the bowl
% depth (sagitta = curv_radius_mm - dist_geom_ep_mm).
%
% Use as:
%   plot_axial_profile(acoustic_Ipa, acoustic_pressure, medium_masks, ...
%                      parameters, exitplane_pos, focus_pos, highlighted_pos)
%
% Input:
%   acoustic_Ipa      - [Nx×Ny×Nz] Isppa map [W/cm²]
%   acoustic_pressure - [Nx×Ny×Nz] peak pressure map [Pa]
%   medium_masks      - [Nx×Ny×Nz] tissue label map (from preproc_medium_mask)
%   parameters        - PRESTUS config struct
%   exitplane_pos     - [1×3] transducer exit-plane voxel position (axis origin)
%   focus_pos         - [1×3] intended focus voxel position
%   highlighted_pos   - [1×3] actual peak intensity voxel position

arguments
    acoustic_Ipa      (:,:,:) {mustBeNumeric}
    acoustic_pressure (:,:,:) {mustBeNumeric}
    medium_masks      (:,:,:) {mustBeNumeric}
    parameters        (1,1)   struct
    exitplane_pos     (1,3)   {mustBeNumeric}
    focus_pos         (1,3)   {mustBeNumeric}
    highlighted_pos   (1,3)   {mustBeNumeric}
end

res_mm = parameters.grid.resolution_mm;
Nz     = size(acoustic_Ipa, 3);

% Lateral position: use focus_pos (x,y), sweep all z
ix = max(1, min(round(focus_pos(1)), size(acoustic_Ipa,1)));
iy = max(1, min(round(focus_pos(2)), size(acoustic_Ipa,2)));

pressure_profile  = squeeze(acoustic_pressure(ix, iy, :));   % Pa
intensity_profile = squeeze(acoustic_Ipa(ix, iy, :));        % W/cm²
mask_profile      = squeeze(medium_masks(ix, iy, :));        % tissue label

% Distance from transducer exit plane [mm].
% exitplane_pos(3) is the exit-plane voxel (bowl apex shifted toward the focus
% by the bowl depth); positive z goes deeper into tissue.
z_vox  = (1:Nz)';
z_mm   = (z_vox - exitplane_pos(3)) * res_mm;   % mm from exit plane

% Focus and peak markers in mm
z_focus_mm       = (focus_pos(3)       - exitplane_pos(3)) * res_mm;
z_highlight_mm   = (highlighted_pos(3) - exitplane_pos(3)) * res_mm;

%% Save CSV
medium_labels = fieldnames(parameters.medium_properties);
tissue_names  = repmat({''}, Nz, 1);
for k = 1:Nz
    lbl = mask_profile(k);
    if lbl >= 1 && lbl <= numel(medium_labels)
        tissue_names{k} = medium_labels{lbl};
    else
        tissue_names{k} = 'water';
    end
end

T = table(z_mm, pressure_profile, intensity_profile, tissue_names, ...
    'VariableNames', {'depth_from_exit_plane_mm', 'pressure_Pa', ...
                      'intensity_Wcm2', 'tissue'});

csv_path = fullfile(parameters.io.dir_img, ...
    sprintf('sub-%03d_%s_axial_profile%s.csv', ...
    parameters.subject_id, parameters.simulation.medium, ...
    parameters.io.output_affix));
writetable(T, csv_path);

%% Build tissue shading patches
% Layer colours (matched to plot_overlay convention where possible)
layer_colors = struct( ...
    'water',           [0.85 0.93 1.00], ...
    'brain',           [0.95 0.80 0.80], ...
    'skin',            [0.98 0.88 0.70], ...
    'skull',           [0.80 0.80 0.80], ...
    'skull_cortical',  [0.72 0.72 0.72], ...
    'skull_trabecular',[0.90 0.90 0.90]);

default_color = [0.90 0.90 0.90];

% Find contiguous runs of the same tissue label
[run_labels, run_start, run_len] = tissue_runs(mask_profile);

%% Figure
h = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 800 520]);

ax1 = subplot(2,1,1);
draw_tissue_bands(ax1, run_labels, run_start, run_len, z_mm, Nz, ...
    medium_labels, layer_colors, default_color);
hold(ax1, 'on');
plot(ax1, z_mm, pressure_profile / 1e3, 'k-', 'LineWidth', 1.5);
xline(ax1, z_focus_mm,     '--b', 'LineWidth', 1.2, 'Label', 'target');
if abs(z_highlight_mm - z_focus_mm) > res_mm
    xline(ax1, z_highlight_mm, '--r', 'LineWidth', 1.2, 'Label', 'peak');
end
hold(ax1, 'off');
ylabel(ax1, 'Pressure (kPa)');
title(ax1, sprintf('Axial profile — sub-%03d  %s%s', ...
    parameters.subject_id, parameters.simulation.medium, ...
    strrep(parameters.io.output_affix, '_', ' ')));
set(ax1, 'XTickLabel', []);
xlim(ax1, [z_mm(1) z_mm(end)]);
grid(ax1, 'on');
box(ax1, 'on');

ax2 = subplot(2,1,2);
draw_tissue_bands(ax2, run_labels, run_start, run_len, z_mm, Nz, ...
    medium_labels, layer_colors, default_color);
hold(ax2, 'on');
plot(ax2, z_mm, intensity_profile, 'k-', 'LineWidth', 1.5);
xline(ax2, z_focus_mm,     '--b', 'LineWidth', 1.2);
if abs(z_highlight_mm - z_focus_mm) > res_mm
    xline(ax2, z_highlight_mm, '--r', 'LineWidth', 1.2);
end
hold(ax2, 'off');
xlabel(ax2, 'Depth from exit plane (mm)');
ylabel(ax2, 'I_{sppa} (W/cm²)');
xlim(ax2, [z_mm(1) z_mm(end)]);
grid(ax2, 'on');
box(ax2, 'on');

% Shared legend for tissue layers
legend_handles = build_layer_legend(run_labels, medium_labels, layer_colors, default_color);
if ~isempty(legend_handles)
    legend(ax2, legend_handles, 'Location', 'northeast', 'FontSize', 7);
end

%% Save
png_path = fullfile(parameters.io.dir_img, ...
    sprintf('sub-%03d_%s_axial_profile%s.png', ...
    parameters.subject_id, parameters.simulation.medium, ...
    parameters.io.output_affix));
set(h, 'InvertHardcopy', 'off');
saveas(h, png_path, 'png');
close(h);

fprintf('[plot_axial_profile] Saved: %s\n', png_path);
fprintf('[plot_axial_profile] Saved: %s\n', csv_path);
end

% -------------------------------------------------------------------------
function [labels, starts, lens] = tissue_runs(mask_profile)
% Find contiguous runs of the same value in a 1D label vector.
    d      = [true; diff(mask_profile) ~= 0; true];
    starts = find(d(1:end-1));
    ends   = find(d(2:end));
    lens   = ends - starts;
    labels = mask_profile(starts);
end

% -------------------------------------------------------------------------
function draw_tissue_bands(ax, run_labels, run_start, run_len, z_mm, Nz, ...
                            medium_labels, layer_colors, default_color)
% Draw background patch bands for each tissue run.
    ylims = [-Inf Inf];   % will be reset after data is plotted; use patch ydata trick
    for r = 1:numel(run_labels)
        lbl  = run_labels(r);
        iz1  = run_start(r);
        iz2  = run_start(r) + run_len(r) - 1;
        iz2  = min(iz2, Nz);
        z1   = z_mm(iz1);
        z2   = z_mm(iz2);
        if lbl >= 1 && lbl <= numel(medium_labels)
            tname = medium_labels{lbl};
        else
            tname = 'water';
        end
        if isfield(layer_colors, tname)
            fc = layer_colors.(tname);
        else
            fc = default_color;
        end
        patch(ax, [z1 z2 z2 z1], [0 0 1 1], fc, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.45, ...
            'Tag', tname, 'YLimInclude', 'off');
    end
    ylims = ylim(ax);   %#ok<NASGU>  y-extent updated by data plot after this call
end

% -------------------------------------------------------------------------
function handles = build_layer_legend(run_labels, medium_labels, layer_colors, default_color)
% Build one patch handle per unique tissue type present, for the legend.
    seen     = {};
    handles  = [];
    unique_labels = unique(run_labels);
    for u = 1:numel(unique_labels)
        lbl = unique_labels(u);
        if lbl >= 1 && lbl <= numel(medium_labels)
            tname = medium_labels{lbl};
        else
            tname = 'water';
        end
        if ismember(tname, seen); continue; end
        seen{end+1} = tname; %#ok<AGROW>
        if isfield(layer_colors, tname)
            fc = layer_colors.(tname);
        else
            fc = default_color;
        end
        h = patch(NaN, NaN, fc, 'EdgeColor', 'none', 'FaceAlpha', 0.45, ...
                  'DisplayName', strrep(tname, '_', ' '));
        handles = [handles, h]; %#ok<AGROW>
    end
end
