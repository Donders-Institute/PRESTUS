function plot_axial_profile(acoustic_Ipa, acoustic_pressure, medium_masks, ...
                             parameters, exitplane_pos, focus_pos, highlighted_pos, options)
% PLOT_AXIAL_PROFILE  Save axial pressure/intensity profile with tissue-layer shading
%
% Extracts the 1D profile along the beam axis (z-dimension) through the
% focus position, saves the numerical data as a CSV, and saves a two-panel
% figure (pressure + intensity) with coloured background bands indicating
% tissue layers.
%
% The axial-distance origin is the transducer EXIT PLANE (not the bowl apex),
% so distances match the target focal_distance_ep and the O'Neil/Rayleigh
% analytical convention.  The caller (acoustic_analysis) passes the exit-plane
% voxel position, which is the bowl apex shifted toward the focus by
% curv_radius_mm - dist_geom_ep_mm (ANNULAR_BOWL_DEPTH_MM).
%
% The intensity panel also shows the free-water on-axis intensity of the
% transducer settings used in the simulation, as a check of the realised
% focal intensity and depth against the calibrated transducer:
%   - O'Neil analytical model (single annular transducer);
%   - k-Wave free-water profile from WATER_BASELINE, when the baseline ran
%     (options.FreefieldAxial).
% Both are drawn at the drive level of the plotted field (options.IsppaScale)
% and placed on the plot's depth axis via their distance from the bowl apex,
% which assumes the beam runs along grid z (as the simulated profile does).
%
% Use as:
%   plot_axial_profile(acoustic_Ipa, acoustic_pressure, medium_masks, ...
%                      parameters, exitplane_pos, focus_pos, highlighted_pos)
%   plot_axial_profile(..., 'IsppaScale', s, 'FreefieldAxial', ff, 'FreeWater', tf)
%
% Input:
%   acoustic_Ipa      - [Nx×Ny×Nz] Isppa map [W/cm²]
%   acoustic_pressure - [Nx×Ny×Nz] peak pressure map [Pa]
%   medium_masks      - [Nx×Ny×Nz] tissue label map (from preproc_medium_mask)
%   parameters        - PRESTUS config struct
%   exitplane_pos     - [1×3] transducer exit-plane voxel position (axis origin)
%   focus_pos         - [1×3] intended focus voxel position
%   highlighted_pos   - [1×3] actual peak intensity voxel position
%
% Options (name-value):
%   IsppaScale     - intensity factor applied to the field after the
%                    simulation (APPLY_ISPPA_SCALING: scale_p^2); applied to
%                    the free-water references (default 1)
%   FreefieldAxial - struct with axial_distance_bowl [mm] and axial_intensity
%                    [W/cm²] from WATER_BASELINE (default: none)
%   FreeWater      - draw the free-water references at all (default true);
%                    false for fields whose drive level the references do not
%                    match (combined async/FOV fields, BabelBrain)
%
% See also: ACOUSTIC_ANALYSIS, WATER_BASELINE, FOCUSEDANNULUSONEIL

arguments
    acoustic_Ipa      (:,:,:) {mustBeNumeric}
    acoustic_pressure (:,:,:) {mustBeNumeric}
    medium_masks      (:,:,:) {mustBeNumeric}
    parameters        (1,1)   struct
    exitplane_pos     (1,3)   {mustBeNumeric}
    focus_pos         (1,3)   {mustBeNumeric}
    highlighted_pos   (1,3)   {mustBeNumeric}
    options.IsppaScale     (1,1) {mustBeNumeric, mustBePositive} = 1
    options.FreefieldAxial (1,1) struct = struct()
    options.FreeWater      (1,1) logical = true
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
% by the bowl depth). Grid z does not always increase toward the focus (the
% transducer can sit on either the low-z or high-z face), so the sign of the
% axis is set from the exit-plane -> focus direction rather than assumed.
prop_dir = sign(focus_pos(3) - exitplane_pos(3));
if prop_dir == 0
    prop_dir = 1;
end
z_vox  = (1:Nz)';
z_mm   = prop_dir * (z_vox - exitplane_pos(3)) * res_mm;   % mm from exit plane

% Focus and peak markers in mm
z_focus_mm       = prop_dir * (focus_pos(3)       - exitplane_pos(3)) * res_mm;
z_highlight_mm   = prop_dir * (highlighted_pos(3) - exitplane_pos(3)) * res_mm;

%% Free-water references (calibrated transducer settings)
[oneil_profile, kwave_profile] = free_water_references(parameters, z_mm, options);

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
if ~isempty(oneil_profile)
    T.freewater_oneil_Wcm2 = oneil_profile;
end
if ~isempty(kwave_profile)
    T.freewater_kwave_Wcm2 = kwave_profile;
end

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
stretch_tissue_bands(ax1);
ylabel(ax1, 'Pressure (kPa)');
title(ax1, sprintf('Axial profile — sub-%03d  %s%s', ...
    parameters.subject_id, parameters.simulation.medium, ...
    strrep(parameters.io.output_affix, '_', ' ')));
set(ax1, 'XTickLabel', []);
xlim(ax1, [min(z_mm) max(z_mm)]);
grid(ax1, 'on');
box(ax1, 'on');

ax2 = subplot(2,1,2);
draw_tissue_bands(ax2, run_labels, run_start, run_len, z_mm, Nz, ...
    medium_labels, layer_colors, default_color);
hold(ax2, 'on');
h_sim = plot(ax2, z_mm, intensity_profile, 'k-', 'LineWidth', 1.5, 'DisplayName', 'simulated');
ref_handles = gobjects(1, 0);
if ~isempty(oneil_profile)
    ref_handles(end+1) = plot(ax2, z_mm, oneil_profile, '--', 'LineWidth', 1.2, ...
        'Color', [0.2 0.2 0.2], 'DisplayName', 'free water (O''Neil)');
end
if ~isempty(kwave_profile)
    ref_handles(end+1) = plot(ax2, z_mm, kwave_profile, '-.', 'LineWidth', 1.2, ...
        'Color', [0.85 0.45 0], 'DisplayName', 'free water (k-Wave baseline)');
end
xline(ax2, z_focus_mm,     '--b', 'LineWidth', 1.2);
if abs(z_highlight_mm - z_focus_mm) > res_mm
    xline(ax2, z_highlight_mm, '--r', 'LineWidth', 1.2);
end
hold(ax2, 'off');
stretch_tissue_bands(ax2);
xlabel(ax2, 'Depth from exit plane (mm)');
ylabel(ax2, 'I_{sppa} (W/cm²)');
xlim(ax2, [min(z_mm) max(z_mm)]);
grid(ax2, 'on');
box(ax2, 'on');

% Shared legend for tissue layers
legend_handles = build_layer_legend(run_labels, medium_labels, layer_colors, default_color);
if ~isempty(ref_handles)
    legend_handles = [h_sim, ref_handles, legend_handles];
end
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
function [oneil_profile, kwave_profile] = free_water_references(parameters, z_mm, options)
% Free-water on-axis intensity [W/cm²] at the plot's depths (mm from the exit
% plane); [] when unavailable. Both references are functions of the distance
% from the bowl apex, which lies a bowl depth behind the exit plane.
    oneil_profile = [];
    kwave_profile = [];
    if ~options.FreeWater
        return;
    end
    % The references are optional overlays: without a (typed) transducer the
    % O'Neil curve is skipped and depths are taken as distances from the bowl
    n_tr = 0; tr = struct(); tr_type = '';
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        n_tr = numel(parameters.transducer);
        tr   = parameters.transducer(1);
        if isfield(tr, 'type') && ~isempty(tr.type); tr_type = char(tr.type); end
    end
    bowl_depth_mm = 0;
    if ~isempty(tr_type) && isfield(tr, tr_type)
        bowl_depth_mm = annular_bowl_depth_mm(tr.(tr_type));
    end
    d_bowl_mm = z_mm + bowl_depth_mm;

    % O'Neil model of the simulated annular transducer (elem_p_amp, phases).
    % A coherent multi-transducer field is not one transducer's free-water field.
    if n_tr == 1 && strcmp(tr_type, 'annular') && ...
            isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'water')
        try
            an  = tr.annular;
            rho = parameters.medium_properties.water.density;
            c   = parameters.medium_properties.water.sound_speed;
            if isfield(an, 'elem_phase_rad') && ~isempty(an.elem_phase_rad)
                phase_rad = an.elem_phase_rad(:)';
            else
                phase_rad = deg2rad(an.elem_phase_deg(:)');
            end
            vel = an.elem_p_amp(:)' / (rho * c);
            if isscalar(vel); vel = repmat(vel, 1, an.elem_n); end
            in_front = d_bowl_mm > 0;
            if isinf(an.curv_radius_mm)
                % Flat element(s): on-axis plane-piston solution per ring,
                % p = rho c v (exp(-ik z) - exp(-ik sqrt(z^2 + a^2))), summed coherently
                k  = 2 * pi * tr.freq_hz / c;
                zz = d_bowl_mm(in_front)' * 1e-3;
                pc = zeros(size(zz));
                for el = 1:an.elem_n
                    ri = an.elem_id_mm(el) / 2e3; ro = an.elem_od_mm(el) / 2e3;
                    pc = pc + rho * c * vel(el) * exp(1i * phase_rad(el)) .* ...
                        (exp(-1i * k * sqrt(zz.^2 + ri^2)) - exp(-1i * k * sqrt(zz.^2 + ro^2)));
                end
                p_ref = abs(pc);
            else
                p_ref = focusedAnnulusONeil(an.curv_radius_mm / 1e3, ...
                    [an.elem_id_mm(:)'; an.elem_od_mm(:)'] / 1e3, vel, phase_rad, ...
                    tr.freq_hz, c, rho, d_bowl_mm(in_front)' * 1e-3);
            end
            oneil_profile = NaN(size(z_mm));
            oneil_profile(in_front) = options.IsppaScale * p_ref(:).^2 / (2 * c * rho) * 1e-4;
        catch ME
            warning('plot_axial_profile:oneil', ...
                'plot_axial_profile: O''Neil free-water reference not computed (%s).', ME.message);
            oneil_profile = [];
        end
    end

    % k-Wave free-water profile from the baseline run
    ff = options.FreefieldAxial;
    if isfield(ff, 'axial_intensity') && numel(ff.axial_intensity) > 1
        kwave_profile = options.IsppaScale * interp1(ff.axial_distance_bowl(:), ...
            ff.axial_intensity(:), d_bowl_mm, 'linear', NaN);
    end
end

% -------------------------------------------------------------------------
function [labels, starts, lens] = tissue_runs(mask_profile)
% Find contiguous runs of the same value in a 1D label vector.
    d      = [true; diff(mask_profile) ~= 0; true];
    starts = find(d(1:end-1));
    ends   = find(d(2:end));
    lens   = ends - starts + 1;
    labels = mask_profile(starts);
end

% -------------------------------------------------------------------------
function draw_tissue_bands(ax, run_labels, run_start, run_len, z_mm, Nz, ...
                            medium_labels, layer_colors, default_color)
% Draw background patch bands for each tissue run (y extent set by
% stretch_tissue_bands once the data are plotted).
    for r = 1:numel(run_labels)
        lbl  = run_labels(r);
        iz1  = run_start(r);
        iz2  = run_start(r) + run_len(r) - 1;
        iz2  = min(iz2, Nz);
        % band edges half a voxel beyond the run's first/last sample (z_mm
        % may run either way), so neighbouring bands meet without a gap
        dz   = (z_mm(end) - z_mm(1)) / max(Nz - 1, 1);
        z1   = z_mm(iz1) - dz / 2;
        z2   = z_mm(iz2) + dz / 2;
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
end

% -------------------------------------------------------------------------
function stretch_tissue_bands(ax)
% Extend the tissue bands over the full y-range of the plotted data. The bands
% are drawn first with a placeholder extent and excluded from autoscaling.
    yl = ylim(ax);
    ylim(ax, yl);
    bands = findobj(ax, 'Type', 'patch', '-not', 'Tag', '');
    set(bands, 'YData', [yl(1); yl(1); yl(2); yl(2)]);
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
