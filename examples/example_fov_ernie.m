%% EXAMPLE_FOV_ERNIE  Partial acoustic field of view against the full grid
% Runs the tutorial's acoustic simulation twice: on the full grid (the reference, affix
% 'fullgrid') and on a partial field of view (example_fov_ernie.yaml, affix 'fov'). Then it
% prints the metrics both runs write, so you can judge whether the box is large enough for
% your head and target before using it for many runs. See doc/doc_partial-fov.md.
%
% Acoustic only (as in tutorial_ernie.yaml): the thermal run uses the full grid by default, so
% the partial field of view saves time in the acoustic stage.
%
% Before running: set the paths in tutorial_ernie.yaml and run its segmentation step
% (section 1 of tutorial_ernie.m). For your own image, change path.* there; nothing here.
% Not yet run end to end.

prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'config')));
addpath(genpath(fullfile(prestus_root, 'external')));

cfg_dir = fullfile(prestus_root, 'examples');

%% 1. Reference: full grid
p_full = load_parameters('tutorial_ernie.yaml', cfg_dir, ...
                         {'subject_id', 1; 'io.output_affix', 'fullgrid'});
t = tic; prestus_pipeline_start(p_full); t_full = toc(t);

%% 2. Partial field of view
p_fov = load_parameters('tutorial_ernie.yaml', cfg_dir, 'example_fov_ernie.yaml', cfg_dir, ...
                        {'subject_id', 1});
t = tic; prestus_pipeline_start(p_fov); t_fov = toc(t);

%% 3. Compare
% Wall times are meaningful with platform: matlab only (a cluster run returns at submission).
csv_of = @(p) fullfile(get_output_dir(p), sprintf('sub-%03d_%s%s.csv', ...
    p.subject_id, p.simulation.medium, p.io.output_affix));
full = readtable(csv_of(p_full));
fov  = readtable(csv_of(p_fov));

cols = {'Isppa_brain', 'Isppa_skull', 'Isppa_skin', 'Psptp_brain', 'MI_tc', ...
        'real_focal_distance_ep_mm', 'focal_fwhm_axial_mm', 'focal_fwhm_lateral_max_mm', ...
        'focal_vol_6dB_mm3'};
fprintf('\n%-28s %12s %12s %9s\n', 'Metric', 'Full grid', 'Partial FOV', 'Diff (%)');
for i = 1:numel(cols)
    c = cols{i};
    if ~ismember(c, full.Properties.VariableNames) || ~ismember(c, fov.Properties.VariableNames)
        continue
    end
    a = full.(c)(1); b = fov.(c)(1);
    fprintf('%-28s %12.4g %12.4g %9.1f\n', c, a, b, 100 * (b - a) / a);
end
if strcmp(p_fov.platform, 'matlab')
    fprintf('\nWall time: full grid %.0f s, partial FOV %.0f s\n', t_full, t_fov);
end
% A drop in Isppa_skull or Isppa_skin, or a shifted focus, suggests the box cuts the beam or the
% near field: set grid.acoustic_fov_diameter_mm / acoustic_fov_length_mm larger than 'auto' gives
% (doc/doc_partial-fov.md, "Choosing FOV Dimensions").
