% DEMO_MATRIX_WATER  Free-water run of the tutorial matrix (Clover) transducer
%
% Simulates the matrix-array block of config/config_tutorial.yaml (two
% 52-element Fibonacci leaves, 6 mm bowl elements, ROC 100 mm, 250 kHz) in
% free water, for a set of requested foci: the natural focus and electronically
% steered targets. This is an empirical run of the setup, not a validation —
% there are no hydrophone measurements of this array.
%
% What this demo does:
%   1. Builds the element layout exactly as the pipeline does
%      (convert_to_element_pos + create_clover_array + create_matrix_karray)
%      and plots it.
%   2. Analytic preview (seconds, no k-Wave): sums the free-field contribution
%      of every element (point sources with the pipeline's steering phases) in
%      the x-z plane, and reports where the peak lands for each requested focus.
%   3. Optional (RUN_KWAVE = true): runs the PRESTUS water pipeline in 3D
%      (k-Wave, kWaveArray source) per focus and reports peak position,
%      focal error and ISPPA from the output tables. This is a 3D simulation
%      (~144 x 96 x 160 voxels, 104 elements); run it on a workstation/cluster.
%
% Outputs go to OUT_DIR. Compare the k-Wave peak positions with the analytic
% preview: the preview ignores element directivity and is only a sanity check
% of the geometry and phases.
%
% See also: CREATE_CLOVER_ARRAY, CREATE_MATRIX_KARRAY, CONVERT_TO_ELEMENT_POS,
%           doc/doc_transducer.md (matrix array)

%% ---- Settings ---------------------------------------------------------
RUN_KWAVE  = false;                 % true: run the 3D k-Wave water simulations
CODE_TYPE  = 'matlab_gpu';          % k-Wave backend for RUN_KWAVE ('matlab_cpu' | 'matlab_gpu' | 'cuda' ...)
OUT_DIR    = fullfile(tempdir, 'prestus_demo_matrix_water');
RES_MM     = 1.0;                   % grid resolution [mm] (6 points per wavelength at 250 kHz)
DIMS       = [144 96 160];          % grid [x y z] in voxels; leaves span about +-62 mm in x, foci up to 120 mm
TRANS_POS  = [72 48 5];             % parent bowl apex [voxels]; beam along +z
% Requested foci relative to TRANS_POS [mm]: natural focus (ROC 100 mm),
% axially and laterally steered targets
FOCI_MM    = [ 0 0 100;
               0 0  80;
               0 0 120;
              10 0 100;
               0 10 100];

%% ---- Paths and transducer ---------------------------------------------
prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'external')));
if ~isfolder(OUT_DIR); mkdir(OUT_DIR); end

tut = yaml.loadFile(fullfile(prestus_root, 'config', 'config_tutorial.yaml'), 'ConvertToArray', true);

parameters = load_parameters();
parameters.subject_id = 1;
parameters.path.sim  = OUT_DIR;
parameters.path.anat = OUT_DIR;
parameters.path.seg  = OUT_DIR;
parameters.simulation.medium      = 'water';
parameters.simulation.code_type   = CODE_TYPE;
parameters.simulation.interactive = 0;
parameters.simulation.debug       = 0;
parameters.io.overwrite_files     = 'always';
parameters.grid.resolution_mm     = RES_MM;
parameters.grid.default_dims      = DIMS;
parameters.grid.axisymmetric      = 0;
parameters.grid.use_kWaveArray    = 1;
parameters.modules.run_heating_sims     = 0;
parameters.modules.run_thermal_analysis = 0;
parameters.modules.run_posthoc_water_sims = 0;
parameters.modules.generate_report      = 0;
parameters.transducer.type    = 'matrix';
parameters.transducer.freq_hz = str2double(string(tut.transducer.freq_hz));
parameters.transducer.matrix  = tut.transducer.matrix;
parameters.transducer.trans_pos = TRANS_POS;
parameters.transducer.focus_pos = TRANS_POS + round(FOCI_MM(1,:) / RES_MM);
parameters = load_transducer_parameters(parameters);
tr = parameters.transducer;

%% ---- 1. Element layout (as in the pipeline) ----------------------------
kgrid = kWaveGrid(DIMS(1), RES_MM*1e-3, DIMS(2), RES_MM*1e-3, DIMS(3), RES_MM*1e-3);
vox2m = @(p) [kgrid.x_vec(p(1)); kgrid.y_vec(p(2)); kgrid.z_vec(p(3))];
tp_m  = vox2m(TRANS_POS);
[base_pos_m, tr] = convert_to_element_pos(parameters, tr, tp_m, tp_m);
[elem_pos_m, aim_m] = create_clover_array(parameters, tr.matrix, base_pos_m, tp_m, tp_m);
fprintf('Elements: %d per leaf x %d leaves = %d\n', size(base_pos_m, 2), ...
    tr.matrix.clover.n_leaves, size(elem_pos_m, 2));

h = figure('Color', 'w');
scatter3(elem_pos_m(1,:)*1e3, elem_pos_m(2,:)*1e3, elem_pos_m(3,:)*1e3, 20, 'filled'); hold on
pc = aim_m(:,1) * 1e3;
plot3(pc(1), pc(2), pc(3), 'r*', 'MarkerSize', 10);
axis equal; grid on; xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
title('Clover element layout (red: parent bowl centre)'); view(30, 20);
saveas(h, fullfile(OUT_DIR, 'matrix_layout.png')); close(h);

%% ---- 2. Analytic preview ------------------------------------------------
c0  = parameters.medium_properties.water.sound_speed;
k   = 2*pi*tr.freq_hz / c0;
xs  = kgrid.x_vec';  zs = kgrid.z_vec';
[X, Z] = ndgrid(xs, zs);
n_f = size(FOCI_MM, 1);
preview = table('Size', [n_f 6], 'VariableTypes', repmat({'double'}, 1, 6), ...
    'VariableNames', {'req_x_mm','req_z_mm','peak_x_mm','peak_z_mm','err_x_mm','err_z_mm'});   % mm from the bowl apex
h = figure('Color', 'w', 'Position', [100 100 300*n_f 320]);
for f = 1:n_f
    fp_vox = TRANS_POS + round(FOCI_MM(f,:) / RES_MM);
    [~, trf] = create_matrix_karray(kgrid, kWaveArray, parameters, tr, base_pos_m, TRANS_POS, fp_vox);
    phi = trf.matrix.elem_phase_rad(:);
    fp_m = vox2m(fp_vox);
    Y  = fp_m(2) * ones(size(X));
    p  = zeros(size(X));
    for e = 1:size(elem_pos_m, 2)
        r = sqrt((X - elem_pos_m(1,e)).^2 + (Y - elem_pos_m(2,e)).^2 + (Z - elem_pos_m(3,e)).^2);
        p = p + exp(1i * (phi(e) - k * r)) ./ max(r, 1e-3);
    end
    I = abs(p).^2;
    I(Z < tp_m(3) + 0.03) = 0;                  % ignore the element near field
    [~, idx] = max(I(:));
    [ix, iz] = ind2sub(size(I), idx);
    rq = [fp_m(1) - tp_m(1), fp_m(3) - tp_m(3)];
    pk = [xs(ix) - tp_m(1), zs(iz) - tp_m(3)];
    preview{f, :} = [rq, pk, pk - rq] * 1e3;
    subplot(1, n_f, f);
    imagesc(zs*1e3, xs*1e3, I / max(I(:))); axis image; hold on
    plot(fp_m(3)*1e3, fp_m(1)*1e3, 'w+', 'MarkerSize', 10, 'LineWidth', 1.5);
    title(sprintf('focus [%g %g %g] mm', FOCI_MM(f,:))); xlabel('z [mm]'); ylabel('x [mm]');
end
saveas(h, fullfile(OUT_DIR, 'matrix_preview_xz.png')); close(h);
disp('Analytic preview (point-source sum, x-z plane through each focus; mm from the bowl apex):');
disp(preview);
writetable(preview, fullfile(OUT_DIR, 'matrix_preview.csv'));

%% ---- 3. k-Wave free-water runs (optional) ----------------------------------
if ~RUN_KWAVE
    fprintf('RUN_KWAVE = false: skipping the 3D k-Wave runs. Outputs in %s\n', OUT_DIR);
    return
end
results = table();
for f = 1:n_f
    p_f = parameters;
    p_f.transducer.focus_pos = TRANS_POS + round(FOCI_MM(f,:) / RES_MM);
    p_f.io.output_affix = sprintf('_focus%02d', f);
    prestus_pipeline(p_f);
    csv = dir(fullfile(OUT_DIR, '**', sprintf('sub-%03d_water%s.csv', p_f.subject_id, p_f.io.output_affix)));
    if isempty(csv); warning('No output table for focus %d.', f); continue; end
    T = readtable(fullfile(csv(1).folder, csv(1).name));
    T.focus_idx = f;
    T.req_focus_mm = FOCI_MM(f,:);
    results = [results; T]; %#ok<AGROW>
end
writetable(results, fullfile(OUT_DIR, 'matrix_kwave_results.csv'));
fprintf('k-Wave results written to %s\n', fullfile(OUT_DIR, 'matrix_kwave_results.csv'));
