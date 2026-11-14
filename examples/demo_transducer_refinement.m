% DEMO_TRANSDUCER_REFINEMENT  Optional refinement of the transducer library to your own hydrophone data
%
% PRESTUS ships a pre-calibrated library (config/calibration/). If you have
% measured axial profiles of YOUR transducer in a water tank, one call fits
% them and writes a user library entry; the shipped library is not modified.
%
% What this demo does:
%   1. Takes a real CSV if you set PROFILE_CSV below; otherwise it creates a
%      demonstration CSV in the standard format by simulating a transducer
%      whose hardware differs from the library (one weaker ring, a small
%      phase error) at three focal depths, with measurement noise.
%   2. Refines the library entry at the measured depths
%      (refine_transducer_library, per-depth mode) and shows the report.
%   3. Shows what PRESTUS then uses at a refined depth, between refined
%      depths and far from them (phase_model: refined / blend / static).
%   4. Prints the study-config lines that activate the refinement.
%
% Runtime: a few seconds per depth (analytical Rayleigh model, no k-Wave).
%
% See also: REFINE_TRANSDUCER_LIBRARY, FIT_TRANSDUCER_MODEL, CALIBRATION_GUI,
%           doc/doc_transducer_library.md (Level 2)

%% ---- Settings ---------------------------------------------------------
COMBO       = 'CTX_250_001_SC_105_010';   % transducer_driving-system key of the library entry
PROFILE_CSV = '';                         % '' = create a demonstration CSV; or your own file
DEPTHS_MM   = [25 35 45];                 % depths to refine (must be columns of the CSV)
OUT_LIBRARY = fullfile(tempdir, 'prestus_demo_library');   % your library folder

%% ---- Paths ------------------------------------------------------------
prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'external')));
equip = load_equipment_config();
lib0  = yaml.loadFile(library_entry_file(fullfile(prestus_root, 'config', 'calibration'), COMBO), 'ConvertToArray', true);
tran  = equip.trans.(char(lib0.meta.tran_serial));
ann   = tran.transducer.annular;

%% ---- 1. Profiles ------------------------------------------------------
if isempty(PROFILE_CSV)
    % "Your" transducer: the library model with ring 2 at 70 % and a 25° phase error on ring 3
    truth = library_static_model(lib0);
    n = ann.elem_n;
    g = ones(1, n); g(min(2, n)) = 0.7;
    ph = zeros(1, n); ph(min(3, n)) = 25;
    if isfield(truth, 'elem_amp_weight');           g  = g  .* truth.elem_amp_weight; end
    if isfield(truth, 'elem_phase_correction_deg'); ph = ph +  truth.elem_phase_correction_deg; end
    z_ep = (2:1:110)';
    H    = axial_transfer_matrix(tran, z_ep + ann.curv_radius_mm - ann.dist_geom_ep_mm, 1500, 994);
    k    = 2 * pi * tran.transducer.freq_hz / 1500;
    dist = steering_path_lengths(tran, DEPTHS_MM);
    rng(7);
    I = zeros(numel(z_ep), numel(DEPTHS_MM));
    for j = 1:numel(DEPTHS_MM)
        p = H * (0.05 * g(:) .* exp(1i * (k * dist(j, :)' + ph(:) * pi/180)));
        Ij = abs(p).^2;
        Ij = Ij / max(Ij(z_ep >= 10)) * 30;                     % matched focal ISPPA, 30 W/cm²
        I(:, j) = max(Ij .* (1 + 0.03 * randn(size(Ij))), 0);   % 3 % noise
    end
    PROFILE_CSV = fullfile(tempdir, sprintf('Axial_profiles_%s_demo.csv', COMBO));
    fid = fopen(PROFILE_CSV, 'w');
    fprintf(fid, ',Focal depths (mm)%s\n', repmat(',', 1, numel(DEPTHS_MM) - 1));
    fprintf(fid, 'Distance WRT Exit plane (mm)%s\n', sprintf(',%.1f', DEPTHS_MM));
    fprintf(fid, ['%.1f' repmat(',%.4f', 1, numel(DEPTHS_MM)) '\n'], [z_ep, I]');
    fclose(fid);
    fprintf('Demonstration profiles written to %s\n', PROFILE_CSV);
end

%% ---- 2. Refine ----------------------------------------------------------
[entry, fits, report] = refine_transducer_library(COMBO, PROFILE_CSV, OUT_LIBRARY, ...
    struct('Mode', 'per_depth', 'Depths', DEPTHS_MM, 'Plot', true));

%% ---- 3. What PRESTUS uses ------------------------------------------------
fprintf('\nPhases PRESTUS uses from %s:\n', entry);
for d = [DEPTHS_MM(2), mean(DEPTHS_MM(1:2)), DEPTHS_MM(end) + 3, DEPTHS_MM(end) + 20]
    out = load_transducer_from_library(COMBO, d, 30, equip, OUT_LIBRARY);
    fprintf('  %5.1f mm  →  %-8s (weight %.2f)\n', d, out.transducer.tx_calibration.phase_model, ...
        out.transducer.tx_calibration.refinement_weight);
end

%% ---- 4. Use it in a study ------------------------------------------------
fprintf(['\nAdd to your study config:\n' ...
         '  transducer:\n    serial: %s\n    combo:\n      ds_serial: %s\n' ...
         '    library_path: %s\n    focal_distance_ep: %g\n    target_isppa_wcm2: 30\n'], ...
    char(lib0.meta.tran_serial), char(lib0.meta.driving_system), OUT_LIBRARY, DEPTHS_MM(2));
fprintf('Report figure: %s\n', fullfile(OUT_LIBRARY, [COMBO '_refinement.png']));
