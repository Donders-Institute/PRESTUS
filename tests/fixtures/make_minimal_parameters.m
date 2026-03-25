function parameters = make_minimal_parameters(sim_path)
% MAKE_MINIMAL_PARAMETERS  Build a minimal valid parameters struct for water simulations.
%
%   parameters = make_minimal_parameters(sim_path)
%
%   Returns a parameters struct sufficient to run a complete water-medium
%   pipeline without any MRI data, SimNIBS, or GPU. Suitable for integration
%   tests and CI.
%
%   Input:
%     sim_path  - Directory for simulation outputs (created if missing).

    if nargin < 1 || isempty(sim_path)
        sim_path = fullfile(tempdir, 'prestus_test_outputs');
    end
    if ~exist(sim_path, 'dir'); mkdir(sim_path); end

    % --- Load defaults, then override ---
    parameters = load_parameters();

    % Identity
    parameters.subject_id = 1;
    parameters.platform   = 'matlab';

    % Paths
    parameters.path.sim  = sim_path;
    parameters.path.anat = sim_path;   % unused for water but must be set
    parameters.path.seg  = sim_path;

    % Simulation type
    parameters.simulation.medium      = 'water';
    parameters.simulation.code_type   = 'matlab_cpu';
    parameters.simulation.interactive = 0;
    parameters.simulation.precision   = 'single';

    % I/O
    parameters.io.overwrite_files   = 'always';
    parameters.io.overwrite_simnibs = 0;
    parameters.io.save_matrices     = 0;
    parameters.io.save_heatingvideo = 0;
    parameters.io.output_affix      = '_test';

    % Grid — small for speed
    parameters.grid.resolution_mm  = 1.0;
    parameters.grid.default_dims   = [72, 72, 128];
    parameters.grid.pml_size       = 10;
    parameters.grid.axisymmetric   = 0;
    parameters.grid.use_kWaveArray = 1;

    % Modules — water pipeline only, no thermal
    parameters.modules.run_grid_setup         = 1;
    parameters.modules.run_medium_setup       = 1;
    parameters.modules.run_source_setup       = 1;
    parameters.modules.run_acoustic_sims      = 1;
    parameters.modules.run_acoustic_analysis  = 1;
    parameters.modules.run_heating_sims       = 0;
    parameters.modules.run_thermal_analysis   = 0;
    parameters.modules.run_nifti_creation     = 1;
    parameters.modules.run_posthoc_water_sims = 0;
    parameters.modules.generate_report        = 0;

    % Transducer — single-element bowl (Ernie tutorial values, scaled down)
    parameters.transducer.n_elements         = 1;
    parameters.transducer.Elements_ID_mm     = 0;
    parameters.transducer.Elements_OD_mm     = 32.0;
    parameters.transducer.curv_radius_mm     = 63.2;
    parameters.transducer.dist_to_plane_mm   = 52.38;
    parameters.transducer.source_freq_hz     = 250e3;
    parameters.transducer.source_amp         = 91590;
    parameters.transducer.source_phase_deg   = 0;
    parameters.transducer.source_phase_rad   = 0;

    % Transducer & focus position (centre of grid, pointing along z)
    half = round(parameters.grid.default_dims / 2);
    parameters.transducer.trans_pos = [half(1), half(2), 5];
    parameters.transducer.focus_pos = [half(1), half(2), half(3)];

end
