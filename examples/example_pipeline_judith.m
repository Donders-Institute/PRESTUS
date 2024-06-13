function example_pipeline_judith(subject_id)
    %% Initialization
    % start from tusim/code folder
    currDir = cd;
    cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/

    % add paths
    
    addpath('functions')
    addpath(genpath('toolboxes/'));
    
    parameters = load_parameters('judith_config.yaml');

    parameters.simulation_medium = 'water_and_skull';
    parameters.overwrite_files = 'never';
    parameters.overwrite_simnibs = 0;
    parameters.interactive = 0;
    parameters.use_kWaveArray = 1;
    parameters.results_filename_affix = '_kwave_array';
    if subject_id == 12
        parameters.csf_mask_expansion_factor = 15;
    end
    parameters.run_posthoc_water_sims = 1;

    single_subject_pipeline(subject_id, parameters)
    cd(currDir);
end