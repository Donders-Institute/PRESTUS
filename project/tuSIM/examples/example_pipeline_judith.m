function example_pipeline_judith(subject_id)
    %% Initialization
    % start from tusim/code folder
    cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/

    % add paths
    addpath('.');
    addpath('functions')
    addpath('toolboxes/kwave') % set your kwave path here
    addpath('toolboxes/Colormaps') % set your path to Colormaps files here
    addpath('toolboxes/export_fig') % set your path to export_fig files here
    addpath('toolboxes/yaml') % set your path to yaml files here

    parameters = load_parameters('judith_config.yaml');

    parameters.simulation_medium = 'water_and_skull';
    parameters.overwrite_files = 'always';
    parameters.overwrite_simnibs = 0;

    parameters.run_posthoc_water_sims = 1;

    single_subject_pipeline(subject_id, parameters)
end