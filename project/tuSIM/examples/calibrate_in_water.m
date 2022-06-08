function output_pressure_file = calibrate_in_water(extra_config, source_amp, source_phase_deg)
    
    % start from tusim/code folder
    cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/

    % add paths
    addpath('functions')
    addpath('toolboxes/kwave') % set your kwave path here
    addpath('toolboxes/Colormaps') % set your path to Colormaps files here
    addpath('toolboxes/export_fig') % set your path to export_fig files here
    addpath('toolboxes/yaml') % set your path to yaml files here

    parameters = load_parameters(extra_config);
    parameters.transducer.source_amp = repmat(source_amp, [1 parameters.transducer.n_elements]);
    parameters.simulation_medium = 'water';
    parameters.interactive = 0;
    parameters.overwrite_files = 'never';
    parameters.transducer.source_phase_deg = source_phase_deg;
	parameters.transducer.source_phase_rad = parameters.transducer.source_phase_deg/180*pi;

    parameters.results_filename_affix = sprintf('_amp_%s_phase_%s', string(source_amp), strjoin(string(source_phase_deg),'_'));

    subject_id = 999; % subject id here does not really matter as the calibrations are the same

    output_pressure_file = single_subject_pipeline(subject_id, parameters);
    close all
end