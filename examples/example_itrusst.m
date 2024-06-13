cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/

% add paths
addpath('functions')
addpath('toolboxes/kwave') % set your kwave path here
addpath('toolboxes/kwaveArray') % set your kwave path here
addpath('toolboxes/Colormaps') % set your path to Colormaps files here
addpath('toolboxes/export_fig') % set your path to export_fig files here
addpath('toolboxes/yaml') % set your path to yaml files here

path_to_itrusst_results = '/home/visual/andche/STAFF_SCI/andche_sandbox/transcranial-ultrasound-benchmarks-master/KWAVE/PH1-BM1-SC1_KWAVE.mat';

parameters = load_parameters('itrusst_config.yaml');
parameters.data_path = '/project/3015999.02/andche_sandbox/TUS_sims/tusim/data/images_itrusst'; % path to save the data at

parameters.simulation_medium = 'water';
parameters.interactive = 1;
parameters.overwrite_files = 'ask';
subject_id = 999; % subject id here does not really matter as the calibrations are the same

% we start with 2d simulations, using the existing config that already has
% 2d grid and transducer position
parameters.results_filename_affix = '_2d';

output_pressure_file = single_subject_pipeline(subject_id, parameters);

% then the 3d simulations using 0.5 grid

parameters.transducer.pos_grid = [71 71 1];
parameters.default_grid_dims = [141 141 241];
parameters.grid_step_mm = 0.5;
parameters.n_sim_dims = 3;
parameters.results_filename_affix = '_3d';

output_pressure_file = single_subject_pipeline(subject_id, parameters);

% then the 3d simulations using 0.25 grid
parameters.transducer.pos_grid = [71 71 1]*2-1;
parameters.default_grid_dims = [141 141 241]*2-1;
parameters.grid_step_mm = 0.25;
parameters.n_sim_dims = 3;
parameters.results_filename_affix = '_3d_big';

output_pressure_file = single_subject_pipeline(subject_id, parameters);

% using kwaveArray and the 0.5 grid

parameters.transducer.pos_grid = [71 71 1];
parameters.default_grid_dims = [141 141 241];
parameters.grid_step_mm = 0.5;
parameters.n_sim_dims = 3;
parameters.use_kWaveArray = 1;
parameters.results_filename_affix = '_3d_kWaveArray';

output_pressure_file = single_subject_pipeline(subject_id, parameters);



% compute O'Neil solution and plot it along with comparisons

% define transducer parameters

velocity    = 0.04;   % [m/s] see p.3 in https://arxiv.org/pdf/2202.04552.pdf
frequency   = parameters.transducer.source_freq_hz;      % [Hz]
sound_speed = 1500;     % [m/s]
density     = 1000;     % [kg/m^3]

% define position vectors
axial_position   = (1:241)*0.5;       % [mm]
lateral_position = -35e-3:1e-4:35e-3;   % [m]

% evaluate pressure
[p_axial_oneil, p_lateral] = focusedBowlONeil(parameters.transducer.curv_radius_mm/1e3, parameters.transducer.Elements_OD_mm/1e3, velocity, ...
    frequency, sound_speed, density, (axial_position-0.5)*1e-3, lateral_position);


% plot
figure(2);
plot(axial_position, p_axial_oneil .* 1e-3);
hold on

affixes_list = ["2d","3d","3d_big","3d_kWaveArray","itrusst"];

for affix_n = 1:length(affixes_list)
    affix = affixes_list(affix_n);
    if affix~="itrusst"
        load(sprintf('%s/sim_outputs/sub-999_water_results_%s.mat', parameters.data_path , affix))
        p_max = gather(sensor_data.p_max_all/1000);
        if ndims(p_max) == 3
            p_max = squeeze(p_max(parameters.transducer.pos_grid(1),:,:));
            parameters.transducer.pos_grid = parameters.transducer.pos_grid(2:3);
        end
        pred_axial_pressure = squeeze(p_max(parameters.transducer.pos_grid(1), :));
        
        grid_step_mm = parameters.grid_step_mm;
        if grid_step_mm ~= 0.5
            pred_axial_pressure = pred_axial_pressure(1:(0.5/grid_step_mm):length(pred_axial_pressure));
        end
    else
        load(path_to_itrusst_results)
        p_max = p_amp'*1e-3;
        pred_axial_pressure = squeeze(p_max(71,:));
        grid_step_mm = 0.5;
    end
    figure(1)
    subplot(1,length(affixes_list), affix_n)
    
    imagesc((1:size(p_max,1))*grid_step_mm, (1:size(p_max,2))*grid_step_mm , p_max')
    title(affix, 'interpreter', 'none')  
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    axis image;
    cb = colorbar;
    title(cb, '[kPa]');
    
    figure(2)
    plot(axial_position, pred_axial_pressure,'--');
end

hold off 
colormap('lines')
xlabel('Axial Position [mm]');
ylabel('Pressure [kPa]');
legend(["O'Neil" affixes_list] , 'interpreter', 'none');

figure(2)
export_fig(sprintf('%s/axial_pressure.png',parameters.data_path))

figure(1)
export_fig(sprintf('%s/2d_pressure.png',parameters.data_path), '-native')

%% heating simulations
load(sprintf('%s/sim_outputs/sub-999_water_results_%s.mat', parameters.data_path , '3d_kWaveArray'))
