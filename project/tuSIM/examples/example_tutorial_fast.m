clear
cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/

% add paths
addpath('functions')
addpath(genpath('toolboxes')) % add toolboxes
addpath('/home/common/matlab/fieldtrip/qsub') % qsub on Donders HPC


subject_id = 1; % subject id - we use the brain of Ernie from SimNIBS example dataset, renamed to sub-001 in the data folder

%% Simulations in free water

parameters = load_parameters('tutorial_config_fast.yaml'); % load the configuration file

parameters.simulation_medium = 'water';

% start the simulations
% single_subject_pipeline(subject_id, parameters);

% load results
load(sprintf('%s/sim_outputs/sub-%03d_water_results%s.mat', parameters.data_path , subject_id, parameters.results_filename_affix));

% get maximum pressure
p_max = gather(sensor_data.p_max_all);
pred_axial_pressure = squeeze(p_max(parameters.transducer.pos_grid(1),parameters.transducer.pos_grid(2),:));

% compute O'Neil solution and plot it along with comparisons
% define transducer parameters

velocity = parameters.transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]

% define position vectors
axial_position   = (1:parameters.default_grid_dims(3))*0.5;       % [mm]

% evaluate pressure
[p_axial_oneil] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,4), ...
    parameters.transducer.source_phase_rad, parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);


% plot
figure;
plot(axial_position, p_axial_oneil .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');
hold on
plot(axial_position-(parameters.transducer.pos_grid(3)-1)*0.5, pred_axial_pressure.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4,'--');
hold off
legend('Analytic solution','Simulated results')

simulated_grid_adj_factor = max(pred_axial_pressure(:))/max(p_axial_oneil(:));


figure

imagesc((1:size(p_max,1))*parameters.grid_step_mm, (1:size(p_max,3))*parameters.grid_step_mm , squeeze(p_max(:,parameters.transducer.pos_grid(2),:))')
axis image;
colormap(getColorMap);
xlabel('Lateral Position [mm]');
ylabel('Axial Position [mm]');
axis image;
cb = colorbar;

fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil==max(p_axial_oneil)))

%% Optimize for a given distance and pressure

% reload the parameters
parameters = load_parameters('tutorial_config_fast.yaml'); 
parameters.focus_pos_t1_grid = [128, 142, 75]; % target we want to hit
parameters.transducer.pos_t1_grid = [128, 139, 21]; % approximate transducer location
parameters.expected_focal_distance_mm = norm(parameters.focus_pos_t1_grid- parameters.transducer.pos_t1_grid); % this T1 voxel is 1 mm^3, so the distance in grid units equals distance in mm
fprintf('Expected focal distance %.2f mm \n', parameters.expected_focal_distance_mm);

gs = GlobalSearch;
optimize_phases = @(phases) phase_optimization_annulus(phases, parameters, velocity, axial_position, parameters.expected_focal_distance_mm);

problem = createOptimProblem('fmincon','x0',[333 13 20]/180*pi,...
    'objective',optimize_phases,'lb',zeros(1,3),'ub',2*pi*ones(1,3));
[opt_phases, min_err] = run(gs,problem);

intensity_to_pressure_in_water_factor = 2*parameters.medium.water.sound_speed*parameters.medium.water.density/1e-4;

[p_axial_oneil_opt] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,4), ...
    [0 opt_phases], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);

max_pressure = max(p_axial_oneil_opt);
max_intensity = max_pressure.^2/intensity_to_pressure_in_water_factor;
desired_intensity = 20;
desired_pressure = sqrt(desired_intensity*intensity_to_pressure_in_water_factor);

opt_velocity = desired_pressure/max_pressure*velocity;

[p_axial_oneil_opt] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
    [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(opt_velocity,1,4), ...
    [0 opt_phases], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
    parameters.medium.water.density, (axial_position-0.5)*1e-3);

opt_sim_pressure = round(parameters.transducer.source_amp*desired_pressure/max_pressure/simulated_grid_adj_factor);

figure;
plot(axial_position, p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');
hold on
plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
hold off
xline(parameters.expected_focal_distance_mm, '--');
yline(desired_intensity, '--');
legend('Original', sprintf('Optimized for %2.f mm distance', parameters.expected_focal_distance_mm))

fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil_opt==max(p_axial_oneil_opt)))


%% Simulate again in water to check the optimization results
opt_parameters = parameters;
opt_parameters.transducer.source_amp = opt_sim_pressure;
opt_parameters.transducer.source_phase_rad = [0 opt_phases];
opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
opt_parameters.results_filename_affix = '_optimized_fast';

opt_parameters.simulation_medium = 'water';

opt_parameters.interactive = 1;
opt_parameters.overwrite_files = 'ask';

% addpath '/home/common/matlab/fieldtrip/qsub'
% qsubfeval(@single_subject_pipeline_wrapper, subject_id, opt_parameters, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
% 
single_subject_pipeline(subject_id, opt_parameters)

opt_res = load(sprintf('%s/sim_outputs/sub-%03d_water_results%s.mat', opt_parameters.data_path , subject_id, opt_parameters.results_filename_affix),'sensor_data','parameters');

% get maximum pressure
p_max = gather(opt_res.sensor_data.p_max_all);
pred_axial_pressure_opt = squeeze(p_max(opt_parameters.transducer.pos_grid(1), opt_parameters.transducer.pos_grid(2),:));

figure;
plot(axial_position, p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');
hold on
plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);

plot(axial_position-(opt_parameters.transducer.pos_grid(3)-1)*0.5, ...
    pred_axial_pressure_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);

hold off
xline(opt_parameters.expected_focal_distance_mm, '--')
legend('Original', sprintf('Optimized for %2.f mm distance, analytical', opt_parameters.expected_focal_distance_mm), sprintf('Optimized for %2.f mm distance, simulated', opt_parameters.expected_focal_distance_mm),'Max intensity')

fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(pred_axial_pressure==max(pred_axial_pressure)))


%% Simulations using skull and brain
opt_parameters = load_parameters('tutorial_config_fast.yaml'); 
opt_parameters.focus_pos_t1_grid = [128, 142, 75]; % target we want to hit
opt_parameters.transducer.pos_t1_grid = [128, 139, 21]; % approximate transducer location
opt_parameters.expected_focal_distance_mm = norm(opt_parameters.focus_pos_t1_grid - opt_parameters.transducer.pos_t1_grid); % this T1 voxel is 1 mm^3, so the distance in grid units equals distance in mm

opt_parameters.transducer.source_amp = opt_sim_pressure;
opt_parameters.transducer.source_phase_rad = [0 opt_phases];
opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
opt_parameters.results_filename_affix = '_optimized_fast';

opt_parameters.simulation_medium = 'layered'; % see default config for the list of mediums possible
opt_parameters.transducer = rmfield(opt_parameters.transducer, 'pos_grid');

opt_parameters.interactive = 0;
opt_parameters.overwrite_files = 'never';
opt_parameters.run_heating_sims = 1;

% single_subject_pipeline(subject_id, opt_parameters);

qsubfeval(@single_subject_pipeline_wrapper, subject_id, opt_parameters, 'timreq',  60*60*7,  'memreq',  80*(1024^3),  ...
     'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');

%% Postprocessing skull & brain data

sim_res = struct('type', {'water','water_and_skull','layered'},'isppa_max',[],'p_max',[],'parameters',[],'pred_axial_pressure',[],'pred_axial_intensity',[]);

for sim_i = 1:length(sim_res)
    cur_sim = sim_res(sim_i);
    sim_type = cur_sim.type;
    res = load(sprintf('%s/sim_outputs/sub-%03d_%s_results%s.mat', opt_parameters.data_path , subject_id, sim_type, '_optimized_fast'),'sensor_data','parameters','kwave_medium');
    cur_sim.p_max = gather(res.sensor_data.p_max_all);
    cur_sim.isppa_max = cur_sim.p_max.^2./(2*(res.kwave_medium.sound_speed.*res.kwave_medium.density)).*1e-4; 

    cur_sim.parameters = res.parameters;
    cur_sim.pred_axial_pressure = squeeze(cur_sim.p_max(res.parameters.transducer.pos_grid(1), res.parameters.transducer.pos_grid(2),:));
    cur_sim.pred_axial_intensity = squeeze(cur_sim.isppa_max(res.parameters.transducer.pos_grid(1), res.parameters.transducer.pos_grid(2),:));
    sim_res(sim_i) = cur_sim;

end


figure;
xlabel('Axial Position [mm]');
ylabel('Intensity [W/cm^2]');

hold on
for sim_i = 1:length(sim_res)
    cur_sim = sim_res(sim_i);
    axial_position = (1:cur_sim.parameters.grid_dims(3))*0.5;
    plot(axial_position-(cur_sim.parameters.transducer.pos_grid(3)-1)*0.5, cur_sim.pred_axial_intensity, 'DisplayName',cur_sim.parameters.simulation_medium);
end
hold off
legend show
xline(opt_parameters.expected_focal_distance_mm, '--')

%% Heating simulations

res = load(sprintf('%s/sim_outputs/sub-%03d_%s_results%s.mat', opt_parameters.data_path , subject_id, 'layered', '_optimized_fast'),'sensor_data','parameters','kwave_medium','kgrid','sensor','source');
load(sprintf('%s/sim_outputs/sub-%03d_%s_heating_res%s.mat', opt_parameters.data_path , subject_id, 'layered', '_optimized_fast'));
load(fullfile(opt_parameters.data_path, sprintf('sub-%03d_%s_after_cropping_and_smoothing%s.mat', subject_id, res.parameters.simulation_medium, res.parameters.results_filename_affix)), 'medium_masks')

heating_window_dims = ones(2,3);
for i = 1:2
    heating_window_dims(:,i) = [max(1, -opt_parameters.thermal.sensor_xy_halfsize + res.parameters.transducer.pos_grid(i)), min(res.parameters.grid_dims(i), opt_parameters.thermal.sensor_xy_halfsize + res.parameters.transducer.pos_grid(i))];
end

heating_window_dims(2,3) = res.parameters.grid_dims(3);

imshow(mat2gray(squeeze(medium_masks(:,res.parameters.transducer.pos_grid(2),:))))

plot_heating_sims(kwaveDiffusion, time_status_seq, res.parameters, heating_window_dims, medium_masks)


