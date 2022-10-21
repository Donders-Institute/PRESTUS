clear

% Delete if you have rights to add paths to Matlab
cd /home/mrphys/kenvdzee/Documents/MATLAB/
addpath(genpath('SimNIBS-3.2'))
addpath(genpath('k-wave'))

% Change path to tuSIM folder
cd /home/mrphys/kenvdzee/orca-lab/project/tuSIM

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% Read table with old and new names
%readtable(parameters.data_path + 'localite_new_filenames.csv', 'Delimiter',',');

% Set config files
config_left_transducer = 'sjoerd_config_opt_CTX250-011_64.5mm.yaml';
config_right_transducer = 'sjoerd_config_opt_CTX250-001_64.5mm.yaml';
overwrite_option = 'never';

parameters = load_parameters(config_left_transducer);
out_folder_gen = parameters.data_path+'sim_outputs/';
reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

% Create list of files in datafolder (location can be changed in config file)
files = dir(parameters.data_path);
notsubjects = (strlength(extractfield(files,'name')) == 7);
files(~notsubjects) = [];
subject_list = [];
for i = 1:length(files)
    fname = files(i).name;
    if regexp(fname, 'sub-\d+')
        subject_list = [subject_list str2num(fname(5:7))];
    end
end
% [1, 3, 4, 5, 8, 9, 10, 11, 14, 17, 18, 19];
subject_list = [5, 8, 9, 10, 11, 14, 17, 18, 19];

for subject_id = subject_list

% Setting folder locations
subj_folder = fullfile(parameters.data_path,sprintf('sub-%1$03d/', subject_id));
filename_t1 = dir(sprintf(fullfile(parameters.data_path,parameters.t1_path_template), subject_id));
t1_header = niftiinfo(fullfile(filename_t1.folder,filename_t1.name));
t1_image = niftiread(fullfile(filename_t1.folder,filename_t1.name));

if parameters.output_folder == 1
    out_folder = fullfile(out_folder_gen, sprintf('sub-%1$03d/', subject_id));
else
    out_folder = out_folder_gen;
end

% Left transducer localite files
%trig_mark_files = dir(sprintf('%ssub-%03d/TriggerMarkers_Coil0_*.xml',parameters.data_path,subject_id));
trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_left*.xml',parameters.data_path, subject_id, subject_id));

% sort by datetime
extract_dt = @(x) datetime(x.name(28:end-4),'InputFormat','yyyyMMddHHmmssSSS'); %22-28
[~,idx] = sort([arrayfun(extract_dt,trig_mark_files)],'descend');
trig_mark_files = trig_mark_files(idx);

% Translate transducer trigger markers to raster positions
[left_trans_ras_pos, left_amygdala_ras_pos] = get_trans_pos_from_trigger_markers(fullfile(trig_mark_files(1).folder, trig_mark_files(1).name), 5, ...
    reference_to_transducer_distance, parameters.expected_focal_distance_mm);
left_trans_pos = ras_to_grid(left_trans_ras_pos, t1_header);
left_amygdala_pos = ras_to_grid(left_amygdala_ras_pos, t1_header);

% Right transducer localite file
%trig_mark_files = dir(sprintf('%ssub-%03d/TriggerMarkers_Coil1_*.xml',parameters.data_path,subject_id));
trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_right*.xml',parameters.data_path, subject_id, subject_id));

% sort by datetime
extract_dt = @(x) datetime(x.name(29:end-4),'InputFormat','yyyyMMddHHmmssSSS'); %22-29
[~,idx] = sort([arrayfun(extract_dt,trig_mark_files)],'descend');
trig_mark_files = trig_mark_files(idx);

% Translate transducer trigger markers to raster positions
[right_trans_ras_pos, right_amygdala_ras_pos] = get_trans_pos_from_trigger_markers(fullfile(trig_mark_files(1).folder, trig_mark_files(1).name), 5, ...
    reference_to_transducer_distance, parameters.expected_focal_distance_mm);
right_trans_pos = ras_to_grid(right_trans_ras_pos, t1_header);
right_amygdala_pos = ras_to_grid(right_amygdala_ras_pos, t1_header);

% Generate plots with both transducers and targets separately
imshowpair(plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), left_trans_pos, left_amygdala_pos, parameters), plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), right_trans_pos, right_amygdala_pos, parameters),'montage');

% Index transducer locations for simulation selection
transducers = [left_trans_pos right_trans_pos];
targets = [left_amygdala_pos right_amygdala_pos];
target_names = {'left_amygdala', 'right_amygdala'};

% Simulations for left amygdala
% Loading parameters
parameters = load_parameters(config_left_transducer);
out_folder_gen = parameters.data_path+'sim_outputs/';
parameters.overwrite_files = overwrite_option;

% Select 'layered' when simulating the transmission in a skull
parameters.simulation_medium = 'layered';
reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

target_id = 1;
parameters.transducer.pos_t1_grid = transducers(:,target_id)';
parameters.focus_pos_t1_grid = targets(:,target_id)';
parameters.results_filename_affix = sprintf('_target_%s', target_names{target_id});
parameters.interactive = 0;
single_subject_pipeline(subject_id, parameters)

% Simulations for right amygdala
parameters = load_parameters(config_right_transducer);
out_folder_gen = parameters.data_path+'sim_outputs/';
parameters.overwrite_files = overwrite_option;

% Select 'layered' when simulating the transmission in a skull
parameters.simulation_medium = 'layered';
reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

target_id = 2;
parameters.transducer.pos_t1_grid = transducers(:,target_id)';
parameters.focus_pos_t1_grid = targets(:,target_id)';
parameters.results_filename_affix = sprintf('_target_%s', target_names{target_id});
parameters.interactive = 0;
single_subject_pipeline(subject_id, parameters)

end