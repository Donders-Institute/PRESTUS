% Delete if you have rights to add paths to Matlab
%cd /home/mrphys/kenvdzee/Documents/
cd /home/action/elecar/
addpath(genpath('SimNIBS-4.0'))
%cd /home/mrphys/kenvdzee/Documents/MATLAB/
cd /home/action/elecar/Documents/
addpath(genpath('k-Wave'))

% Change path to tuSIM folder
%cd /home/mrphys/kenvdzee/orca-lab/project/tuSIM
cd /home/action/elecar/PRESTUS

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
%addpath /home/action/elecar/orca-lab/project/tuSIM/toolboxes
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

% Set config files and export location
config_left_transducer = 'sjoerd_config_opt_CTX250-011_64.5mm.yaml';
%config_left_transducer = 'new_sjoerd_config_opt_CTX250-001_105_60.9mm.yaml';
%config_left_transducer = 'sjoerd_config_opt_CTX500-024_72.6mm.yaml';
config_right_transducer = 'sjoerd_config_opt_CTX250-001_64.5mm.yaml';
%config_right_transducer = 'new_sjoerd_config_opt_CTX250-001_105_60.9mm.yaml';
%config_right_transducer = 'sjoerd_config_opt_CTX500-026_73.5mm.yaml';

overwrite_option = 'always'; %change to never if i don't want to overwrite (if i have to stop the simulation
%for example)
parameters = load_parameters(config_left_transducer);
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
%subject_list = [1,3,4,5,8,9,10,14,17,18,19]; % Temporary, selects subjects with complete files
%subject_list = [8,14]; % Thickest skulls
%subject_list = [3, 17]; % Most permissable skulls
%subject_list = [23,24,25,26,27,28,29];
%subject_list = [21,24,26,27,28,31,33,34,36,37,38,39,41]  % New one!   %For 27 check the
%transducers positions!check also 36 (try to invert); 31 has a super high
%Isppa on one side!
%subject_list = [59];
subject_list = [43,45,51];

correctly_named_transducers = [1,5,14,24,28,33,34,36,38,39,41,43,45,47,48,49,50,51,52,53,55,56,57,59,60,61]; %37 should be correctly named!

for subject_id = subject_list

    % Setting folder locations
    subj_folder = fullfile(parameters.data_path,sprintf('sub-%1$03d/', subject_id));
    filename_t1 = dir(sprintf(fullfile(parameters.data_path,parameters.t1_path_template), subject_id));
    t1_header = niftiinfo(fullfile(filename_t1.folder,filename_t1.name));
    t1_image = niftiread(fullfile(filename_t1.folder,filename_t1.name));
    
    % Left transducer localite files
    if ismember(subject_id, correctly_named_transducers)
        %trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_left*.xml',parameters.data_path, subject_id, subject_id));
        trig_mark_files = dir(sprintf('%ssub-%03d/ses-tus01/simu/localite_sub%03d_ses01_left*.xml',parameters.data_path, subject_id, subject_id));
        extract_dt = @(x) datetime(x.name(28:end-4),'InputFormat','yyyyMMddHHmmssSSS');
    else
        trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_right*.xml',parameters.data_path, subject_id, subject_id));
        %trig_mark_files = dir(sprintf('%ssub-%03d/ses-tus01/simu/localite_sub%03d_ses01_left*.xml',parameters.data_path, subject_id, subject_id));
        extract_dt = @(x) datetime(x.name(29:end-4),'InputFormat','yyyyMMddHHmmssSSS'); % adjust it
    end
    
    % sort by datetime
    [~,idx] = sort([arrayfun(extract_dt,trig_mark_files)],'descend');
    trig_mark_files = trig_mark_files(idx);
    
    % Translate transducer trigger markers to raster positions
    [left_trans_ras_pos, left_amygdala_ras_pos] = get_trans_pos_from_trigger_markers(fullfile(trig_mark_files(1).folder, trig_mark_files(1).name), 5, ...
        reference_to_transducer_distance, parameters.expected_focal_distance_mm);
    left_trans_pos = ras_to_grid(left_trans_ras_pos, t1_header);
    %left_trans_pos = [18;163;134]
    left_amygdala_pos = ras_to_grid(left_amygdala_ras_pos, t1_header);
    %left_amygdala_pos = [89;163;134]
    
    % Right transducer localite file
    if ismember(subject_id, correctly_named_transducers)
        %trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_right*.xml',parameters.data_path, subject_id, subject_id));
        trig_mark_files = dir(sprintf('%ssub-%03d/ses-tus01/simu/localite_sub%03d_ses01_right*.xml',parameters.data_path, subject_id, subject_id));
        extract_dt = @(x) datetime(x.name(29:end-4),'InputFormat','yyyyMMddHHmmssSSS');
    else
        trig_mark_files = dir(sprintf('%ssub-%03d/localite_sub%03d_ses01_left*.xml',parameters.data_path, subject_id, subject_id));
        %trig_mark_files = dir(sprintf('%ssub-%03d/ses-tus01/simu/localite_sub%03d_ses01_right*.xml',parameters.data_path, subject_id, subject_id));
        extract_dt = @(x) datetime(x.name(28:end-4),'InputFormat','yyyyMMddHHmmssSSS'); % adjust it
    end
    
    % sort by datetime
    [~,idx] = sort([arrayfun(extract_dt,trig_mark_files)],'descend');
    trig_mark_files = trig_mark_files(idx);
    
    % Translate transducer trigger markers to raster positions
    [right_trans_ras_pos, right_amygdala_ras_pos] = get_trans_pos_from_trigger_markers(fullfile(trig_mark_files(1).folder, trig_mark_files(1).name), 5, ...
        reference_to_transducer_distance, parameters.expected_focal_distance_mm);
    right_trans_pos = ras_to_grid(right_trans_ras_pos, t1_header);
    %right_trans_pos = [206;163;134]
    right_amygdala_pos = ras_to_grid(right_amygdala_ras_pos, t1_header);
    %right_amygdala_pos = [135;163;134]

    % Generate plots with both transducers and targets separately
    imshowpair(plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), left_trans_pos, left_amygdala_pos, parameters), plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), right_trans_pos, right_amygdala_pos, parameters),'montage');
    
    % Index transducer locations for simulation selection (and flip if necessary)
    transducers = [left_trans_pos right_trans_pos];
    targets = [left_amygdala_pos right_amygdala_pos];
    target_names = {'left_amygdala', 'right_amygdala'};
    
    % Simulations for left amygdala
    % Loading parameters
    parameters = load_parameters(config_left_transducer);
    parameters.overwrite_files = overwrite_option;
    
    % Select 'layered' when simulating the transmission in a skull
    parameters.simulation_medium = 'layered';
    reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);
    
    target_id = 1;
    parameters.transducer.pos_t1_grid = transducers(:,target_id)';
    parameters.focus_pos_t1_grid = targets(:,target_id)';
    parameters.results_filename_affix = sprintf('_target_%s', target_names{target_id});
    parameters.interactive = 0;
    %single_subject_pipeline(subject_id, parameters);
    single_subject_pipeline_with_qsub(subject_id, parameters);
    %qsubfeval(@single_subject_pipeline_wrapper, subject_id, parameters, 'timreq',  60604,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
    
    % Simulations for right amygdala
    parameters = load_parameters(config_right_transducer);
    parameters.overwrite_files = overwrite_option;
    
    % Select 'layered' when simulating the transmission in a skull
    parameters.simulation_medium = 'layered';
    reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);
    
    target_id = 2;
    parameters.transducer.pos_t1_grid = transducers(:,target_id)';
    parameters.focus_pos_t1_grid = targets(:,target_id)';
    parameters.results_filename_affix = sprintf('_target_%s', target_names{target_id});
    parameters.interactive = 0;
    %single_subject_pipeline(subject_id, parameters);
    single_subject_pipeline_with_qsub(subject_id, parameters);
    %qsubfeval(@single_subject_pipeline_wrapper, subject_id, parameters, 'timreq',  60604,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');

end