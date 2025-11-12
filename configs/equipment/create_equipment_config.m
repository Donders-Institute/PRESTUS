clear; 

config_path = fileparts(mfilename('fullpath'));
main_folder = fileparts(fileparts(config_path));

% add yaml toolbox to path
addpath(fullfile(main_folder,'toolboxes','yaml'));

% go to correct directory to save config in
cd(config_path)

% Initialize the structure
config = struct();
config.gen.charac_folder = '';
config.gen.axial_prof_name = 'Axial_profiles_'; % First part of axial profiles filename, second part is always a combination of [transducer name]~[driving system name]
config.gen.prestus_virt_folder = 'PRESTUS_virtual_parameters/';
config.gen.prestus_virt_name = 'PRESTUS_virtual_'; % First part of PRESTUS virtual parameters filename, second part is always a combination of [transducer name]~[driving system name]
config.gen.phase_table_folder = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Transducers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                   Sonic Concepts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SONIC_CONCEPTS = 'Sonic Concepts';

field_name_tran_1 = 'CTX_250_001';
config.trans.(field_name_tran_1) = struct();
config.trans.(field_name_tran_1).serial = field_name_tran_1;
config.trans.(field_name_tran_1).name = string(strcat('NeuroFUS 4 ch.', {' '}, config.trans.(field_name_tran_1).serial));
config.trans.(field_name_tran_1).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_1).n_elem = 4; % number of elements
config.trans.(field_name_tran_1).min_foc = 13.79; % [mm], focus wrt exit plane
config.trans.(field_name_tran_1).max_foc = 61.48; % [mm], focus wrt exit plane
config.trans.(field_name_tran_1).prestus.transducer.n_elements = 4; % number of virtual elements
config.trans.(field_name_tran_1).prestus.transducer.Elements_ID_mm = [0, 32.9184, 46.1264, 56.0324]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_1).prestus.transducer.Elements_OD_mm = [32.3596, 45.5676, 55.5244, 64.008]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_1).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_1).prestus.transducer.dist_to_plane_mm = 52.38; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_1).prestus.transducer.source_freq_hz = 250000; % [Hz]

field_name_tran_2 = 'CTX_250_009';
config.trans.(field_name_tran_2).serial = field_name_tran_2;
config.trans.(field_name_tran_2).name = string(strcat('NeuroFUS 2 ch.', {' '}, config.trans.(field_name_tran_2).serial));
config.trans.(field_name_tran_2).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_2).n_elem = 2; % number of elements
config.trans.(field_name_tran_2).min_foc = 15.9; % [mm], focus wrt exit plane
config.trans.(field_name_tran_2).max_foc = 46.0; % [mm], focus wrt exit plane
config.trans.(field_name_tran_2).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_2).prestus.transducer.Elements_ID_mm = [0    6.5764   13.1529   19.7293   26.3058, 32.6136   35.2535   37.8935   40.5334   43.1734]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_2).prestus.transducer.Elements_OD_mm = [5.7744   12.3509   18.9273   25.5038   32.0802, 34.9316   37.5716   40.2115   42.8515   45.4914]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_2).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_2).prestus.transducer.dist_to_plane_mm = 56.87;
config.trans.(field_name_tran_2).prestus.transducer.source_freq_hz = 250000; % [Hz]

field_name_tran_3 = 'CTX_250_014';
config.trans.(field_name_tran_3).serial = field_name_tran_3;
config.trans.(field_name_tran_3).name = string(strcat('NeuroFUS 2 ch.', {' '}, config.trans.(field_name_tran_3).serial));
config.trans.(field_name_tran_3).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_3).n_elem = 2; % number of elements
config.trans.(field_name_tran_3).min_foc = 12.6; % [mm], focus wrt exit plane
config.trans.(field_name_tran_3).max_foc = 44.1; % [mm], focus wrt exit plane
config.trans.(field_name_tran_3).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_3).prestus.transducer.Elements_ID_mm = [0    6.5764   13.1529   19.7293   26.3058, 32.6136   35.2535   37.8935   40.5334   43.1734]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_3).prestus.transducer.Elements_OD_mm = [5.7744   12.3509   18.9273   25.5038   32.0802, 34.9316   37.5716   40.2115   42.8515   45.4914]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_3).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_3).prestus.transducer.dist_to_plane_mm = 56.87;
config.trans.(field_name_tran_3).prestus.transducer.source_freq_hz = 250000; % [Hz]

field_name_tran_4 = 'CTX_250_026';
config.trans.(field_name_tran_4).serial = field_name_tran_4;
config.trans.(field_name_tran_4).name = string(strcat('NeuroFUS 4 ch.', {' '}, config.trans.(field_name_tran_4).serial));
config.trans.(field_name_tran_4).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_4).n_elem = 4; % number of elements
config.trans.(field_name_tran_4).min_foc = 21.9; % [mm], focus wrt exit plane
config.trans.(field_name_tran_4).max_foc = 61.5; % [mm], focus wrt exit plane
config.trans.(field_name_tran_4).prestus.transducer.n_elements = 4; % number of virtual elements
config.trans.(field_name_tran_4).prestus.transducer.Elements_ID_mm = [0, 32.9184, 46.1264, 56.0324]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_4).prestus.transducer.Elements_OD_mm = [32.3596, 45.5676, 55.5244, 64.008]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_4).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_4).prestus.transducer.dist_to_plane_mm = 52.38; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_4).prestus.transducer.source_freq_hz = 250000; % [Hz]

field_name_tran_5 = 'CTX_500_006';
config.trans.(field_name_tran_5).serial = field_name_tran_5;
config.trans.(field_name_tran_5).name = string(strcat('NeuroFUS 2 ch.', {' '}, config.trans.(field_name_tran_5).serial));
config.trans.(field_name_tran_5).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_5).n_elem = 2; % number of elements
config.trans.(field_name_tran_5).min_foc = 33.2; % [mm], focus wrt exit plane
config.trans.(field_name_tran_5).max_foc = 79.4; % [mm], focus wrt exit plane
config.trans.(field_name_tran_5).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_5).prestus.transducer.Elements_ID_mm = [0, 4.59969, 9.19937, 13.7991, 18.3987, 22.9984, 27.5981, 32.1978, 36.7975, 41.3972]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_5).prestus.transducer.Elements_OD_mm = [4.09423, 8.69391, 13.2936, 17.8933, 22.493, 27.0927, 31.6923, 36.292, 40.8917, 45.4914]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_5).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_5).prestus.transducer.dist_to_plane_mm = 56.87; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_5).prestus.transducer.source_freq_hz = 500000; % [Hz]

field_name_tran_6 = 'CTX_500_026';
config.trans.(field_name_tran_6).serial = field_name_tran_6;
config.trans.(field_name_tran_6).name = string(strcat('NeuroFUS 4 ch.', {' '}, config.trans.(field_name_tran_6).serial));
config.trans.(field_name_tran_6).manufact = SONIC_CONCEPTS;
config.trans.(field_name_tran_6).n_elem = 4; % number of elements
config.trans.(field_name_tran_6).min_foc = 39.6; % [mm], focus wrt exit plane
config.trans.(field_name_tran_6).max_foc = 79.6; % [mm], focus wrt exit plane
config.trans.(field_name_tran_6).prestus.transducer.n_elements = 4; % number of virtual elements
config.trans.(field_name_tran_6).prestus.transducer.Elements_ID_mm = [0, 33.02, 46.228, 56.0832]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_6).prestus.transducer.Elements_OD_mm = [32.512, 45.7708, 55.626, 64.008]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_6).prestus.transducer.curv_radius_mm = 63.20; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_6).prestus.transducer.dist_to_plane_mm = 52.38; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_6).prestus.transducer.source_freq_hz = 500000; % [Hz]

%                   Imasonic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IMASONIC = 'Imasonic';

field_name_tran_7 = 'IS_PCD15287_01001';
config.trans.(field_name_tran_7).serial = field_name_tran_7;
config.trans.(field_name_tran_7).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_7).serial));
config.trans.(field_name_tran_7).manufact = IMASONIC;
config.trans.(field_name_tran_7).n_elem = 10; % number of elements
config.trans.(field_name_tran_7).min_foc = 5.0; % [mm], focus wrt exit plane
config.trans.(field_name_tran_7).max_foc = 91.7; % [mm], focus wrt exit plane
config.trans.(field_name_tran_7).prestus.transducer.n_elements = 30; % number of virtual elements
config.trans.(field_name_tran_7).prestus.transducer.Elements_ID_mm = [10, 13.7667, 17.5333, 22.3, 24.5667, 26.8333, 30, 31.7667, 33.5333, 36.3, 37.7667, 39.2333, 41.7, 43.0, 44.3, 46.5, 47.6667, 48.8333, 51, 52.0333, 53.0667, 55.1, 56.0667, 57.0333, 58.9, 59.8, 60.7, 62.5, 63.3333, 64.1667]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_7).prestus.transducer.Elements_OD_mm = [13.7667, 17.5333, 21.3, 24.5667, 26.8333, 29.1, 31.7667, 33.5333, 35.3, 37.7667, 39.2333, 40.7, 43.0, 44.3, 45.6, 47.6667, 48.8333, 50, 52.0333, 53.0667, 54.1, 56.0667, 57.0333, 58, 59.8, 60.7, 61.6, 63.3333, 64.1667, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_7).prestus.transducer.curv_radius_mm = 75; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_7).prestus.transducer.dist_to_plane_mm = 65.3; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_7).prestus.transducer.source_freq_hz = 300000; % [Hz]

field_name_tran_8 = 'IS_PCD15287_01002';
config.trans.(field_name_tran_8).serial = field_name_tran_8;
config.trans.(field_name_tran_8).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_8).serial));
config.trans.(field_name_tran_8).manufact = IMASONIC;
config.trans.(field_name_tran_8).n_elem = 10; % number of elements
config.trans.(field_name_tran_8).min_foc = 6.1; % [mm], focus wrt exit plane
config.trans.(field_name_tran_8).max_foc = 93.2; % [mm], focus wrt exit plane
config.trans.(field_name_tran_8).prestus.transducer.n_elements = 30; % number of virtual elements
config.trans.(field_name_tran_8).prestus.transducer.Elements_ID_mm = [10, 13.7667, 17.5333, 22.3, 24.5667, 26.8333, 30, 31.7667, 33.5333, 36.3, 37.7667, 39.2333, 41.7, 43.0, 44.3, 46.5, 47.6667, 48.8333, 51, 52.0333, 53.0667, 55.1, 56.0667, 57.0333, 58.9, 59.8, 60.7, 62.5, 63.3333, 64.1667]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_8).prestus.transducer.Elements_OD_mm = [13.7667, 17.5333, 21.3, 24.5667, 26.8333, 29.1, 31.7667, 33.5333, 35.3, 37.7667, 39.2333, 40.7, 43.0, 44.3, 45.6, 47.6667, 48.8333, 50, 52.0333, 53.0667, 54.1, 56.0667, 57.0333, 58, 59.8, 60.7, 61.6, 63.3333, 64.1667, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_8).prestus.transducer.curv_radius_mm = 75; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_8).prestus.transducer.dist_to_plane_mm = 65.3; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_8).prestus.transducer.source_freq_hz = 300000; % [Hz]

field_name_tran_9 = 'IS_PCD15473_01001';
config.trans.(field_name_tran_9).serial = field_name_tran_9;
config.trans.(field_name_tran_9).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_9).serial));
config.trans.(field_name_tran_9).manufact = IMASONIC;
config.trans.(field_name_tran_9).n_elem = 10; % number of elements
config.trans.(field_name_tran_9).min_foc = 6.7; % [mm], focus wrt exit plane
config.trans.(field_name_tran_9).max_foc = 92.6; % [mm], focus wrt exit plane
config.trans.(field_name_tran_9).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_9).prestus.transducer.Elements_ID_mm = [10, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_9).prestus.transducer.Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_9).prestus.transducer.curv_radius_mm = 100; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_9).prestus.transducer.dist_to_plane_mm = 92.7; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_9).prestus.transducer.source_freq_hz = 300000; % [Hz]

field_name_tran_10 = 'IS_PCD15473_01003';
config.trans.(field_name_tran_10).serial = field_name_tran_10;
config.trans.(field_name_tran_10).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_10).serial));
config.trans.(field_name_tran_10).manufact = IMASONIC;
config.trans.(field_name_tran_10).n_elem = 10; % number of elements
config.trans.(field_name_tran_10).min_foc = 7.2; % [mm], focus wrt exit plane
config.trans.(field_name_tran_10).max_foc = 93.6; % [mm], focus wrt exit plane
config.trans.(field_name_tran_10).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_10).prestus.transducer.Elements_ID_mm = [10, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_10).prestus.transducer.Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_10).prestus.transducer.curv_radius_mm = 100; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_10).prestus.transducer.dist_to_plane_mm = 92.7; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_10).prestus.transducer.source_freq_hz = 300000; % [Hz]

field_name_tran_11 = 'IS_PCD15473_01001_OPM';
config.trans.(field_name_tran_11).serial = field_name_tran_11;
config.trans.(field_name_tran_11).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_11).serial));
config.trans.(field_name_tran_11).manufact = IMASONIC;
config.trans.(field_name_tran_11).n_elem = 10; % number of elements
config.trans.(field_name_tran_11).min_foc = 7.3; % [mm], focus wrt exit plane
config.trans.(field_name_tran_11).max_foc = 92.0; % [mm], focus wrt exit plane
config.trans.(field_name_tran_11).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_11).prestus.transducer.Elements_ID_mm = [10, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_11).prestus.transducer.Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_11).prestus.transducer.curv_radius_mm = 100; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_11).prestus.transducer.dist_to_plane_mm = 92.7; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_11).prestus.transducer.source_freq_hz = 300000; % [Hz]

field_name_tran_12 = 'IS_PCD15473_01003_OPM';
config.trans.(field_name_tran_12).serial = field_name_tran_12;
config.trans.(field_name_tran_12).name = string(strcat(IMASONIC, {' '}, ' 10 ch.', {' '}, config.trans.(field_name_tran_12).serial));
config.trans.(field_name_tran_12).manufact = IMASONIC;
config.trans.(field_name_tran_12).n_elem = 10; % number of elements
config.trans.(field_name_tran_12).min_foc = 6.7; % [mm], focus wrt exit plane
config.trans.(field_name_tran_12).max_foc = 93.2; % [mm], focus wrt exit plane
config.trans.(field_name_tran_12).prestus.transducer.n_elements = 10; % number of virtual elements
config.trans.(field_name_tran_12).prestus.transducer.Elements_ID_mm = [10, 22.1, 29.8, 36.0, 41.4, 46.3, 50.7, 54.9, 58.7, 62.4]; % Inner diameter of each element [mm]
config.trans.(field_name_tran_12).prestus.transducer.Elements_OD_mm = [21.1, 28.8, 35.0, 40.4, 45.3, 49.7, 53.9, 57.8, 61.5, 65]; % Outer diameter of each element [mm]
config.trans.(field_name_tran_12).prestus.transducer.curv_radius_mm = 100; % Radius of curvature of the bowl  [mm]
config.trans.(field_name_tran_12).prestus.transducer.dist_to_plane_mm = 92.7; % Distance to the transducer exit plane from the geometric focus [mm]
config.trans.(field_name_tran_12).prestus.transducer.source_freq_hz = 300000; % [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Driving systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                   Sonic Concepts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_name_ds_1 = 'SC_203_035';
config.ds.(field_name_ds_1).serial = field_name_ds_1;
config.ds.(field_name_ds_1).name = string(strcat('NeuroFUS 1 x 4 ch. or 1 x 2 ch. TPO junior', {' '}, config.ds.(field_name_ds_1).serial));
config.ds.(field_name_ds_1).manufact = SONIC_CONCEPTS;
config.ds.(field_name_ds_1).available_ch = 4;

field_name_ds_2 = 'SC_105_010';
config.ds.(field_name_ds_2).serial = field_name_ds_2;
config.ds.(field_name_ds_2).name = string(strcat('NeuroFUS 1 x 4 ch. or 1 x 2 ch. TPO senior', {' '}, config.ds.(field_name_ds_2).serial));
config.ds.(field_name_ds_2).manufact = SONIC_CONCEPTS;
config.ds.(field_name_ds_2).manufact = SONIC_CONCEPTS;
config.ds.(field_name_ds_2).available_ch = 4;

%                   IGT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IGT = 'IGT';

field_name_ds_3 = 'IGT_128_ch_comb_10_ch';
config.ds.(field_name_ds_3).serial = field_name_ds_3;
config.ds.(field_name_ds_3).name = string(strcat(IGT, {' '}, ' 128 ch. - 10 ch.'));
config.ds.(field_name_ds_3).manufact = IGT;
config.ds.(field_name_ds_3).available_ch = 128;

field_name_ds_4 = 'IGT_32_ch_comb_10_ch';
config.ds.(field_name_ds_4).serial = field_name_ds_4;
config.ds.(field_name_ds_4).name = string(strcat(IGT, {' '}, ' 32 ch. - 10 ch.'));
config.ds.(field_name_ds_4).manufact = IGT;
config.ds.(field_name_ds_4).available_ch = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Combinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

split_combo = string(strcat(field_name_tran_1, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_1;
config.combos.(split_combo).ds_serial = field_name_ds_1; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-001~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-001 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_2, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_2;
config.combos.(split_combo).ds_serial = field_name_ds_1; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-009~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-009 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_3, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_3;
config.combos.(split_combo).ds_serial = field_name_ds_1; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-014~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-014 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_4, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_4;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-026~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-026 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_5, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_5;
config.combos.(split_combo).ds_serial = field_name_ds_1; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-500-006~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-500-006 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_6, '_', field_name_ds_1)); 
config.combos.(split_combo).tran_serial = field_name_tran_6;
config.combos.(split_combo).ds_serial = field_name_ds_1; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-500-026~TPO-203-035.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-500-026 - TPO-203-035 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_1, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_1;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name, 'CTX-250-001~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-001 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_2, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_2;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name, 'CTX-250-009~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-009 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_3, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_3;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-014~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-014 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_4, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_4;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-250-026~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-250-026 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_5, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_5;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-500-006~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-500-006 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_6, '_', field_name_ds_2)); 
config.combos.(split_combo).tran_serial = field_name_tran_6;
config.combos.(split_combo).ds_serial = field_name_ds_2; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'CTX-500-026~TPO-105-010.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'CTX-500-026 - TPO-105-010 - Steer Table.csv');

split_combo = string(strcat(field_name_tran_7, '_', field_name_ds_3)); 
config.combos.(split_combo).tran_serial = field_name_tran_7;
config.combos.(split_combo).ds_serial = field_name_ds_3; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15287_01001~IGT_128_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15287_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_8, '_', field_name_ds_3)); 
config.combos.(split_combo).tran_serial = field_name_tran_8;
config.combos.(split_combo).ds_serial = field_name_ds_3; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15287_01002~IGT_128_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15287_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_9, '_', field_name_ds_3)); 
config.combos.(split_combo).tran_serial = field_name_tran_9;
config.combos.(split_combo).ds_serial = field_name_ds_3; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01001~IGT_128_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_10, '_', field_name_ds_3)); 
config.combos.(split_combo).tran_serial = field_name_tran_10;
config.combos.(split_combo).ds_serial = field_name_ds_3; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01003~IGT_128_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_7, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_7;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15287_01001~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15287_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_8, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_8;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15287_01002~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15287_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_9, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_9;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01001~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_10, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_10;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01003~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_11, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_11;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01001_OPM~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

split_combo = string(strcat(field_name_tran_12, '_', field_name_ds_4)); 
config.combos.(split_combo).tran_serial = field_name_tran_12;
config.combos.(split_combo).ds_serial = field_name_ds_4; 
config.combos.(split_combo).char_data_path = strcat(config.gen.charac_folder, config.gen.axial_prof_name,  'IS_PCD15473_01003_OPM~IGT_32_ch_comb_10_ch.csv');
config.combos.(split_combo).phase_table = strcat(config.gen.phase_table_folder,  'transducer_15473_10_300kHz.ini');

yaml.dumpFile('equipment_config.yaml', config)
