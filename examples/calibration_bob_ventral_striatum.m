clear
close all
cd /project/3015999.02/andche_sandbox/orca-lab/project/tuSIM/ % change path here

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

%parameters = load_parameters('bob_config_opt_CTX250-011_64.5mm.yaml');
parameters = load_parameters('_config_opt_CTX250-011_64.5mm.yaml');
transducer_thickness = parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm; % from the transducer specifications

%% load the real acoustic profile
real_profile = readtable('examples/CTX_250-011_4chan_all_distances.csv');
%real_profile(:, end) = [] ;
steering_opts = table2array(real_profile(1,2:end));
real_profile.Properties.VariableNames = cellfun(@string, table2cell(real_profile(1,:)));
real_profile.Properties.VariableNames(1) = "Distance";
real_profile.Distance = str2double(real_profile.Distance);
real_profile = real_profile(2:(end-1),:);

%real_profile = real_profile(:, [1, end]);
correction = readtable('examples/CTX-250-011_steerTable_2022-08-23_13-44-21.csv');
correction = table2array(correction(:,["Distance","ISPPAFactor"]));

real_profile = table2array(real_profile);
real_profile = real_profile(1:end-5,:);
real_profile(:,2:end) =  real_profile(:,2:end) .* bounded_interp1(correction(:,1),correction(:,2),real_profile(:,1));
real_profile_norm = real_profile;

real_profile_norm(:,2:end) = real_profile_norm(:,2:end) ./ max(real_profile_norm(:,2:end), [],1) * 30;
real_profile_2d = reshape(real_profile_norm(:,2:end),[], 1);
dist_to_peak = arrayfun(@(i) real_profile_norm(real_profile_norm(:,i)==max(real_profile_norm(:,i)),1), 2:size(real_profile_norm, 2));
real_profile_2d = [repmat(real_profile_norm(:,1), [size(real_profile_norm,2)-1, 1]), repelem(steering_opts , size(real_profile_norm, 1))', real_profile_2d ];

x = real_profile_2d(:,1);
y = real_profile_2d(:,2);
z = real_profile_2d(:,3);
xlin = linspace(min(x), 150, 300);
ylin = linspace(min(y), 100 + 30, 300);
[X,Y] = meshgrid(xlin, ylin);
Z = griddata(x,y,z,X,Y,'natural');

xlin_orig = linspace(min(x), max(x), 300);
ylin_orig = linspace(min(y), max(y), 300);
[X_orig_lim,Y_orig_lim] = meshgrid(xlin_orig, ylin_orig);
intrapolated_acoustic_profile = griddata(x,y,z,X_orig_lim,Y_orig_lim,'natural');
%Z = griddata(x,y,z,X,Y,'cubic');
%Z = griddata(x,y,z,X,Y,'v4');

figure('Position', [200 200 900 300]);
subplot(1,3,1)
surf(X,Y,Z,'EdgeColor', 'none')
axis tight; 
xlim([min(xlin), max(xlin)])
ylim([min(ylin), max(ylin)])
view(0,90);
%plot3(x,y,z,'.','MarkerSize',15);
xlabel('Distance from exit plane [mm]');
ylabel('Steering setting (peak center) [mm]');
zlabel('ISPPA [w/cm^2]');

% 
% ft = fittype( 'thinplateinterp' );
opts = fitoptions( 'lowess' );
opts.Normalize = 'off';
% opts.Robust = 'Bisquare';
opts.Span = 0.3;

f = fit([x,y],z, 'lowess',opts);


extrapolated_acoustic_profile = reshape(feval(f,[X(:) Y(:)]), size(X));
extrapolated_acoustic_profile(extrapolated_acoustic_profile<0) = 0;
extrapolated_acoustic_profile = extrapolated_acoustic_profile./max(extrapolated_acoustic_profile, [],2) * 30;
extrapolated_acoustic_profile = imgaussfilt(extrapolated_acoustic_profile, 3);
subplot(1,3,2)
surf(X,Y,extrapolated_acoustic_profile,'EdgeColor', 'none')

xlabel('Distance from exit plane [mm]');
ylabel('Steering setting to peak (peak center) [mm]');
zlabel('ISPPA [w/cm^2]');
view(0,90);

surf(xlin, ylin, extrapolated_acoustic_profile ,'EdgeColor', 'none')
xlim([min(xlin), max(xlin)])
ylim([min(ylin), max(ylin)])
xlabel('Distance from exit plane [mm]');
ylabel('Steering setting to peak (peak center) [mm]');
zlabel('ISPPA [w/cm^2]');

view(0,90);
axis tight; 
subplot(1,3,3)

[~, closest_measured] = min(abs(ylin-max(y)));
[~, closest_measured_to_min] = min(abs(ylin-min(y)));
[~, closest_measured_orig_lim] = min(abs(ylin_orig-max(y)));
[~, closest_measured_to_min_orig_lim] = min(abs(ylin_orig-min(y)));

plot(xlin, extrapolated_acoustic_profile(end,:),'r--',  x(y==min(y)),  z(y==min(y)),'g-',  x(y==max(y)),z(y==max(y)),'b-',...
    xlin, extrapolated_acoustic_profile(closest_measured_to_min,:),'g--',xlin, extrapolated_acoustic_profile(closest_measured,:),'b--',...
    xlin_orig, intrapolated_acoustic_profile(closest_measured_to_min_orig_lim,:),'g-.',xlin_orig, intrapolated_acoustic_profile(closest_measured_orig_lim,:),'b-.')
xlabel('Distance to peak [mm]');
ylabel('ISPPA [w/cm^2]');


%[ylin(1:58); arrayfun(@(i) xlin(Z(i,:)==max(Z(i,:))), 1:58)]'

%%

out_folder = parameters.data_path+'/sim_outputs/';
mni_targets = struct('amygdala',[-26 -4 -20], 'ventral_striatum',[-8 6 -10], 'right_mediodorsal_thalamus', [6.5 -18 9], 'right_pulvinar', [15 -27 6.5] );

good_positions = readtable(out_folder + '/unique_targets.csv')
good_positions = good_positions(good_positions.realistic_position == 2,:)
good_positions.full_idx = strrep(good_positions.all,'_positioning_bob','')
all_targets = fieldnames(mni_targets);
all_subjs = [3,4,8,9,10];
%all_subjs = [4]
%all_targets = all_targets(4)

for subject_id = all_subjs

    headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
    filename_segmented_headreco = fullfile(headreco_folder, sprintf('sub-%03d_masks_contr.nii.gz', subject_id));

    segmented_img_orig = niftiread(filename_segmented_headreco);
    segmented_img_head = niftiinfo(filename_segmented_headreco);
    pixel_size = segmented_img_head.PixelDimensions(1);

    for target_idx = 1:length(all_targets)
        target = all_targets{target_idx};
        tpos_pars = readtable(sprintf('%stpars_subj%03i_%s.csv', parameters.data_path, subject_id, target),'Delimiter',',');

        tpos_pars.dist_from_real_exit_plane_mm = tpos_pars.dist_to_target*pixel_size - (parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm);

        tppf = tpos_pars(tpos_pars.prop_intersect < 0.05,:); % less than 5% intersection with the skin
        tppf = tppf(tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.3),:);
        [~, sort_i] = sort(tppf.var_dist_skull);
        tppf = tppf(sort_i, :);
        tppf_filtered = tppf(1,:);
        for j = 1:10

            dist = pdist2(table2array(tppf_filtered(:,["trans_x","trans_y","trans_z"])),table2array(tppf(:,["trans_x","trans_y","trans_z"])));
            tppf_to_add = tppf(min(dist,[],1)>(10/pixel_size),:);
            if ~isempty(tppf_to_add)
                tppf_filtered = [tppf_filtered;tppf_to_add(1,:)];
            else
                break
            end
        end
        % 
        output_dir = fullfile(parameters.data_path, 'sim_outputs');
        % 
        if (subject_id == 4)
            existing_files = dir([output_dir + sprintf('/sub-%03d*layered*bob_%s*mat',subject_id, target)]);
            if length(existing_files)>0

                idx = arrayfun(@(i) str2num(regexp(existing_files(i).name,'\d+(?=\.)','match','once')), 1:length(existing_files));
                tppf_filtered = tpos_pars(ismember(tpos_pars.idx, union(idx, [tppf_filtered.idx])),:);
            end
        end
        
        tppf_filtered.full_idx = arrayfun(@(i) sprintf('%03d_%s_%i',subject_id,target,tppf_filtered.idx(i)), 1:length(tppf_filtered.idx), 'UniformOutput', false)';
        %tppf_filtered = tppf_filtered(tppf_filtered.idx == 12592136,:)
        tppf_filtered = tppf_filtered(ismember(tppf_filtered.full_idx, good_positions.full_idx),:)
        do_positioning_plots = 0;

        for cur_trans_pos_idx = 1:size(tppf_filtered,1)
            %fprintf('Transducer position: %i\n', cur_trans_pos_idx )
        close all

        best_trans_pos = tppf_filtered(cur_trans_pos_idx,:)

        parameters.transducer.pos_t1_grid = round(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
        parameters.focus_pos_t1_grid = round(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));
        parameters.expected_focal_distance_mm = best_trans_pos.dist_to_target*pixel_size;
        parameters.results_filename_affix = sprintf('_bob_%s_%i',target, best_trans_pos.idx);

        if do_positioning_plots
            output_plot = fullfile(output_dir,sprintf('sub-%03d_positioning%s.png', subject_id,  parameters.results_filename_affix));
            % if exist(output_plot,'file')
            %     continue
            % end

            size_diff_trans_img = parameters.transducer.pos_t1_grid - size(segmented_img_orig);
            size_diff_trans_img(size_diff_trans_img<0) = 0;
            segmented_img_orig_new = segmented_img_orig;

            if any(size_diff_trans_img)
                segmented_img_orig_new = padarray(segmented_img_orig, size_diff_trans_img,'post');
            end
            coord_mesh_xyz = gpuArray(get_xyz_mesh(segmented_img_orig_new));

            [rotated_img, trans_xyz, target_xyz, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
                align_to_focus_axis_and_scale(segmented_img_orig_new, segmented_img_head, parameters.transducer.pos_t1_grid', parameters.focus_pos_t1_grid', 1, parameters);
            %%

            view_pos = [0,0];
            slice_cap = [-1,0,0];
            if parameters.transducer.pos_t1_grid(1) > size(segmented_img_orig_new,1)/2
                view_pos = [-180, 0];
                slice_cap = [1,0,0];
            end

            figure('Position', [200 200 1000 300]);

            subplot(1,3,1)
            show_3d_scalp(segmented_img_orig_new, parameters.focus_pos_t1_grid, parameters.transducer.pos_t1_grid, parameters, segmented_img_head.PixelDimensions(1), coord_mesh_xyz, [0 0 0], view_pos, 0)

            subplot(1,3,2)
            show_3d_scalp(segmented_img_orig_new, parameters.focus_pos_t1_grid, parameters.transducer.pos_t1_grid, parameters, pixel_size, coord_mesh_xyz, slice_cap, view_pos, 0)

            TF = maketform('affine', transformation_matrix);

            ax3 = subplot(1,3,3)

            imagesc(squeeze(rotated_img(:,round(trans_xyz(2)),:)))
            rectangle('Position',[target_xyz([3,1])' - 2, 4 4],...
                      'Curvature',[0,0], 'EdgeColor','r',...
                     'LineWidth',2,'LineStyle','-');

            rectangle('Position',[trans_xyz([3,1])' - 2, 4 4],...
                      'Curvature',[0,0], 'EdgeColor','b',...
                     'LineWidth',2,'LineStyle','-');

            line([trans_xyz(3) target_xyz(3)], [trans_xyz(1) target_xyz(1)], 'Color', 'white')
            get_transducer_box(trans_xyz([1,3]), target_xyz([1,3]), segmented_img_head.PixelDimensions(1), parameters)
            colormap(ax3, [0.3 0.3 0.3; lines(12)])

            export_fig(output_plot, '-native')
        end
        
        desired_intensity = 30;
        expected_focus = parameters.expected_focal_distance_mm;
        expected_focus_rel_to_exit_plane = parameters.expected_focal_distance_mm - transducer_thickness;

        if expected_focus_rel_to_exit_plane > max(dist_to_peak)
            focus_pos = [ylin; arrayfun(@(i) xlin(extrapolated_acoustic_profile(i,:)==max(extrapolated_acoustic_profile(i,:))), 1:size(extrapolated_acoustic_profile, 1))]';

            [~, closest_to_exp_focus] = min(abs(focus_pos(:,2)-expected_focus_rel_to_exit_plane));
            exp_acoustic_profile = [xlin; extrapolated_acoustic_profile(closest_to_exp_focus,:)]';
        else
            focus_pos = [ylin_orig; arrayfun(@(i) xlin_orig(intrapolated_acoustic_profile(i,:)==max(intrapolated_acoustic_profile(i,:))), 1:size(intrapolated_acoustic_profile, 1))]';

            [~, closest_to_exp_focus] = min(abs(focus_pos(:,2)-expected_focus_rel_to_exit_plane));
            exp_acoustic_profile = [xlin; intrapolated_acoustic_profile(closest_to_exp_focus,:)]';
        end
        exp_acoustic_profile(:,1) = exp_acoustic_profile(:,1) + transducer_thickness;

        figure('Position', [10 10 900 500]);
        plot(exp_acoustic_profile(:,1),exp_acoustic_profile(:,2))

        halfMax = (min(exp_acoustic_profile(:,2)) + max(exp_acoustic_profile(:,2))) / 2;
        % Find where the data first drops below half the max.
        index1 = find(exp_acoustic_profile(:,2) >= halfMax, 1, 'first');
        % Find where the data last rises above half the max.
        index2 = find(exp_acoustic_profile(:,2) >= halfMax, 1, 'last');

        [~, max_x] =  max(exp_acoustic_profile(:,2));
        flhm = exp_acoustic_profile(index2,1) - exp_acoustic_profile(index1,1);
        flhm_center_x = (exp_acoustic_profile(index2,1) - exp_acoustic_profile(index1,1))/2+exp_acoustic_profile(index1,1);
        [~,min_i] = min(abs(exp_acoustic_profile(:,1)-flhm_center_x));
        flhm_center_intensity = exp_acoustic_profile(min_i, 2);
        xlabel('Axial Position [mm]');
        ylabel('Intensity [W/cm^2]');
        yline(desired_intensity, '--');
        yline(halfMax, '--');
        xline(exp_acoustic_profile(index1,1), '--');
        xline(exp_acoustic_profile(index2,1), '--');
        xline(flhm_center_x,'r--');

        text(flhm_center_x+0.5, flhm_center_intensity+3, sprintf('FLHM center intensity %.2f [W/cm^2] at %.1f mm',flhm_center_intensity,flhm_center_x), "Color",'r');
        [~, closest_point_to_exp_focus] = min(abs(exp_acoustic_profile(:,1)-expected_focus));
        intensity_at_expected_focus = exp_acoustic_profile(closest_point_to_exp_focus,2);
        xline(expected_focus,'b--');
        xline(exp_acoustic_profile(max_x,1),'--','Color','#7E2F8E');
        text(expected_focus+0.5, intensity_at_expected_focus+4, sprintf('Expected focus intensity %.2f [W/cm^2] at %.1f mm (%.1f mm from the exit plane)',intensity_at_expected_focus,expected_focus, expected_focus_rel_to_exit_plane),"Color",'b');
        [max_intensity, max_x] = max(exp_acoustic_profile(:,2));
        text(exp_acoustic_profile(max_x,1)+0.5, max_intensity+1.5, sprintf('Max intensity %.2f [W/cm^2] at %.1f mm',max_intensity,exp_acoustic_profile(max_x,1)),"Color",'#7E2F8E');
        ylim([0,max(exp_acoustic_profile(:,2))+6])

        output_plot = fullfile(output_dir,sprintf('sub-%03d_acoustic_profile%s.png', subject_id,  parameters.results_filename_affix));
        export_fig(output_plot, '-native')
        
        % continue
        %% do water sims with default settings

        parameters.simulation_medium = 'water'; % indicate that we only want the simulation in the water medium for now
        parameters.run_heating_sims = 0;
        % subject_id = 999; % subject id - for free water sims
        % 
        % parameters.interactive = 0;
        % parameters.overwrite_files = 'always';
        % qsubfeval(@single_subject_pipeline_wrapper, subject_id, parameters, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
        % 
        load(sprintf('%s/sim_outputs/sub-%03d_water_results.mat', parameters.data_path , 999),'sensor_data','parameters');

        % get maximum pressure
        p_max = gather(sensor_data.p_max_all); % transform from GPU array to normal array

        parameters.transducer.pos_t1_grid = round(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
        parameters.focus_pos_t1_grid = round(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));
        parameters.expected_focal_distance_mm = best_trans_pos.dist_to_target*pixel_size;
        parameters.results_filename_affix = sprintf('_bob_%s_%i',target, best_trans_pos.idx);

        %% 
        % Then, we can compare the simulated pressure along the focal axis and the pressure 
        % estimated with an analytic solution based on the equations provided by O'Neil 
        % (O'Neil, H. Theory of focusing radiators. J. Acoust. Soc. Am., 21(5), 516-526, 
        % 1949) and implemented in k-wave |focusedAnnulusONeil()| function.

        % simulated pressure along the focal axis
        pred_axial_pressure = squeeze(p_max(parameters.transducer.pos_grid(1),parameters.transducer.pos_grid(2),:)); % get the values at the focal axis

        % compute O'Neil solution and plot it along with comparisons
        % define transducer parameters

        velocity = parameters.transducer.source_amp(1)/(parameters.medium.water.density*parameters.medium.water.sound_speed);   % [m/s]

        % define position vectors
        axial_position_sim_grid = (1:parameters.default_grid_dims(3))*0.5;       % [mm]
        axial_position   = exp_acoustic_profile(:,1)%(1:parameters.default_grid_dims(3))*0.5;       % [mm]

        % evaluate pressure analytically
        % focusedAnnulusONeil provides an analytic solution for the pressure at the
        % focal (beam) axis
        [p_axial_oneil] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
            [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,parameters.transducer.n_elements), ...
            parameters.transducer.source_phase_rad, parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
            parameters.medium.water.density, (axial_position-0.5)*1e-3);

        % plot focal axis pressure
        figure('Position', [10 10 900 500]);

        plot(axial_position, p_axial_oneil .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
        xlabel('Axial Position [mm]');
        ylabel('Intensity [W/cm^2]');
        hold on
        plot(axial_position_sim_grid-(parameters.transducer.pos_grid(3))*parameters.grid_step_mm, ...
            pred_axial_pressure.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4,'--');
        plot(exp_acoustic_profile(:,1),exp_acoustic_profile(:,2))
        hold off
        xline(parameters.expected_focal_distance_mm, '--');
        legend('Analytic solution','Simulated results','Real profile')
        title('Pressure along the beam axis')
        % what is distance to the maximum pressure?
        fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil==max(p_axial_oneil)))

        % compute the approximate adjustment from simulated (on a grid) to analytic solution
        simulated_grid_adj_factor = max(pred_axial_pressure(:))/max(p_axial_oneil(:));

        %% Optimize for a given distance and pressure
        % So how to find the settings for the simulations that match the desired pressure 
        % and distance? It is easy to do, given that there is an analytic solution. For 
        % our subject, we assume that we know where the transducer is positioned and where 
        % we want to have the maximum pressure, so we know the distance at which the pressure 
        % should be maximal. We need to find the set of phases for transducer elements 
        % that would give the maximum pressure at that distance. We do so by searching 
        % through the parameter space (that is, varying the phases) as to minimize the 
        % error in distance to maximum pressure point. 

        gs = GlobalSearch;
        %opt_velocity = desired_pressure/max_pressure*velocity;
        opt_limits = [min(exp_acoustic_profile(:,1)), max(exp_acoustic_profile(:,1))];
        opt_weights = normpdf(exp_acoustic_profile(:,1), expected_focus, 45);%ones([1 length(exp_acoustic_profile(:,1))]);
        %optimize_phases = @(phases) phase_optimization_annulus(phases, parameters, velocity, axial_position, parameters.expected_focal_distance_mm);
        optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(phases_and_velocity(1:(parameters.transducer.n_elements-1)), ...
            parameters, phases_and_velocity(parameters.transducer.n_elements),...
            exp_acoustic_profile(:,1), exp_acoustic_profile(:,2), 0,opt_limits,opt_weights);

        rng(160,'twister') % setting seed for consistency
        problem = createOptimProblem('fmincon','x0', [randi(360, [1 parameters.transducer.n_elements-1])/180*pi rand(1)],...
            'objective',optimize_phases,'lb',zeros(1,4),'ub',[2*pi*ones(1,parameters.transducer.n_elements-1) 0.2], 'options', optimoptions('fmincon','OptimalityTolerance', 1e-8)); 

        [opt_phases_and_velocity, min_err] = run(gs,problem);

        opt_phases = opt_phases_and_velocity(1:(parameters.transducer.n_elements-1));
        opt_velocity = opt_phases_and_velocity(parameters.transducer.n_elements);
        opt_source_amp = round(opt_velocity/velocity*parameters.transducer.source_amp/simulated_grid_adj_factor);

        phase_optimization_annulus_full_curve(opt_phases_and_velocity(1:(parameters.transducer.n_elements-1)), ...
            parameters, opt_phases_and_velocity(parameters.transducer.n_elements),...
            exp_acoustic_profile(:,1), exp_acoustic_profile(:,2), 1, opt_limits, opt_weights)

        fprintf('Optimal phases: %s deg.; velocity: %.2f; optimization error: %.2f', ...
            mat2str(round(opt_phases_and_velocity(1:(parameters.transducer.n_elements-1))/pi*180)), opt_phases_and_velocity(parameters.transducer.n_elements), min_err)

        %%
        [p_axial_oneil_opt] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
            [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(opt_velocity,1,parameters.transducer.n_elements), ...
            [0 opt_phases], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
            parameters.medium.water.density, (axial_position-0.5)*1e-3);

        figure('Position', [10 10 900 500]);
        plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
        hold on
        plot(exp_acoustic_profile(:,1),exp_acoustic_profile(:,2))
        xlabel('Axial Position [mm]');
        ylabel('Intensity [W/cm^2]');
        hold off
        xline(parameters.expected_focal_distance_mm, '--');
        yline(desired_intensity, '--');
        legend( sprintf('Optimized to match the real profile'),'Real profile')
        title('Pressure along the beam axis')

        fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',axial_position(p_axial_oneil_opt==max(p_axial_oneil_opt)))
        %fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', get_flhm_center_position(axial_position, p_axial_oneil_opt))
        parameters.results_filename_affix = sprintf('_bob_%s_%i',target, best_trans_pos.idx);

        output_plot = fullfile(output_dir,sprintf('sub-%03d_fitted_acoustic_profile%s.png', subject_id,  parameters.results_filename_affix));
        export_fig(output_plot, '-native')
        continue


        %% Double-check in the water

        % 
        % opt_parameters = load_parameters('bob_config_opt_CTX250-011_64.5mm.yaml'); 
        % opt_parameters.transducer.source_amp = opt_source_amp;
        % opt_parameters.transducer.source_phase_rad = [0 opt_phases];
        % opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
        % opt_parameters.transducer.pos_t1_grid = round(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
        % opt_parameters.focus_pos_t1_grid = round(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));
        % opt_parameters.expected_focal_distance_mm = best_trans_pos.dist_to_target*pixel_size;
        % opt_parameters.results_filename_affix = sprintf('_bob_%s_%i',target, best_trans_pos.idx);
        % 
        % opt_parameters.simulation_medium = 'water';
        % 
        % opt_parameters.interactive = 0;
        % opt_parameters.overwrite_files = 'never';
        % single_subject_pipeline(subject_id, opt_parameters)
        % 
        % % If you are using the Donders HPC cluster, you can do the simulations in
        % % a non-interactive session with a qsub. To do so, set the interactive flag
        % % to zero and set overwrite_files to 'always' (if you already have the results and want to recompute them). 
        % % % 
        % opt_parameters.interactive = 0;
        % opt_parameters.overwrite_files = 'always';
        % % 
        % qsubfeval(@single_subject_pipeline_wrapper, subject_id, opt_parameters, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
        % 
        %%
        % opt_res = load(sprintf('%s/sim_outputs/sub-%03d_water_results%s.mat', opt_parameters.data_path , subject_id, opt_parameters.results_filename_affix),'sensor_data','parameters');
        % 
        % % get maximum pressure
        % p_max = gather(opt_res.sensor_data.p_max_all);
        % pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1), opt_res.parameters.transducer.pos_grid(2),:));
        % 
        % figure('Position', [10 10 900 500]);
        % hold on
        % plot(axial_position, p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
        % xlabel('Axial Position [mm]');
        % ylabel('Intensity [W/cm^2]');
        % plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
        % 
        % sim_res_axial_position = axial_position-(opt_res.parameters.transducer.pos_grid(3)-1)*0.5; % axial position for the simulated results, relative to transducer position
        % plot(sim_res_axial_position, ...
        %     pred_axial_pressure_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
        % plot(exp_acoustic_profile(:,1),exp_acoustic_profile(:,2))
        % hold off
        % xline(opt_res.parameters.expected_focal_distance_mm, '--');
        % yline(desired_intensity, '--');
        % legend('Original', sprintf('Optimized for %2.f mm distance, analytical', opt_res.parameters.expected_focal_distance_mm), ...
        %     sprintf('Optimized for %2.f mm distance, simulated', opt_res.parameters.expected_focal_distance_mm),'Real profile','Location', 'best')
        % fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',sim_res_axial_position(pred_axial_pressure_opt==max(pred_axial_pressure_opt)))

        %% run on the head

        opt_parameters = load_parameters('bob_config_opt_CTX250-001_64.5mm.yaml'); 
        opt_parameters.transducer.source_amp = opt_source_amp;
        opt_parameters.transducer.source_phase_rad = [0 opt_phases];
        opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;

        opt_parameters.simulation_medium = 'layered';
        opt_parameters.run_posthoc_water_sims = 1;
        opt_parameters.run_heating_sims = 0;

        opt_parameters.interactive = 0;
        opt_parameters.overwrite_files = 'never';
        if best_trans_pos.idx == 12155115 || best_trans_pos.idx == 12592136
            opt_parameters.overwrite_files = 'always';
        end

        % subject_id = 3;

        opt_parameters.transducer.pos_t1_grid = round(table2array(best_trans_pos(1,["trans_x","trans_y","trans_z"])));
        opt_parameters.focus_pos_t1_grid = round(table2array(best_trans_pos(1,["targ_x","targ_y","targ_z"])));
        opt_parameters.expected_focal_distance_mm = best_trans_pos.dist_to_target*pixel_size;

        opt_parameters.results_filename_affix = sprintf('_bob_%s_%i',target, best_trans_pos.idx);
        opt_parameters.run_heating_sims = 1;
        % single_subject_pipeline(subject_id, opt_parameters)
        qsubfeval(@single_subject_pipeline_wrapper, subject_id, opt_parameters, 'timreq',  60*60*7,  'memreq',  50*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"');
        close all

    end

    end
end