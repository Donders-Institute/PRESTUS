function [opt_source_amp, opt_phases] = transducer_calibration(pn, parameters, transducer_name, subject_id, desired_intensity)

    % TO DO: specify necessary input parameters

    %% Simulations in free water

    % Now, we start with the actual simulations. The first thing with the simulations 
    % is to ensure that you have the right settings. For the real US transducers, 
    % the settings are calibrated by the manufacturer. The company we buy the transducers 
    % from calibrates them to have a given maximum intensity and a given location 
    % of the acoustic peak (specifically, the center of its half-maximum range). 
    % 
    % Below, we read in two columns from the calibration sheet provided by the manufacturer. 
    % Specifically, we take the distance from the transducer axis and the intensity 
    % measured at each of these points. When using transducers with multiple elements, 
    % the depth of the acoustic peak (and thus the focus) can be steered by changing 
    % the set of phases for the individual transducer elements with respect to one-another. 
    
    real_profile = readmatrix(fullfile(pn.configs, [transducer_name, '.csv']));
    
    % if the peak is calibrated to less than the desired intensity, scale
    % profile linearly
    
    scaling = desired_intensity/max(real_profile(:,2));
    real_profile(:,2) = real_profile(:,2).*scaling;
    
    dist_to_exit_plane = round(parameters.transducer.curv_radius_mm-...
        parameters.transducer.dist_to_plane_mm); % from the transducer specifications

    % need to add distance from transducer to exit plane to profiles
    parameters.expected_focal_distance_mm = parameters.expected_focal_distance_mm + dist_to_exit_plane;
    real_profile(:,1) = dist_to_exit_plane + real_profile(:,1);

    h = figure('Position', [10 10 900 500]);
    % calibrations done at the surface, we have to add distance to transducer
    plot(real_profile(:,1),real_profile(:,2), 'linewidth', 2)
    
    halfMax = (min(real_profile(:,2)) + max(real_profile(:,2))) / 2;
    % Find where the data first drops below half the max.
    % Hack to exclude the first few points
    index1 = 10+find(real_profile(11:end,2) >= halfMax, 1, 'first');
    % Find where the data last rises above half the max.
    index2 = find(real_profile(:,2) >= halfMax, 1, 'last');
    
    flhm = real_profile(index2,1) - real_profile(index1,1);
    flhm_center_x = (real_profile(index2,1) - real_profile(index1,1))/2+real_profile(index1,1);
    flhm_center_intensity = real_profile(find(real_profile(:,1)>=flhm_center_x,1,'first'), 2);
    xlabel('Axial Position [mm from transducer center surface]');
    ylabel('Intensity [W/cm^2]');
    xline(dist_to_exit_plane,'k--');
    text(dist_to_exit_plane+0.5, flhm_center_intensity+3, ...
        sprintf('Exit plane at %i mm',dist_to_exit_plane), "Color",'k');
    yline(desired_intensity, '--');
    yline(halfMax, '--');
    xline(real_profile(index1,1), '--');
    xline(real_profile(index2,1), '--');
    xline(flhm_center_x,'r--');
    text(flhm_center_x+0.5, flhm_center_intensity+3, ...
        sprintf('FLHM center intensity %.2f [W/cm^2] at %i mm',flhm_center_intensity,flhm_center_x), "Color",'r');
    expected_focus = parameters.expected_focal_distance_mm;
    intensity_at_expected_focus = mean(real_profile(real_profile(:,1)>=expected_focus-1&...
        real_profile(:,1)<=expected_focus+1,2));
    xline(expected_focus,'b--');
    text(expected_focus+0.5, intensity_at_expected_focus+4, ...
        sprintf('Expected focus intensity %.2f [W/cm^2] at %i mm',intensity_at_expected_focus,expected_focus),"Color",'b');
    [max_intensity, max_x] = max(real_profile(:,2));
    xline(real_profile(max_x,1),'--','Color','#7E2F8E');
    text(real_profile(max_x,1)+0.5, max_intensity+4, ...
        sprintf('Max intensity %.2f [W/cm^2] at %i mm',max_intensity,real_profile(max_x,1)),"Color",'#7E2F8E');
    plotname = fullfile(pn.outputs_folder, 'transducer_profile.png');
    saveas(h, plotname, 'png');
    close(h);
    
    %% 
    % These settings are transformed into a set of phases for the individual transducer 
    % elements. Unfortunately, each transducer is unique, so the manufacturer's settings 
    % are not necessarily the ones you want to use in your simulations. The approach 
    % used here is to have an acoustic profile (i.e., the intensity along the beam 
    % axis of the transducer) so that the maximum point is at the expected distance 
    % and has the expected pressure. To have that, we start with the simulations in 
    % the water medium to find the phase-angles for the transducer elements needed 
    % for the simulation. 
    
    parameters.simulation_medium = 'water'; % indicate that we only want the simulation in the water medium for now
    parameters.interactive = 0;
    parameters.overwrite_files = 'never';
    parameters.run_source_setup = 1;
    parameters.run_acoustic_sims = 1;
    parameters.run_heating_sims = 0;
    parameters.savemat = 1;
    if ismac % if run locally on a mac
        parameters.using_donders_hpc = 0;
        parameters.code_type = 'matlab_cpu';
    else
        parameters.using_donders_hpc = 1;
    end
    
    single_subject_pipeline(subject_id, parameters);
    
    %% 
    % After the simulations have finished, you can load the results. The results 
    % are saved under |sim_outputs| subdirectory of the data path that is defined 
    % in the config file along with some plots and summary statistics.
    % 
    % The data on pressure is saved in the |sensor_data| structure. First, we will 
    % plot a 2d pressure map by slicing the 3d pressure matrix at the focal plane. 
    
    % load results    
    load(fullfile(pn.outputs_folder, sprintf('sub-%03d_water_results%s.mat',...
        subject_id, parameters.results_filename_affix)),'sensor_data','parameters');
    
    % get maximum pressure
    p_max = gather(sensor_data.p_max_all); % transform from GPU array to normal array
    
    % plot 2d intensity map
    h = figure;
    imagesc((1:size(p_max,1))*parameters.grid_step_mm, ...
        (1:size(p_max,3))*parameters.grid_step_mm , ...
        squeeze(p_max(:,parameters.transducer.pos_grid(2),:))')
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    axis image;
    cb = colorbar;
    title('Pressure for the focal plane')
    plotname = fullfile(pn.outputs_folder, 'intensity_2d.png');
    saveas(h, plotname, 'png');
    close(h);
    
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
    axial_position = (1:parameters.default_grid_dims(3)).*parameters.grid_step_mm; % [mm]
    
    % evaluate pressure analytically
    % focusedAnnulusONeil provides an analytic solution for the pressure at the
    % focal (beam) axis
    n_elem = numel(parameters.transducer.Elements_ID_mm);
    [p_axial_oneil] = focusedAnnulusONeil(parameters.transducer.curv_radius_mm/1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, repmat(velocity,1,n_elem), ...
        parameters.transducer.source_phase_rad, parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, (axial_position-0.5)*1e-3);
    
    % plot focal axis pressure
    h = figure('Position', [10 10 900 500]);
    plot(axial_position, p_axial_oneil .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    xlabel('Axial Position [mm]');
    ylabel('Intensity [W/cm^2]');
    hold on
    plot(axial_position-(parameters.transducer.pos_grid(3)-1)*parameters.grid_step_mm, ...
        pred_axial_pressure.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4,'--');
    plot(real_profile(:,1),real_profile(:,2))
    hold off
    xline(parameters.expected_focal_distance_mm, '--');
    legend('Analytic solution','Simulated results','Real profile')
    title('Pressure along the beam axis')
    % what is distance to the maximum pressure?
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',...
        axial_position(p_axial_oneil==max(p_axial_oneil)))
    plotname = fullfile(pn.outputs_folder, 'focal_axis_pressure.png');
    saveas(h, plotname, 'png');
    close(h);

    % We want the simulated results to match the real profile as closely as possible, 
    % but as you can see that the two lines do not match exactly. There are two ways 
    % to deal with it, either you can use a newer kwaveArray class from k-wave (by 
    % setting |use_kWaveArray| flag in the config to 1) or you can introduce a correction 
    % factor as we do here. The kwaveArray solution is in theory better, but in practice 
    % it requires more time and memory. So for this tutorial, we will use a correction 
    % factor.
    
    % compute the approximate adjustment from simulated (on a grid) to analytic solution
    simulated_grid_adj_factor = max(pred_axial_pressure(:))/max(p_axial_oneil(:));
    
    %% Optimize for a given distance and pressure
    
    % So how to find the settings for the simulations that match the desired pressure 
    % and distance? It is easy to do, given that there is an analytic solution. For 
    % our subject, we assume that we know where the transducer is positioned and where 
    % we want to have the maximum pressure, so we know the distance at which the pressure 
    % should be maximal. We need to find the set of phases for transducer elements 
    % that would give the maximum pressure at the distance given in the manufacturers 
    % calibration sheet. We do so by searching through the parameter space (that is, 
    % varying the phases) as to minimize the error in distance to maximum pressure 
    % point. 
    
    opt_limits = [real_profile(2,1), real_profile(end,1)];

    [flhm_center, flhm_center_index] = get_flhm_center_position(real_profile(:,1), real_profile(:,2));

    % narrowness = .25; % the lower, the broader the gaussian weight filter
    % weights = normpdf(real_profile(:,1), real_profile(flhm_center_index,1), real_profile(flhm_center_index,1)/narrowness);
    weights = 1;

    % alternative: use Margely's implementation (MATLAB native functions)
    %[opt_phases, opt_velocity, min_err] = perform_global_search(parameters, dist_exit_plane, adjusted_profile_focus, velocity)

    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(...,
        phases_and_velocity(1:parameters.transducer.n_elements-1), ...
        parameters, ...
        phases_and_velocity(parameters.transducer.n_elements), ...
        real_profile(:,1), ...
        real_profile(:,2),...
        0, ...
        opt_limits, ...
        weights);
    
    rng(100,'twister') % setting seed for consistency
    
    func = optimize_phases;
    x0 = [randi(360, [1 n_elem-1])/180*pi velocity];
    lb = zeros(1,n_elem);
    ub = [2*pi*ones(1,n_elem-1) 0.2];
    options = setoptimoptions('popsize',5000, 'FinDiffType', 'central'); % number of initializations [, 'TolCon', 1e-8]
    [opt_phases_and_velocity, min_err] = minimize(func, x0, [],[],[],[],lb, ub, [], options);
    
    % plot optimization results
    [~, ~, ~, h] = phase_optimization_annulus_full_curve(...
        opt_phases_and_velocity(1:end-1), ...
        parameters, ...
        opt_phases_and_velocity(end),...
        real_profile(:,1), ...
        real_profile(:,2), ...
        1, ...
        opt_limits, ...
        weights);
    
    plotname = fullfile(pn.outputs_folder, 'optimized_profile.png');
    saveas(h, plotname, 'png');
    close(h);

    fprintf('Optimal phases: %s deg.; velocity: %.2f; optimization error: %.2f', ...
        mat2str(round(opt_phases_and_velocity(1:end-1)/pi*180)), opt_phases_and_velocity(end), min_err)

    % The left plot above shows the real and the fitted profiles along with the 
    % cost function used for fitting, while the right shows the error in fitting (the 
    % squared difference between the real and the fitted profile weighted by the cost 
    % function). 
    % 
    % The phase-angles given above should be copied into the config under the parameter 
    % "source_phase_deg". Keep in mind that different settings should be calculated 
    % and used for every transducer and every focus distance.
    % 
    % We also need to know the pressure for the simulated transducer so that in 
    % the water, the intensity (ISPPA) would be 30 W/cm^2. To do so, we first compute 
    % the analytic solution with the new phases but the original pressure, and then 
    % use it to adjust the simulated pressure and recompute the analytic solution 
    % (the latter is just for plotting). 
    
    opt_phases = opt_phases_and_velocity(1:end-1);
    opt_velocity = opt_phases_and_velocity(end);
    
    [p_axial_oneil_opt] = focusedAnnulusONeil(...
        parameters.transducer.curv_radius_mm/1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm]/1e3, ...
        repmat(opt_velocity,1,n_elem), ...
        [0 opt_phases], ...
        parameters.transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
        (axial_position-parameters.grid_step_mm)*1e-3);
    
    h = figure('Position', [10 10 900 500]);
    plot(axial_position, p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    xlabel('Axial Position [mm]');
    ylabel('Intensity [W/cm^2]');
    hold on
    plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    plot(real_profile(:,1),real_profile(:,2))
    hold off
    xline(parameters.expected_focal_distance_mm, '--');
    yline(30, '--');
    legend('Original simulation', sprintf('Optimized to match the real profile'),'Real profile')
    title('Pressure along the beam axis')
    plotname = fullfile(pn.outputs_folder, 'real_fitted_profiles.png');
    saveas(h, plotname, 'png');
    close(h);
    
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', ...
        axial_position(p_axial_oneil_opt==max(p_axial_oneil_opt)))
    try
        fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', ...
            get_flhm_center_position(axial_position, p_axial_oneil_opt))
    catch
        disp('FWHM cannot be estimated...');
    end
    
    opt_source_amp = round(opt_velocity/velocity*parameters.transducer.source_amp/simulated_grid_adj_factor);
    
    %% Simulate again in water to check the optimization results

    % Now we will update the settings with the optimized parameters, and rerun the simulations. 
    
    opt_parameters = load_parameters(['config_',transducer_name,'.yaml']); % load the configuration file

    opt_parameters.run_source_setup = 1;
    opt_parameters.run_acoustic_sims = 1;
    opt_parameters.run_heating_sims = 0;
    opt_parameters.data_path = parameters.data_path;
    opt_parameters.seg_path = parameters.seg_path;
    opt_parameters.sim_path = pn.data_sims;
    opt_parameters.transducer.source_amp = opt_source_amp;
    opt_parameters.transducer.source_phase_rad = [0 opt_phases];
    opt_parameters.transducer.source_phase_deg = [0 opt_phases]/pi*180;
    opt_parameters.expected_focal_distance_mm = parameters.expected_focal_distance_mm;
    opt_parameters.results_filename_affix = '_optimized';
    opt_parameters.simulation_medium = 'water';
    opt_parameters.overwrite_files = 'always';
    opt_parameters.interactive = 0;
    opt_parameters.usepseudoCT = 0;
    opt_parameters.savemat = 1;

    single_subject_pipeline(subject_id, opt_parameters);

    % After the simulations are complete, we plot them alongside the analytic solution. 
    
    opt_res = load(fullfile(pn.outputs_folder, ...
        sprintf('sub-%03d_water_results%s.mat', subject_id, opt_parameters.results_filename_affix)),...
        'sensor_data','parameters');
    
    % get maximum pressure
    p_max = gather(opt_res.sensor_data.p_max_all);
    pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1), opt_res.parameters.transducer.pos_grid(2),:));
    
    h = figure('Position', [10 10 900 500]);
    hold on
    plot(axial_position, p_axial_oneil.^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    xlabel('Axial Position [mm]');
    ylabel('Intensity [W/cm^2]');
    plot(axial_position, p_axial_oneil_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    
    sim_res_axial_position = axial_position-(opt_res.parameters.transducer.pos_grid(3)-1)*parameters.grid_step_mm; % axial position for the simulated results, relative to transducer position
    plot(sim_res_axial_position, ...
        pred_axial_pressure_opt .^2/(2*parameters.medium.water.sound_speed*parameters.medium.water.density) .* 1e-4);
    plot(real_profile(:,1),real_profile(:,2))
    hold off
    xline(opt_res.parameters.expected_focal_distance_mm, '--');
    yline(desired_intensity, '--');
    legend('Original simulation', sprintf('Optimized for %2.f mm distance, analytical', opt_res.parameters.expected_focal_distance_mm), ...
        sprintf('Optimized for %2.f mm distance, simulated', opt_res.parameters.expected_focal_distance_mm),'Real profile','Location', 'best')
    plotname = fullfile(pn.outputs_folder, 'simulation_analytic.png');
    saveas(h, plotname, 'png');
    close(h);
    
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n',...
        sim_res_axial_position(pred_axial_pressure_opt==max(pred_axial_pressure_opt)))

    % Now the simulated acoustic profile should have the desired distance 
    % to maximum and intensity, which means we can start with the simulations for 
    % the actual brain.
end