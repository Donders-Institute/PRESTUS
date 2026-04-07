function thermal_plot_sim(focal_planeT, time_status_seq, parameters, trans_pos, medium_masks, CEM43, timeseries)

% THERMAL_PLOT_SIM Visualizes heating simulation results over time.
%
% This function generates plots and visualizations for heating simulations, including:
%   1. Temperature profile over time.
%   2. Temperature rise over time.
%   3. CEM43 (Cumulative Equivalent Minutes at 43°C) values over time.
% The function also creates a video showing heating effects on a brain slice during the experiment.
%
% Input:
%   focal_planeT      - [Nx x Ny x Nt] matrix representing temperature values on the focal plane over time.
%   time_status_seq   - Struct array containing time points and status of the simulation (e.g., 'on', 'off').
%   parameters        - Struct containing simulation parameters (e.g., output directory, transducer position).
%   trans_pos         - [1x3] array specifying the transducer position in grid coordinates.
%   medium_masks      - [Nx x Ny x Nz] matrix representing the brain labels
%   CEM43             - [Nx x Ny x Nt] matrix representing CEM43 values over time.
%
% Output:
%   None. The function saves plots and a video in the specified output directory.
%   Note: The timeseries capture metrics in the focal axis of the
%   transducer, not the maximum of observed values!

    %% Define output file paths for plots
    output_plot = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_thermal%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    output_plot_rise = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_thermalrise%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    output_plot_CEM = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_CEM%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    %% Convert GPU arrays to CPU if necessary
    if gpuDeviceCount == 0
        focal_planeT = gather(focal_planeT);
    end

    %% Plot data for the focal axis
    
    %% Extract temperature profiles along the focal axis (note: not max.!)
    % Change definition to the location of the max. temperature increase in
    % the sagittal transducer plane

    focal_axis_temperature = squeeze(focal_planeT(trans_pos(1), :, :));
    focal_axis_trise = squeeze(focal_planeT(trans_pos(1), :, :) - focal_planeT(trans_pos(1), :, 1));
    focal_CEM = squeeze(CEM43(trans_pos(1),:,:));

    %% Extract tissue indices along the focal axis

    if ndims(medium_masks) == 3
        medium_masks_focal = squeeze(medium_masks(trans_pos(1),trans_pos(2),:));
    else
        medium_masks_focal = squeeze(medium_masks(trans_pos(1),:));
    end
    % extract position of focus
    target_pos_focal = parameters.transducer(1).position.focus_pos(end);

    mask = tissuemask_binary(parameters, medium_masks_focal);

    %% [focal axis] temperature over time

    HEAT.recordedtime = [time_status_seq(find([time_status_seq.recorded]==1)).time];

    h = figure;
    hold on;
    % create dummy values to avoid crashes when mask is empty
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(find(mask.skull))
        data2plot = max(focal_axis_temperature(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.temp.Skull, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.skin))
        data2plot = max(focal_axis_temperature(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.temp.Skin, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.brain))
        data2plot = max(focal_axis_temperature(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.temp.Brain, 'LineWidth', 2);
    % add target
    HEAT.temp.TARGET = focal_axis_temperature(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.temp.TARGET, 'LineWidth', 2);
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([0, HEAT.recordedtime(end)]);
    hold off
    legend([p1, p2, p3, p4], {'skull'; 'skin'; 'brain'; 'target'});
    legend('boxoff');
    title('Temperature in focal plane (NOT overall tissue max.)');
    saveas(h, output_plot, 'png');
    close(h);

    %% [focal axis] temperature rise over time

    h = figure;
    hold on;
%     p1 = plot(HEAT.recordedtime, nanmean(focal_axis_trise(mask.skull, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
%     p2 = plot(HEAT.recordedtime, nanmean(focal_axis_trise(mask.skin, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
%     p3 = plot(HEAT.recordedtime, nanmean(focal_axis_trise(mask.brain, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
    % create dummy values to avoid crashes when mask is empty
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(find(mask.skull))
        data2plot = max(focal_axis_trise(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.tempRise.Skull, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.skin))
        data2plot = max(focal_axis_trise(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.tempRise.Skin, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.brain))
        data2plot = max(focal_axis_trise(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.tempRise.Brain, 'LineWidth', 2);
    % add target
    HEAT.tempRise.TARGET = focal_axis_trise(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.tempRise.TARGET, 'LineWidth', 2);
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([0, HEAT.recordedtime(end)]);
    hold off
    legend([p1, p2, p3, p4], {'skull'; 'skin'; 'brain'; 'target'}, 'location', 'NorthWest');
    legend('boxoff');
    title('Temperature rise in focal plane (NOT overall tissue max.)');
    saveas(h, output_plot_rise, 'png');
    close(h);

    %% [focal axis] CEM43 values over time

    g = figure;
    hold on;
    % create dummy values to avoid crashes when mask is empty
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(find(mask.skull))
        data2plot = max(focal_CEM(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.CEM43.Skull, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.skin))
        data2plot = max(focal_CEM(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.CEM43.Skin, 'LineWidth', 2);
    data2plot = NaN(numel(HEAT.recordedtime), 1); % create dummy values to avoid crashes
    if any(find(mask.brain))
        data2plot = max(focal_CEM(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.CEM43.Brain, 'LineWidth', 2);
    % add target
    HEAT.CEM43.TARGET = focal_CEM(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.CEM43.TARGET, 'LineWidth', 2);
    colormap('lines')
    xlabel('Time [s]');
    ylabel(sprintf('CEM43'));
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([0, HEAT.recordedtime(end)]);
    hold off
    legend([p1, p2, p3, p4], {'skull'; 'skin'; 'brain'; 'target'}, 'location', 'NorthWest');
    legend('boxoff');
    title('CEM43 in focal plane (NOT overall tissue max.)');
    saveas(g, output_plot_CEM, 'png');
    close(g);

    %% Plot timeseries for maximum in medium

    output_plot = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_thermal_max%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    output_plot_rise = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_thermalrise_max%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    output_plot_CEM = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_CEM_max%s.png', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));

    % [layer-max] temperature over time

    HEAT.recordedtime = [time_status_seq(find([time_status_seq.recorded]==1)).time];
    HEAT.recordedtime = HEAT.recordedtime(2:end);
    available_layers = fieldnames(timeseries.T);

    h = figure; hold on;
    for i_layer = 1:numel(available_layers)
        plot(HEAT.recordedtime, timeseries.T.(available_layers{i_layer}), 'LineWidth', 2);
    end
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([1, HEAT.recordedtime(end)]);
    hold off
    legend(available_layers);
    legend('boxoff');
    title('Temperature (overall tissue max.)');
    saveas(h, [output_plot], 'png');
    close(h);

    % [layer-max] temperature change over time

    h = figure; hold on;
    for i_layer = 1:numel(available_layers)
        plot(HEAT.recordedtime, timeseries.Tdiff.(available_layers{i_layer}), 'LineWidth', 2);
    end
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature Change [°C]');
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([1, HEAT.recordedtime(end)]);
    hold off
    legend(available_layers);
    legend('boxoff');
    title('Temperature change (overall tissue max.)');
    saveas(h, [output_plot_rise], 'png');
    close(h);

    % [layer-max] thermal index over time

    h = figure; hold on;
    for i_layer = 1:numel(available_layers)
        plot(HEAT.recordedtime, timeseries.CEM43.(available_layers{i_layer}), 'LineWidth', 2);
    end
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Thermal Index [CEM43]');
    y_range = ylim();
    for i = 2:(length(time_status_seq)-1)
        if ~strcmp(time_status_seq(i).status,'on')
            continue
        end
        x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
        y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
        a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
        a.FaceAlpha = 0.1;
    end
    xlim([1, HEAT.recordedtime(end)]);
    hold off
    legend(available_layers);
    legend('boxoff');
    title('Thermal Index CEM43 (overall tissue max.)');
    saveas(h, [output_plot_CEM], 'png');
    close(h);

    %% Save values for post-hoc group analysis

    output_HEAT = fullfile(parameters.io.output_dir, sprintf('sub-%03d_%s_HEAT%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    save(output_HEAT, 'HEAT', 'timeseries');

    %% Create a video of sagittal in-plane heating (if requested)

    if ~isfield(parameters.io, 'save_heatingvideo')
        parameters.io.save_heatingvideo = 1; % default: create & save video
    end
    if parameters.io.save_heatingvideo == 1
        color_limits = [min(focal_planeT(:)), max(focal_planeT(:))];
        if ndims(medium_masks) == 3
            brain_slice = mat2gray(squeeze(medium_masks(:,parameters.transducer.position.trans_pos(2),:)));
        elseif ndims(medium_masks) == 2
            brain_slice = mat2gray(squeeze(medium_masks));
        end
        output_video_name = fullfile(parameters.io.output_dir,...
            sprintf('sub-%03d_%s_heating_animation%s.avi', ...
            parameters.subject_id, ...
            parameters.simulation.medium, ...
            parameters.io.output_affix));
        v = VideoWriter(output_video_name,'Uncompressed AVI');
        v.FrameRate = 2; % frames per second
        % Create a video writer object for the output video file and open the object for writing.
        open(v);
        g = figure('Position', [100 100 304 156]);
        tss_recorded = time_status_seq(find([time_status_seq.recorded]==1));
        %Generate a set of frames, get the frame from the figure, and then write each frame to the file.
        gray_img_brain = repmat(mat2gray(brain_slice),[1 1 3]);
        for k = 1:size(focal_planeT,3) 
            cur_temp = squeeze(focal_planeT(:, :,  k));
            rgbImage = ind2rgb(round((cur_temp-color_limits(1))/diff(color_limits)*256), viridis(256));
            rgbImage(cur_temp <=37.01) = NaN;
            imshowpair(gray_img_brain, rgbImage, 'blend'),
            %set(h, 'AlphaData', 0.2)
            try
                title(sprintf('sample: %i, time: %.2f s\nmax temp: %.2f C, status: %s', ...
                    k, tss_recorded(k).time, max(cur_temp(:)), tss_recorded(k).status), 'FontSize', 8);
                frame = getframe(gcf);
                writeVideo(v,frame);
            catch
                warning("Some parameter is out of bounds");
                disp(['k:', num2str(k)]);
                disp(['dim focal plane (3):', num2str(size(focal_planeT,3))]);
                disp(['dim tss_recorded:', num2str(size(tss_recorded) )]);
            end
        end
        close(v);
    end

end