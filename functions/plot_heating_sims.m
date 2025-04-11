function plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos, medium_masks, CEM43)

% PLOT_HEATING_SIMS Visualizes heating simulation results over time.
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

    %% Define output file paths for plots
    output_plot = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_heating_by_time%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    output_plot_rise = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_heatrise_by_time%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    output_plotCEM = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_CEM_by_time%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));

    %% Convert GPU arrays to CPU if necessary
    if gpuDeviceCount == 0
        focal_planeT = gather(focal_planeT);
    end

    %% Extract temperature profiles along the focal axis (note: not max.)
    focal_axis_temperature = squeeze(focal_planeT(trans_pos(1), :, :));
    focal_axis_trise = squeeze(focal_planeT(trans_pos(1), :, :) - focal_planeT(trans_pos(1), :, 1));
    focal_CEM = squeeze(CEM43(trans_pos(1),:,:));

    %% Extract tissue indices along the focal axis

    medium_masks_focal = squeeze(medium_masks(trans_pos(1),trans_pos(2),:));
    labels = fieldnames(parameters.layer_labels);
    skull_i = find(strcmp(labels, 'skull_cortical'));
    trabecular_i = find(strcmp(labels, 'skull_trabecular'));
    all_skull_ids = [skull_i, trabecular_i];
    mask_skull = ismember(medium_masks_focal,all_skull_ids);
    brain_i = find(strcmp(labels, 'brain'));
    mask_brain = ismember(medium_masks_focal,brain_i);
    skin_i = find(strcmp(labels, 'skin'));
    mask_skin = ismember(medium_masks_focal,skin_i);

    %% Plot temperature profile over time

    recordedtime = [time_status_seq(find([time_status_seq.recorded]==1)).time];

    h = figure;
    hold on;
    p1 = plot(recordedtime, nanmean(focal_axis_temperature(mask_skull, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
    p2 = plot(recordedtime, nanmean(focal_axis_temperature(mask_skin, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
    p3 = plot(recordedtime, nanmean(focal_axis_temperature(mask_brain, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
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
    hold off
    legend([p1, p2, p3], {'skull'; 'skin'; 'brain'});
    legend('boxoff');
    title('Temperature in focal plane (avg. in tissue)');
    saveas(h, output_plot, 'png');
    close(h);

    %% Plot temperature rise over time

    h = figure;
    hold on;
    p1 = plot(recordedtime, nanmean(focal_axis_trise(mask_skull, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
    p2 = plot(recordedtime, nanmean(focal_axis_trise(mask_skin, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
    p3 = plot(recordedtime, nanmean(focal_axis_trise(mask_brain, 1:size(time_status_seq,2)),1), 'LineWidth', 2);
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
    hold off
    legend([p1, p2, p3], {'skull'; 'skin'; 'brain'}, 'location', 'NorthWest');
    legend('boxoff');
    title('Temperature rise in focal plane (avg. in tissue)');
    saveas(h, output_plot_rise, 'png');
    close(h);

    %% Plot CEM43 values over time

    g = figure;
    hold on;
    p1 = plot(recordedtime, max(focal_CEM(mask_skull, 1:size(time_status_seq,2)),[],1), 'LineWidth', 2);
    p2 = plot(recordedtime, max(focal_CEM(mask_skin, 1:size(time_status_seq,2)),[],1), 'LineWidth', 2);
    p3 = plot(recordedtime, max(focal_CEM(mask_brain, 1:size(time_status_seq,2)),[],1), 'LineWidth', 2);
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
    hold off
    legend([p1, p2, p3], {'skull'; 'skin'; 'brain'}, 'location', 'NorthWest');
    legend('boxoff');
    title('CEM43 in focal plane (max. in tissue)');
    saveas(g, output_plotCEM, 'png');
    close(g);

    %% Create a video of heating effects over the course of the experiment (if requested)

    if ~isfield(parameters, 'heatingvideo')
        parameters.heatingvideo = 1; % default: create & save video
    end
    if parameters.heatingvideo == 1
        color_limits = [min(focal_planeT(:)), max(focal_planeT(:))];
        brain_slice = mat2gray(squeeze(medium_masks(:,parameters.transducer.pos_grid(2),:)));
        output_video_name = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_animation%s.avi', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        v = VideoWriter(output_video_name,'Uncompressed AVI');
        v.FrameRate = 2; % frames per second
        % Create a video writer object for the output video file and open the object for writing.
        open(v);
        tss_recorded = time_status_seq(find([time_status_seq.recorded]==1));
        %Generate a set of frames, get the frame from the figure, and then write each frame to the file.
        gray_img_brain = repmat(mat2gray(brain_slice),[1 1 3]);
        for k = 1:size(focal_planeT,3) 
            cur_temp = squeeze(focal_planeT(:, :,  k));
            rgbImage = ind2rgb(round((cur_temp-color_limits(1))/diff(color_limits)*256), viridis(256));
            rgbImage(cur_temp <=37.01) = nan;
            imshowpair(gray_img_brain, rgbImage, 'blend')
            %set(h, 'AlphaData', 0.2)
            try
                title(sprintf('Timepoint: %i, time: %.2f s, max temperature: %.2f C, transducer: %s', ...
                    k, tss_recorded(k).time, max(cur_temp(:)), tss_recorded(k).status), 'FontSize', 10)
            catch
                warning("Some parameter is out of bounds")
                disp(['k:', num2str(k)]);
                disp(['dim focal plane (3):', num2str(size(focal_planeT,3))]);
                disp(['dim tss_recorded:', num2str(size(tss_recorded) )]);
            end
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
    end

end