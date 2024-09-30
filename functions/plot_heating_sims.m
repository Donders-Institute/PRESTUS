function plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos, brain_img, CEM43)
    
    %% Plots the results of the heating simulations in a line graph

    arguments
        focal_planeT (:,:,:)
        time_status_seq
        parameters struct
        trans_pos (1,3) 
        brain_img (:,:,:) double
        CEM43 (:,:,:)
    end
    %new_transducer_pos = parameters.transducer.pos_grid - heating_window_dims(1,:) + 1;
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_by_time%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    output_plot_rise = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heatrise_by_time%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    output_plotCEM = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_CEM_by_time%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));

    % assuming that transducer and focus are aligned the z-axis
    % transducer_plane_final_temperature = squeeze(kwaveDiffusion.T(:, new_transducer_pos(2), :));
    if gpuDeviceCount==0
        focal_planeT = gather(focal_planeT);
    end
    focal_axis_temperature = squeeze(focal_planeT(trans_pos(1),:,:));
    focal_axis_trise = squeeze(focal_planeT(trans_pos(1),:,:)-focal_planeT(trans_pos(1),:,1));

    %% plot temperature

    h = figure;
    plot([time_status_seq(find([time_status_seq.recorded]==1)).time], ...
        focal_axis_temperature(1:10:size(focal_axis_temperature, 1), 1:size(time_status_seq,2))'); % plot every 10th voxel   
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    hold on;
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
    saveas(h, output_plot, 'png');
    close(h);

    %% plot temperature rise

    h = figure;
    plot([time_status_seq(find([time_status_seq.recorded]==1)).time], ...
        focal_axis_trise(1:10:size(focal_axis_trise, 1), 1:size(time_status_seq,2))'); % plot every 10th voxel   
    colormap('lines')
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    hold on;
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
    saveas(h, output_plot_rise, 'png');
    close(h);
    
    %% plot CEM43

    g = figure;
    focal_CEM = squeeze(CEM43(trans_pos(1),:,1:numel(time_status_seq)));

    plot([time_status_seq(find([time_status_seq.recorded]==1)).time], ...
        max(focal_CEM(:, 1:size(time_status_seq,2)),[],1)); 
    colormap('lines')
    xlabel('Time [s]');
    ylabel(sprintf('CEM43'));

    y_range = ylim();
    hold on;
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
    saveas(g, output_plotCEM, 'png');
    close(g);

    %% Create a video of heating effects over the course of the experiment (request by default)

    if ~isfield(parameters, 'heatingvideo')
        parameters.heatingvideo = 1; % default: create & save video
    end
    if parameters.heatingvideo == 1
        color_limits = [min(focal_planeT(:)), max(focal_planeT(:))];
        brain_slice = mat2gray(squeeze(brain_img(:,parameters.transducer.pos_grid(2),:)));
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
