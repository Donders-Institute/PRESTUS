function plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos, brain_img)
    
    %% Plots the results of the heating simulations in a line graph

    arguments
        focal_planeT (:,:,:)
        time_status_seq
        parameters struct
        trans_pos (1,3) 
        brain_img (:,:,:) double
    end
    %new_transducer_pos = parameters.transducer.pos_grid - heating_window_dims(1,:) + 1;
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_by_time%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    % assuming that transducer and focus are aligned the z-axis
    % transducer_plane_final_temperature = squeeze(kwaveDiffusion.T(:, new_transducer_pos(2), :));
    
    focal_axis_temperature = squeeze(focal_planeT(trans_pos(1),:,:));
    
    figure
    plot([time_status_seq(find([time_status_seq.recorded]==1)).time], focal_axis_temperature(1:10:size(focal_axis_temperature, 1), :)'); % plot every 10th voxel
    
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
    export_fig(output_plot, '-native')
    close

    %% Creates a video of heating effects over the course of the experiment

    color_limits = [min(focal_planeT(:)), max(focal_planeT(:))];
    brain_slice = mat2gray(squeeze(brain_img(:,parameters.transducer.pos_grid(2),:)));

    output_video_name = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating%s.avi', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
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
        title(sprintf('Timepoint: %i, time: %.2f s, max temperature: %.2f C, transducer: %s', k, tss_recorded(k).time, max(cur_temp(:)), tss_recorded(k).status), 'FontSize', 10)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    close(v);
    
end