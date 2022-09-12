function plot_heating_sims(kwaveDiffusion, time_status_seq, parameters, heating_window_dims, brain_img)
    arguments
        kwaveDiffusion
        time_status_seq
        parameters struct
        heating_window_dims (2,3) double
        brain_img (:,:,:) double
    end
    temperature = reshape(kwaveDiffusion.sensor_data, [diff(heating_window_dims)+1 size(kwaveDiffusion.sensor_data, 2)]);
    temperature = padarray(temperature, [0 0 0 1], parameters.thermal.temp_0, 'pre');
    new_transducer_pos = parameters.transducer.pos_grid - heating_window_dims(1,:) + 1;
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_by_time%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    % assuming that transducer and focus are aligned the z-axis
    transducer_plane_final_temperature = squeeze(kwaveDiffusion.T(:, new_transducer_pos(2), :));
    
    focal_axis_temperature = squeeze(temperature(new_transducer_pos(1), new_transducer_pos(2),:,:));
    plot([time_status_seq.time], focal_axis_temperature(1:10:size(focal_axis_temperature, 1), :)');
    
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


    %% create video of heating
    color_limits = [min(temperature(:)), max(temperature(:))];
    brain_slice = mat2gray(squeeze(brain_img(:,parameters.transducer.pos_grid(2),:)));
    
    % expand temperature data to brain slice dimensions
    temperature_slice = squeeze(temperature(:, new_transducer_pos(2), :,  :));

    heating_window_slice = heating_window_dims(:,[1,3]);
    for i = 1:2
        temperature_slice = padarray(temperature_slice, heating_window_slice(1,i)-1, parameters.thermal.temp_0, 'pre');
        temperature_slice = padarray(temperature_slice, size(brain_slice, i) - heating_window_slice(2,i), parameters.thermal.temp_0, 'post');
    end
    
    output_video_name = fullfile(parameters.data_path,sprintf('sim_outputs/sub-%03d_%s_heating%s.avi', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    v = VideoWriter(output_video_name,'Uncompressed AVI');
    v.FrameRate = 2; % frames per second
    
    
    % Create a video writer object for the output video file and open the object for writing.
    open(v);

    %Generate a set of frames, get the frame from the figure, and then write each frame to the file.
    gray_img_brain = repmat(mat2gray(brain_slice),[1 1 3]);
    for k = 1:size(temperature,4) 
        cur_temp = squeeze(temperature_slice(:, :,  k));
        rgbImage = ind2rgb(round((cur_temp-color_limits(1))/diff(color_limits)*256), viridis(256));
        rgbImage(cur_temp <=37.01) = nan;
        imshowpair(gray_img_brain, rgbImage, 'blend')
        %set(h, 'AlphaData', 0.2)
        title(sprintf('Timepoint: %i, time: %.2f s, max temperature: %.2f C, transducer: %s', k, time_status_seq(k).time, max(cur_temp(:)), time_status_seq(k).status), 'FontSize', 10)
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    close(v);
end