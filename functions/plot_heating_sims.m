function plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos, brain_img, CEM43)

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
%   brain_img         - [Nx x Ny x Nz] matrix representing the brain image used for visualization.
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

    %% Extract temperature profiles along the focal axis
    focal_axis_temperature = squeeze(focal_planeT(trans_pos(1), :, :));
    focal_axis_trise = squeeze(focal_planeT(trans_pos(1), :, :) - focal_planeT(trans_pos(1), :, 1));

    %% Plot temperature profile over time
    h = figure;
    plot([time_status_seq([time_status_seq.recorded] == 1).time], ...
        focal_axis_temperature(1:10:end, :)); % Plot every 10th voxel
    colormap('lines');
    xlabel('Time [s]');
    ylabel('Temperature [°C]');
    y_range = ylim();
    hold on;
    
    % Highlight 'on' periods in the simulation timeline
    for i = 2:length(time_status_seq) - 1
        if ~strcmp(time_status_seq(i).status, 'on')
            continue;
        end
        x_points = [time_status_seq(i-1).time, time_status_seq(i-1).time, ...
                    time_status_seq(i).time, time_status_seq(i).time];
        y_points = [y_range(1), y_range(2), y_range(2), y_range(1)];
        fill(x_points, y_points, [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    hold off;
    saveas(h, output_plot, 'png');
    close(h);

    %% Plot temperature rise over time
    h = figure;
    plot([time_status_seq([time_status_seq.recorded] == 1).time], ...
        focal_axis_trise(1:10:end, :)); % Plot every 10th voxel
    colormap('lines');
    xlabel('Time [s]');
    ylabel('Temperature Rise [°C]');
    y_range = ylim();
    hold on;
    
    % Highlight 'on' periods in the simulation timeline
    for i = 2:length(time_status_seq) - 1
        if ~strcmp(time_status_seq(i).status, 'on')
            continue;
        end
        x_points = [time_status_seq(i-1).time, time_status_seq(i-1).time, ...
                    time_status_seq(i).time, time_status_seq(i).time];
        y_points = [y_range(1), y_range(2), y_range(2), y_range(1)];
        fill(x_points, y_points, [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    hold off;
    saveas(h, output_plot_rise, 'png');
    close(h);

    %% Plot CEM43 values over time
    g = figure;
    focal_CEM = squeeze(CEM43(trans_pos(1), :, :));
    
    plot([time_status_seq([time_status_seq.recorded] == 1).time], ...
         max(focal_CEM(:, :), [], 1)); % Maximum CEM43 values over axial positions at each time point
         
    colormap('lines');
    xlabel('Time [s]');
    ylabel(sprintf('CEM43'));
    
    y_range = ylim();
    hold on;
    
    % Highlight 'on' periods in the simulation timeline
    for i = 2:length(time_status_seq) - 1
        if ~strcmp(time_status_seq(i).status, 'on')
            continue;
        end
        x_points = [time_status_seq(i-1).time, time_status_seq(i-1).time, ...
                    time_status_seq(i).time, time_status_seq(i).time];
        y_points = [y_range(1), y_range(2), y_range(2), y_range(1)];
        fill(x_points, y_points, [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    hold off;
    
    saveas(g, output_plotCEM, 'png');
    close(g);

end