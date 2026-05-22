function thermal_plot_sim(focal_planeT, time_status_seq, parameters, trans_pos, medium_masks, CEM43, timeseries, CEM43_iso)
% THERMAL_PLOT_SIM  Visualise heating simulation results over time and save figures
%
% Generates six or more PNG files and one AVI video (if requested):
%   focal-axis temperature, temperature rise, and CEM43 (kWave & ISO) for
%   the plane through the transducer; and overall tissue-maximum versions
%   of each metric using the layer timeseries. Tissue masks are derived
%   from medium_masks via tissuemask_binary. Grey-shaded regions mark ON
%   periods in the sonication protocol. The video encodes focal-plane
%   temperature frames coloured by viridis, blended over the tissue mask.
%   Note: focal-axis plots track the single axial line through trans_pos,
%   not the spatial maximum — see the timeseries plots for true maxima.
%
% Use as:
%   thermal_plot_sim(focal_planeT, time_status_seq, parameters, trans_pos, medium_masks, CEM43, timeseries)
%   thermal_plot_sim(focal_planeT, time_status_seq, parameters, trans_pos, medium_masks, CEM43, timeseries, CEM43_iso)
%
% Input:
%   focal_planeT    - temperature in focal plane over time [°C]
%   time_status_seq - struct array with fields time [s], status, recorded
%   parameters      - PRESTUS config; must contain io.dir_output, io.output_affix,
%                     subject_id, simulation.medium, transducer(1).trans_pos,
%                     transducer(1).focus_pos, thermal.temp_0, io.save_heatingvideo
%   trans_pos       - [1x3] transducer position in grid indices
%   medium_masks    - layer label map
%   CEM43           - k-Wave CEM43 in focal plane over time [min]
%   timeseries      - struct from THERMAL_UPDATE_TIMESERIES with fields
%                     T, Tdiff, CEM43, CEM43_iso (each a struct of layer arrays)
%   CEM43_iso       - ISO CEM43 in focal plane over time [min] (optional, default: [])
%
% See also: THERMAL_SIMULATION, THERMAL_ANALYSIS, THERMAL_UPDATE_TIMESERIES

arguments
    focal_planeT    {mustBeNumeric}
    time_status_seq (1,:) struct
    parameters      (1,1) struct
    trans_pos       (1,:) {mustBeNumeric}
    medium_masks    {mustBeNumericOrLogical}
    CEM43           {mustBeNumeric}
    timeseries      (1,1) struct
    CEM43_iso       {mustBeNumeric} = []
end

    %% Define output file paths for plots
    td = parameters.io.dir_img;
    sid = parameters.subject_id;
    med = parameters.simulation.medium;
    aff = parameters.io.output_affix;
    output_plot      = fullfile(td, sprintf('sub-%03d_%s_thermal%s.png',     sid, med, aff));
    output_plot_rise = fullfile(td, sprintf('sub-%03d_%s_thermalrise%s.png', sid, med, aff));
    output_plot_CEM  = fullfile(td, sprintf('sub-%03d_%s_CEM%s.png',         sid, med, aff));

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
    target_pos_focal = parameters.transducer(1).focus_pos(end);

    mask = tissuemask_binary(parameters, medium_masks_focal);

    %% [focal axis] temperature over time

    HEAT.recordedtime = [time_status_seq(find([time_status_seq.recorded]==1)).time];

    h = figure;
    hold on;
    leg_handles = []; leg_labels = {};
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skull)
        data2plot = max(focal_axis_temperature(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.temp.Skull, 'LineWidth', 2);
    if any(mask.skull), leg_handles(end+1) = p1; leg_labels{end+1} = 'skull'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skin)
        data2plot = max(focal_axis_temperature(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.temp.Skin, 'LineWidth', 2);
    if any(mask.skin), leg_handles(end+1) = p2; leg_labels{end+1} = 'skin'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.brain)
        data2plot = max(focal_axis_temperature(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.temp.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.temp.Brain, 'LineWidth', 2);
    if any(mask.brain), leg_handles(end+1) = p3; leg_labels{end+1} = 'brain'; end
    HEAT.temp.TARGET = focal_axis_temperature(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.temp.TARGET, 'LineWidth', 2);
    leg_handles(end+1) = p4; leg_labels{end+1} = 'target';
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
    legend(leg_handles, leg_labels);
    legend('boxoff');
    title('Temperature in focal plane (NOT overall tissue max.)');
    saveas(h, output_plot, 'png');
    close(h);

    %% [focal axis] temperature rise over time

    h = figure;
    hold on;
    leg_handles = []; leg_labels = {};
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skull)
        data2plot = max(focal_axis_trise(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.tempRise.Skull, 'LineWidth', 2);
    if any(mask.skull), leg_handles(end+1) = p1; leg_labels{end+1} = 'skull'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skin)
        data2plot = max(focal_axis_trise(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.tempRise.Skin, 'LineWidth', 2);
    if any(mask.skin), leg_handles(end+1) = p2; leg_labels{end+1} = 'skin'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.brain)
        data2plot = max(focal_axis_trise(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.tempRise.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.tempRise.Brain, 'LineWidth', 2);
    if any(mask.brain), leg_handles(end+1) = p3; leg_labels{end+1} = 'brain'; end
    HEAT.tempRise.TARGET = focal_axis_trise(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.tempRise.TARGET, 'LineWidth', 2);
    leg_handles(end+1) = p4; leg_labels{end+1} = 'target';
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
    legend(leg_handles, leg_labels, 'location', 'NorthWest');
    legend('boxoff');
    title('Temperature rise in focal plane (NOT overall tissue max.)');
    saveas(h, output_plot_rise, 'png');
    close(h);

    %% [focal axis] CEM43 values over time

    g = figure;
    hold on;
    leg_handles = []; leg_labels = {};
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skull)
        data2plot = max(focal_CEM(mask.skull, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Skull = data2plot;
    p1 = plot(HEAT.recordedtime, HEAT.CEM43.Skull, 'LineWidth', 2);
    if any(mask.skull), leg_handles(end+1) = p1; leg_labels{end+1} = 'skull'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.skin)
        data2plot = max(focal_CEM(mask.skin, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Skin = data2plot;
    p2 = plot(HEAT.recordedtime, HEAT.CEM43.Skin, 'LineWidth', 2);
    if any(mask.skin), leg_handles(end+1) = p2; leg_labels{end+1} = 'skin'; end
    data2plot = NaN(numel(HEAT.recordedtime), 1);
    if any(mask.brain)
        data2plot = max(focal_CEM(mask.brain, 1:size(time_status_seq,2)),[],1);
    end
    HEAT.CEM43.Brain = data2plot;
    p3 = plot(HEAT.recordedtime, HEAT.CEM43.Brain, 'LineWidth', 2);
    if any(mask.brain), leg_handles(end+1) = p3; leg_labels{end+1} = 'brain'; end
    HEAT.CEM43.TARGET = focal_CEM(target_pos_focal, 1:size(time_status_seq,2));
    p4 = plot(HEAT.recordedtime, HEAT.CEM43.TARGET, 'LineWidth', 2);
    leg_handles(end+1) = p4; leg_labels{end+1} = 'target';
    colormap('lines')
    xlabel('Time [s]');
    ylabel('CEM43');
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
    legend(leg_handles, leg_labels, 'location', 'NorthWest');
    legend('boxoff');
    title('CEM43 (kWave) in focal plane (NOT overall tissue max.)');
    saveas(g, output_plot_CEM, 'png');
    close(g);

    %% [focal axis] ISO CEM43 values over time

    if nargin >= 8 && ~isempty(CEM43_iso)
        output_plot_CEM_iso = fullfile(td, sprintf('sub-%03d_%s_CEM_iso%s.png', sid, med, aff));
        focal_CEM_iso = squeeze(CEM43_iso(trans_pos(1),:,:));
        g = figure; hold on;
        leg_handles = []; leg_labels = {};
        data2plot = NaN(numel(HEAT.recordedtime), 1);
        if any(mask.skull)
            data2plot = max(focal_CEM_iso(mask.skull, 1:size(time_status_seq,2)),[],1);
        end
        p1 = plot(HEAT.recordedtime, data2plot, 'LineWidth', 2);
        if any(mask.skull), leg_handles(end+1) = p1; leg_labels{end+1} = 'skull'; end
        data2plot = NaN(numel(HEAT.recordedtime), 1);
        if any(mask.skin)
            data2plot = max(focal_CEM_iso(mask.skin, 1:size(time_status_seq,2)),[],1);
        end
        p2 = plot(HEAT.recordedtime, data2plot, 'LineWidth', 2);
        if any(mask.skin), leg_handles(end+1) = p2; leg_labels{end+1} = 'skin'; end
        data2plot = NaN(numel(HEAT.recordedtime), 1);
        if any(mask.brain)
            data2plot = max(focal_CEM_iso(mask.brain, 1:size(time_status_seq,2)),[],1);
        end
        p3 = plot(HEAT.recordedtime, data2plot, 'LineWidth', 2);
        if any(mask.brain), leg_handles(end+1) = p3; leg_labels{end+1} = 'brain'; end
        p4 = plot(HEAT.recordedtime, focal_CEM_iso(target_pos_focal, 1:size(time_status_seq,2)), 'LineWidth', 2);
        leg_handles(end+1) = p4; leg_labels{end+1} = 'target';
        colormap('lines')
        xlabel('Time [s]'); ylabel('CEM43 (ISO)');
        y_range = ylim();
        for i = 2:(length(time_status_seq)-1)
            if ~strcmp(time_status_seq(i).status,'on'), continue; end
            x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
            y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
            a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
            a.FaceAlpha = 0.1;
        end
        xlim([0, HEAT.recordedtime(end)]); hold off
        legend(leg_handles, leg_labels, 'location', 'NorthWest');
        legend('boxoff');
        title('CEM43 ISO in focal plane (NOT overall tissue max.)');
        saveas(g, output_plot_CEM_iso, 'png');
        close(g);
    end

    %% Plot timeseries for maximum in medium

    output_plot      = fullfile(td, sprintf('sub-%03d_%s_thermal_max%s.png',     sid, med, aff));
    output_plot_rise = fullfile(td, sprintf('sub-%03d_%s_thermalrise_max%s.png', sid, med, aff));
    output_plot_CEM  = fullfile(td, sprintf('sub-%03d_%s_CEM_max%s.png',         sid, med, aff));

    % [layer-max] temperature over time

    HEAT.recordedtime = [time_status_seq(find([time_status_seq.recorded]==1)).time];
    available_layers = fieldnames(timeseries.T);
    n_ts = numel(timeseries.T.(available_layers{1}));
    if numel(HEAT.recordedtime) == n_ts + 1
        HEAT.recordedtime = HEAT.recordedtime(2:end);
    elseif numel(HEAT.recordedtime) ~= n_ts
        warning('thermal_plot_sim:timeseriesLengthMismatch', ...
            'recordedtime length (%d) does not match timeseries length (%d); truncating to min.', ...
            numel(HEAT.recordedtime), n_ts);
        n_min = min(numel(HEAT.recordedtime), n_ts);
        HEAT.recordedtime = HEAT.recordedtime(end-n_min+1:end);
        for il = 1:numel(available_layers)
            timeseries.T.(available_layers{il}) = timeseries.T.(available_layers{il})(end-n_min+1:end);
        end
    end
    available_layer_labels = strrep(available_layers, '_', ' ');

    h = figure; hold on;
    layer_handles_T = gobjects(0);
    for i_layer = 1:numel(available_layers)
        layer_handles_T(end+1) = plot(HEAT.recordedtime, timeseries.T.(available_layers{i_layer}), 'LineWidth', 2);
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
    if ~isempty(layer_handles_T)
        legend(layer_handles_T, available_layer_labels);
        legend('boxoff');
    end
    title('Temperature (overall tissue max.)');
    saveas(h, [output_plot], 'png');
    close(h);

    % [layer-max] temperature change over time

    h = figure; hold on;
    layer_handles_Tdiff = gobjects(0);
    for i_layer = 1:numel(available_layers)
        layer_handles_Tdiff(end+1) = plot(HEAT.recordedtime, timeseries.Tdiff.(available_layers{i_layer}), 'LineWidth', 2);
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
    if ~isempty(layer_handles_Tdiff)
        legend(layer_handles_Tdiff, available_layer_labels);
        legend('boxoff');
    end
    title('Temperature change (overall tissue max.)');
    saveas(h, [output_plot_rise], 'png');
    close(h);

    % [layer-max] thermal index over time

    h = figure; hold on;
    layer_handles_CEM = gobjects(0);
    for i_layer = 1:numel(available_layers)
        layer_handles_CEM(end+1) = plot(HEAT.recordedtime, timeseries.CEM43.(available_layers{i_layer}), 'LineWidth', 2);
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
    if ~isempty(layer_handles_CEM)
        legend(layer_handles_CEM, available_layer_labels);
        legend('boxoff');
    end
    title('Thermal Index CEM43 kWave (overall tissue max.)');
    saveas(h, [output_plot_CEM], 'png');
    close(h);

    % [layer-max] ISO CEM43 over time
    if isfield(timeseries, 'CEM43_iso') && ~isempty(fieldnames(timeseries.CEM43_iso))
        output_plot_CEM_iso_max = fullfile(td, sprintf('sub-%03d_%s_CEM_iso_max%s.png', sid, med, aff));
        h = figure; hold on;
        layer_handles_CEM_iso = gobjects(0);
        layer_labels_CEM_iso = {};
        for i_layer = 1:numel(available_layers)
            if isfield(timeseries.CEM43_iso, available_layers{i_layer})
                layer_handles_CEM_iso(end+1) = plot(HEAT.recordedtime, timeseries.CEM43_iso.(available_layers{i_layer}), 'LineWidth', 2);
                layer_labels_CEM_iso{end+1} = available_layer_labels{i_layer};
            end
        end
        colormap('lines')
        xlabel('Time [s]'); ylabel('Thermal Index [CEM43 ISO]');
        y_range = ylim();
        for i = 2:(length(time_status_seq)-1)
            if ~strcmp(time_status_seq(i).status,'on'), continue; end
            x_points = [time_status_seq(i-1).time time_status_seq(i-1).time time_status_seq(i).time time_status_seq(i).time];
            y_points = [y_range(1) y_range(2) y_range(2) y_range(1)];
            a = fill(x_points, y_points, [0.5 0.5 0.5],'EdgeColor','none', 'FaceAlpha', 0.1);
            a.FaceAlpha = 0.1;
        end
        xlim([1, HEAT.recordedtime(end)]); hold off
        if ~isempty(layer_handles_CEM_iso)
            legend(layer_handles_CEM_iso, layer_labels_CEM_iso); legend('boxoff');
        end
        title('Thermal Index CEM43 ISO (overall tissue max.)');
        saveas(h, output_plot_CEM_iso_max, 'png');
        close(h);
    end

    %% Save values for post-hoc group analysis

    output_HEAT = fullfile(parameters.io.dir_cache, sprintf('sub-%03d_%s_HEAT%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
    save(output_HEAT, 'HEAT', 'timeseries');

    %% Create a video of sagittal in-plane heating (if requested)

    if ~isfield(parameters.io, 'save_heatingvideo')
        parameters.io.save_heatingvideo = 1; % default: create & save video
    end
    if parameters.io.save_heatingvideo == 1
        color_limits = [min(focal_planeT(:)), max(focal_planeT(:))];
        if ndims(medium_masks) == 3
            brain_slice = mat2gray(squeeze(medium_masks(:,parameters.transducer.trans_pos(2),:)));
        elseif ndims(medium_masks) == 2
            brain_slice = mat2gray(squeeze(medium_masks));
        end
        output_video_name = fullfile(parameters.io.dir_output,...
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
                warn("Some parameter is out of bounds");
                disp(['k:', num2str(k)]);
                disp(['dim focal plane (3):', num2str(size(focal_planeT,3))]);
                disp(['dim tss_recorded:', num2str(size(tss_recorded) )]);
            end
        end
        close(v);
    end

end