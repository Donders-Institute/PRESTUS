function [profile_sim] = extract_simulated_profile(initial_res, parameters)
%--------------------------------------------------------------------------
%
% PURPOSE:
%   This function extracts and visualizes simulated acoustic pressure data
%   (e.g., from a k-Wave simulation). It generates a 2D intensity map at
%   the focal plane and extracts the axial pressure profile along the
%   acoustic axis passing through the transducer center.
%
% INPUTS:
%   initial_res       - simulation results structure (contains sensor_data)
%   parameters    - structure containing simulation and grid parameters
%       .calibration.desired_focal_distance_ep - desired focal distance (for labeling)
%       .calibration.desired_intensity - numeric target intensity (for output naming)
%       .calibration.equipment_name    - string identifying the experimental setup
%       .calibration.prefix            - prefix for figure output
%
% OUTPUT:
%   profile_sim.axial_intensity - simulated intensity values along the focal axis (from transducer bowl) [W/cm^2]
%   profile_sim.axial_distance_bowl - corresponding distances from transducer bowl [mm]
%   profile_sim.velocity - simulated particle velocity [m/s]
%
% NOTE:
%   The code works for both 2D and 3D simulations
%--------------------------------------------------------------------------

    %% Retrieve and prepare pressure data

    % The simulation may output data as a GPU array, so gather converts it
    % into a standard MATLAB CPU array for post-processing.
    p_max = gather(initial_res.sensor_data.p_max_all);

    %% Prepare and plot a 2D map of the focal-plane pressure distribution

    figure;

    % Convert grid indices to physical distances in millimeters.
    p_distance = (1:size(p_max, 1)) * parameters.grid.resolution_mm;
    p_width = (1:size(p_max, numel(parameters.grid.dims))) * parameters.grid.resolution_mm;

    % Depending on the dimensionality, select an appropriate plane:
    if numel(parameters.grid.dims) == 2
        % 2D simulation -> directly available pressure field.
        p_axialprofile = squeeze(p_max(:, :))';
    elseif numel(parameters.grid.dims) == 3
        % 3D simulation -> take the central slice at the transducer’s lateral position.
        p_axialprofile = squeeze(p_max(:, parameters.transducer.position.trans_pos(2), :))';
    else
        error('Unsupported simulation dimensionality: expected 2 or 3.');
    end

    % Location of the transducer bowl, exit plane, and focus (in mm)
    i_bowl = parameters.transducer.position.trans_pos(end)*parameters.grid.resolution_mm;
    i_ep = round(parameters.transducer.position.trans_pos(end)*parameters.grid.resolution_mm+...
        parameters.transducer.focal_distance_offset);
    i_focus = parameters.transducer.position.focus_pos(end)*parameters.grid.resolution_mm;

    % Plot the pressure map
    imagesc(p_distance, p_width, p_axialprofile);
    axis image;
    hold on;
    yline(i_bowl, 'Color',[0 0 0])
    yline(i_ep, 'Color',[1 0 0])
    yline(i_focus, 'Color',[1 1 1])
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position in water medium (incl. PML) [mm]');
    colorbar;
    title('Pressure for the focal plane');

    clear p_distance p_width p_axialprofile;

    fig_path = fullfile(parameters.outputs_folder, ...
        strcat(parameters.calibration.prefix, 'Intensity_map_2D_at_F_', ...
            num2str(parameters.calibration.desired_focal_distance_ep), ...
            '_at_I_', num2str(parameters.calibration.desired_intensity), ...
            '_', parameters.calibration.equipment_name, '.png'));

    % Save and close the figure.
    saveas(gcf, fig_path);
    close(gcf);

    %% Extract the simulated pressure profile along the focal axis

    % This section isolates the 1D pressure distribution along the acoustic beam.
    if numel(parameters.grid.dims) == 2
        % For 2D: extract the axial profile along the transducer's horizontal position.
        i_x = parameters.transducer.position.trans_pos(1);
        i_axis = parameters.transducer.position.trans_pos(end);
        pred_axial_pressure = squeeze(p_max(i_x, i_axis:end));
        clear i_x;

    elseif numel(parameters.grid.dims) == 3
        % For 3D: extract along both lateral center coordinates.
        i_x = parameters.transducer.position.trans_pos(1);
        i_y = parameters.transducer.position.trans_pos(2);
        i_axis = parameters.transducer.position.trans_pos(end);
        pred_axial_pressure = squeeze(p_max(i_x, i_y, i_axis:end));
        clear i_x i_y;
    end

    % Convert pressure to intensities [W/cm^2]
    pred_axial_intensity = pred_axial_pressure .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;

    % Compute particle velocity [m/s]
    velocity = parameters.transducer.source_amp(1) / ...
               (parameters.medium_properties.water.density * parameters.medium_properties.water.sound_speed);

    %% Collect profile and distance

    % caluclated distances in mm in the grid
    axial_trans_pos_mm = parameters.transducer.position.trans_pos(end)* parameters.grid.resolution_mm;
    axial_end_pos_mm = parameters.grid.default_dims(end) * parameters.grid.resolution_mm;
    axial_position_sim_mm = axial_trans_pos_mm:parameters.grid.resolution_mm:axial_end_pos_mm;
    % for distance from bowl, remove PML/transducer position
    axial_distance_bowl = axial_position_sim_mm-axial_position_sim_mm(1);

    profile_sim.axial_intensity = pred_axial_intensity;
    profile_sim.axial_distance_bowl = axial_distance_bowl;
    profile_sim.velocity = velocity;

end
