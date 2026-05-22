function [profile_sim] = extract_simulated_profile(initial_res, parameters)
% EXTRACT_SIMULATED_PROFILE  Extract axial intensity profile from a k-Wave simulation result
%
% Gathers the p_max_all field, generates a 2D focal-plane intensity map,
% and extracts the axial pressure profile through the transducer centre.
% Works for both 2D and 3D simulations.
%
% Use as:
%   [profile_sim] = extract_simulated_profile(initial_res, parameters)
%
% Input:
%   initial_res - simulation results struct with initial_res.sensor_data.p_max_all [Pa]
%   parameters  - PRESTUS config; uses grid.resolution_mm [mm], grid.dims,
%                 transducer(1).trans_pos, calibration.desired_focal_distance_ep [mm],
%                 calibration.desired_intensity [W/cm²], calibration.equipment_name,
%                 calibration.prefix
%
% Output:
%   profile_sim - struct with fields:
%                 axial_intensity [W/cm²], axial_distance_bowl [mm], velocity [m/s]
%
% See also: EXTRACT_REAL_INTENSITY_PROFILE, COMPUTE_ONEIL_SOLUTION

arguments
    initial_res  (1,1) struct
    parameters   (1,1) struct
end

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
        p_axialprofile = squeeze(p_max(:, parameters.transducer(1).trans_pos(2), :))';
    else
        error('Unsupported simulation dimensionality: expected 2 or 3.');
    end

    % Location of the transducer bowl, exit plane, and focus (in mm)
    i_bowl = parameters.transducer(1).trans_pos(end)*parameters.grid.resolution_mm;
    i_ep = round(parameters.transducer(1).trans_pos(end)*parameters.grid.resolution_mm+...
        parameters.transducer(1).focal_distance_offset);
    i_focus = parameters.transducer(1).focus_pos(end)*parameters.grid.resolution_mm;

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

    img_folder = fullfile(parameters.io.outputs_folder, 'img_calibration');
    if ~exist(img_folder, 'dir'); mkdir(img_folder); end
    fig_path = fullfile(img_folder, ...
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
        i_x = parameters.transducer(1).trans_pos(1);
        i_axis = parameters.transducer(1).trans_pos(end);
        pred_axial_pressure = squeeze(p_max(i_x, i_axis:end));
        clear i_x;

    elseif numel(parameters.grid.dims) == 3
        % For 3D: extract along both lateral center coordinates.
        i_x = parameters.transducer(1).trans_pos(1);
        i_y = parameters.transducer(1).trans_pos(2);
        i_axis = parameters.transducer(1).trans_pos(end);
        pred_axial_pressure = squeeze(p_max(i_x, i_y, i_axis:end));
        clear i_x i_y;
    end

    % Convert pressure to intensities [W/cm^2]
    pred_axial_intensity = pred_axial_pressure .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;

    % Compute particle velocity [m/s]
    velocity = parameters.transducer(1).annular.elem_amp(1) /...
               (parameters.medium_properties.water.density * parameters.medium_properties.water.sound_speed);

    %% Collect profile and distance

    % caluclated distances in mm in the grid
    axial_trans_pos_mm = parameters.transducer(1).trans_pos(end)* parameters.grid.resolution_mm;
    axial_end_pos_mm = parameters.grid.default_dims(end) * parameters.grid.resolution_mm;
    axial_position_sim_mm = axial_trans_pos_mm:parameters.grid.resolution_mm:axial_end_pos_mm;
    % for distance from bowl, remove PML/transducer position
    axial_distance_bowl = axial_position_sim_mm-axial_position_sim_mm(1);

    profile_sim.axial_intensity = pred_axial_intensity;
    profile_sim.axial_distance_bowl = axial_distance_bowl;
    profile_sim.velocity = velocity;

end
