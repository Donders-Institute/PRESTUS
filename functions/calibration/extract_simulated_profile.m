function [pred_axial_pressure] = extract_simulated_profile(initial_res, initial_params, profile_empirical, desired_intensity, equipment_name)
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
%   initial_params    - structure containing simulation and grid parameters
%   profile_empirical - structure with empirical profile info (for labeling)
%   desired_intensity - numeric target intensity (for output naming)
%   equipment_name    - string identifying the experimental setup
%
% OUTPUT:
%   pred_axial_pressure - simulated pressure values along the focal axis
%
% NOTE:
%   The code works for both 2D and 3D simulations, depending on 
%   initial_params.n_sim_dims (2 or 3).
%--------------------------------------------------------------------------

    %% Retrieve and prepare pressure data

    % The simulation may output data as a GPU array, so gather converts it
    % into a standard MATLAB CPU array for post-processing.
    p_max = gather(initial_res.sensor_data.p_max_all);

    %% Prepare and plot a 2D map of the focal-plane pressure distribution

    figure;

    % Convert grid indices to physical distances in millimeters.
    p_distance = (1:size(p_max, 1)) * initial_params.grid_step_mm;
    p_width = (1:size(p_max, initial_params.n_sim_dims)) * initial_params.grid_step_mm;

    % Depending on the dimensionality, select an appropriate plane:
    if initial_params.n_sim_dims == 2
        % 2D simulation -> directly available pressure field.
        p_axialprofile = squeeze(p_max(:, :))';
    elseif initial_params.n_sim_dims == 3
        % 3D simulation -> take the central slice at the transducer’s lateral position.
        p_axialprofile = squeeze(p_max(:, initial_params.transducer.trans_pos(2), :))';
    else
        error('Unsupported simulation dimensionality: expected 2 or 3.');
    end

    % Plot the pressure map.
    imagesc(p_distance, p_width, p_axialprofile);
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    cb = colorbar; %#ok<NASGU> % keep reference in case of future annotations
    title('Pressure for the focal plane');

    clear p_distance p_width p_axialprofile;

    fig_path = fullfile(initial_params.outputs_folder, ...
        strcat('Initial_Intensity_map_2D_at_F_', ...
            num2str(profile_empirical.focus_wrt_exit_plane), ...
            '_at_I_', num2str(desired_intensity), ...
            '_', equipment_name, '.png'));

    % Save and close the figure.
    saveas(gcf, fig_path);
    close(gcf);

    %% Extract the simulated pressure profile along the focal axis

    % This section isolates the 1D pressure distribution along the acoustic beam.
    if initial_params.n_sim_dims == 2
        % For 2D: extract the axial profile along the transducer's horizontal position.
        i_x = initial_params.transducer.trans_pos(1);
        pred_axial_pressure = squeeze(p_max(i_x, :));
        clear i_x;

    elseif initial_params.n_sim_dims == 3
        % For 3D: extract along both lateral center coordinates.
        i_x = initial_params.transducer.trans_pos(1);
        i_y = initial_params.transducer.trans_pos(2);
        pred_axial_pressure = squeeze(p_max(i_x, i_y, :));
        clear i_x i_y;
    end

end
