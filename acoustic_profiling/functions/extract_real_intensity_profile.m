function [norm_profile_focus, max_intens] = extract_real_intensity_profile(available_foci_wrt_exit_plane, focus_wrt_exit_plane, intens_data, prestus_dir, equipment_name, dist_from_tran, skip_front_peak_mm, local_output_dir, save_in_general_folder)        
    % Extracts or interpolates the intensity profile at a specific focal depth.
    %
    % Arguments:
    % - available_foci_wrt_exit_plane: Array of available focal depths relative to the exit plane [mm].
    % - focus_wrt_exit_plane: Desired focal depth relative to the exit plane [mm].
    % - intens_data: Matrix containing intensity profiles for different focal depths.
    % - prestus_dir: Directory path for saving results.
    % - equipment_name: Name of the equipment for labeling plots.
    % - dist_from_tran: Distance vector from the transducer [mm].
    % - skip_front_peak_mm: Distance to skip near-field peaks when finding the maximum intensity [mm].
    % - local_output_dir: local output directory if data is not saved in general PRESTUS output folder.
    % - save_in_general_folder: Option to save data in the general PRESTUS output folder.
    %
    % Returns:
    % - profile_focus: Extracted or interpolated intensity profile at the desired focal depth.
    % - max_intens: Maximum intensity in the profile beyond the specified skip distance.

    % Check if the exact focal depth is available
    col_index = find(available_foci_wrt_exit_plane == focus_wrt_exit_plane);

    if isempty(col_index)
        % Perform linear interpolation if the exact focus is not available
        [~, closestIndex] = min(abs(available_foci_wrt_exit_plane - focus_wrt_exit_plane));
        closest_foci_wrt_exit_plane = available_foci_wrt_exit_plane(closestIndex);

        % Determine neighboring focal depths for interpolation
        if closest_foci_wrt_exit_plane > focus_wrt_exit_plane
            closestIndex2 = closestIndex; % Higher focus
            closestIndex1 = closestIndex2 - 1; % Lower focus
        else
            closestIndex1 = closestIndex; % Lower focus
            closestIndex2 = closestIndex1 + 1; % Higher focus
        end

        % Validate indices for interpolation boundaries
        if closestIndex1 < 1 || closestIndex2 > length(available_foci_wrt_exit_plane)
            error('Focus is outside the range of available axial profiles. Interpolation not possible.')
        end

        % Retrieve intensity profiles for the neighboring focal depths
        focus_wrt_exit_plane_1 = round(available_foci_wrt_exit_plane(closestIndex1), 2);
        focus_wrt_exit_plane_2 = round(available_foci_wrt_exit_plane(closestIndex2), 2);
        profile_1 = intens_data(:, closestIndex1)';
        profile_2 = intens_data(:, closestIndex2)';

        % Normalize profiles by aligning their peaks to 0 (max of profile 1 and profile 2)
        [~, idx1] = max(profile_1);
        [~, idx2] = max(profile_2);

        x1_norm = dist_from_tran - dist_from_tran(idx1); % Align peak of profile 1 to 0
        x2_norm = dist_from_tran - dist_from_tran(idx2); % Align peak of profile 2 to 0

        % Define normalized common x-array
        x_common_norm = linspace(min(min(x1_norm), min(x2_norm)), ...
                                 max(max(x1_norm), max(x2_norm)), length(dist_from_tran));
     
        % Interpolate profiles in normalized space
        y1_interp_norm = interp1(x1_norm, profile_1, x_common_norm, 'spline', 0);
        y2_interp_norm = interp1(x2_norm, profile_2, x_common_norm, 'spline', 0);
        
        % Calculate weight (alpha) for interpolation based on the relative focal depths
        alpha = (focus_wrt_exit_plane - focus_wrt_exit_plane_1) / ...
                (focus_wrt_exit_plane_2 - focus_wrt_exit_plane_1);

        % Interpolate the profiles with the weight alpha
        norm_profile_focus  = (1-alpha) * y1_interp_norm + alpha * y2_interp_norm;

        % Map back to the original dist_from_tran
        x1_2_norm = x_common_norm + (alpha * x2_norm(idx2) + (1-alpha) * x1_norm(idx1));

        % Interpolate the final focused profile in the normalized space back to the original space
        mapped_profile_focus = interp1(x_common_norm, norm_profile_focus, x1_2_norm, 'spline', 0);

        % Calculate the offset (max_loc) to align the profile focus with dist_from_tran
        max_loc = abs(mean(x1_2_norm - dist_from_tran'));

        profile_focus = interp1(x1_2_norm + max_loc, mapped_profile_focus, dist_from_tran, 'spline', 0);

        % Plot the profiles and the interpolated result
        figure;
        plot(dist_from_tran, profile_1, '-o', 'DisplayName', ...
            ['Measurement 1, focus at ' num2str(focus_wrt_exit_plane_1)]);
        hold on;
        plot(dist_from_tran, profile_2, '-o', 'DisplayName', ...
            ['Measurement 2, focus at ' num2str(focus_wrt_exit_plane_2)]);
        plot(dist_from_tran, profile_focus, '-o', 'DisplayName', ...
            ['Interpolated, focus at ' num2str(focus_wrt_exit_plane)]);
        legend;
        xlabel('Distance wrt exit plane [mm]');
        ylabel('Intensity [W/cm^2]');
        title('Interpolation Input and Results');
    
    else
        % Use the exact profile if the focus matches an available value
        norm_profile_focus = intens_data(:, col_index)';

        % Plot the exact profile
        figure;
        plot(dist_from_tran, norm_profile_focus, '-o');
        xlabel('Distance wrt exit plane [mm]');
        ylabel('Intensity [W/cm^2]');
        title(['Axial Profile at Focus wrt Exit Plane: ' num2str(focus_wrt_exit_plane) ' [mm]']);
    end

    % Save the plot to the specified directory
    if save_in_general_folder
        fig_path = fullfile(prestus_dir, strcat('Interpolation_at_F_', num2str(focus_wrt_exit_plane), '_', equipment_name, '.png'));
    else
        fig_path = fullfile(local_output_dir, strcat('Interpolation_at_F_', num2str(focus_wrt_exit_plane), '_', equipment_name, '.png'));
    end

    saveas(gcf, fig_path);

    % Determine the maximum intensity beyond the specified skip distance 
    % for scaling to prevent catching max peak in near field peak.
    [~, closestIndex] = min(abs(dist_from_tran - skip_front_peak_mm));
    max_intens = max(norm_profile_focus(closestIndex:end));
    
end