function [norm_profile_focus, max_intens] = extract_real_intensity_profile(...
    parameters,...
    available_foci_wrt_exit_plane, ...
    desired_focal_distance_ep, ...
    intens_data, ...
    equipment_name, ...
    dist_from_exit_plane)

% EXTRACT_REAL_INTENSITY_PROFILE  Extract or interpolate an empirical axial intensity profile
%
% Selects the intensity profile for the desired focal depth from the
% available measured profiles. At an unmeasured depth the two neighbouring
% profiles are aligned on their focus and blended (INTERPOLATE_AXIAL_PROFILES);
% outside the measured range the nearest profile is shifted, up to
% calibration.max_extrapolation_mm (default 0: error). The profile is
% plotted and saved as a PNG, and a derived (interpolated/extrapolated)
% profile is also written as a profile CSV + sidecar (WRITE_AXIAL_PROFILES),
% both to calibration.path_output_profiles.
%
% Use as:
%   [norm_profile_focus, max_intens] = extract_real_intensity_profile( ...
%       parameters, available_foci_wrt_exit_plane, desired_focal_distance_ep, ...
%       intens_data, equipment_name, dist_from_exit_plane)
%
% Input:
%   parameters                  - PRESTUS config; uses calibration.skip_front_peak_mm [mm],
%                                 calibration.path_output_profiles and optionally
%                                 calibration.profile_align ('label' | 'peak') and
%                                 calibration.max_extrapolation_mm [mm]
%   available_foci_wrt_exit_plane - available focal depths relative to exit plane [mm]
%   desired_focal_distance_ep   - desired focal depth relative to exit plane [mm]
%   intens_data                 - matrix of intensity profiles per focal depth [W/cm²]
%   equipment_name              - equipment name for plot labelling
%   dist_from_exit_plane        - distance vector from transducer [mm]
%
% Output:
%   norm_profile_focus - [1 x N_z] extracted or interpolated intensity profile [W/cm²]
%                        on dist_from_exit_plane
%   max_intens         - maximum intensity beyond the skip-front-peak distance [W/cm²]
%
% See also: INTERPOLATE_AXIAL_PROFILES, WRITE_AXIAL_PROFILES, SCALE_REAL_INTENSITY_PROFILE,
%           CALIBRATION_TRANSDUCER

    cal = parameters.calibration;
    io  = struct('Align', get_or(cal, 'profile_align', 'label'), ...
                 'MaxExtrapolationMm', get_or(cal, 'max_extrapolation_mm', 0), ...
                 'NearFieldSkipMm', get_or(cal, 'skip_front_peak_mm', 10));
    measured = struct('depths_ep_mm', available_foci_wrt_exit_plane(:)', ...
                      'z_ep_mm', dist_from_exit_plane(:), 'I', intens_data);
    prof = interpolate_axial_profiles(measured, desired_focal_distance_ep, io);
    type = prof.profile_type{1};
    if strcmp(type, 'unavailable')
        error(['Focus %.2f mm is outside the range of available axial profiles [%.2f, %.2f] mm ' ...
            '(calibration.max_extrapolation_mm = %g). Interpolation not possible.'], ...
            desired_focal_distance_ep, min(available_foci_wrt_exit_plane), ...
            max(available_foci_wrt_exit_plane), io.MaxExtrapolationMm);
    end
    norm_profile_focus = prof.I(:, 1)';

    if ~exist(cal.path_output_profiles, 'dir'); mkdir(cal.path_output_profiles); end
    name = char(string(equipment_name));
    stem = strcat('F_', num2str(desired_focal_distance_ep), '_', name);

    % Plot the profile (and its sources, if derived) and save it
    fh = figure('Visible', 'off');
    hold on;
    if ~strcmp(type, 'measured')
        ok = find(isfinite(prof.source_depths(1, :)));
        for k = ok
            j = find(abs(available_foci_wrt_exit_plane - prof.source_depths(1, k)) < 1e-6, 1);
            plot(dist_from_exit_plane, intens_data(:, j), '-x', 'DisplayName', ...
                sprintf('Measurement, focus at %g (weight %.2f)', prof.source_depths(1, k), prof.weights(1, k)));
        end
    end
    plot(dist_from_exit_plane, norm_profile_focus, '-o', 'DisplayName', ...
        sprintf('%s, focus at %g', [upper(type(1)) type(2:end)], desired_focal_distance_ep));
    legend; xlabel('Distance wrt exit plane [mm]'); ylabel('Intensity [W/cm^2]');
    title(sprintf('Axial profile at %g mm from the exit plane (%s)', desired_focal_distance_ep, type));
    saveas(fh, fullfile(cal.path_output_profiles, ['Interpolation_at_' stem '.png']));
    close(fh);

    % A derived profile is what the calibration is fitted to: keep it
    if ~strcmp(type, 'measured')
        write_axial_profiles(prof, fullfile(cal.path_output_profiles, ['Profile_' type '_at_' stem '.csv']), ...
            struct('equipment', name));
    end

    % Determine the maximum intensity beyond the specified skip distance
    % for scaling to prevent catching max peak in near field peak.
    [~, closestIndex] = min(abs(dist_from_exit_plane - cal.skip_front_peak_mm));
    max_intens = max(norm_profile_focus(closestIndex:end));
end

function v = get_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
