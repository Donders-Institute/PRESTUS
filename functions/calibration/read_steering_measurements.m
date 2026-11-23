function T = read_steering_measurements(file)
% READ_STEERING_MEASUREMENTS  Read hydrophone focal peaks of a steered matrix array
%
% One CSV row per steering target. Required columns:
%   target_x_mm, target_y_mm, target_z_mm  steering target, transducer frame
%                                          (origin at the apex / array centre,
%                                          beam along +z) [mm]
%   peak_pressure_mpa                      largest pressure amplitude in the
%                                          plane through the target, normal to
%                                          the beam axis [MPa]
% Optional columns: drive_v (drive voltage of that measurement [V]),
% peak_x_mm, peak_y_mm (position of the peak), model_pressure_mpa (a
% modelled peak at the same drive, used instead of MATRIX_FOCAL_PEAK).
% Further columns are kept as they are.
%
% Use as:
%   T = read_steering_measurements(file)
%
% Input:
%   file - path to the CSV file
%
% Output:
%   T - table with the columns above
%
% See also: BUILD_STEERING_CORRECTION, FIT_STEERING_CORRECTION

arguments
    file {mustBeTextScalar}
end

    T = readtable(file, 'VariableNamingRule', 'preserve');
    need = {'target_x_mm', 'target_y_mm', 'target_z_mm', 'peak_pressure_mpa'};
    miss = setdiff(need, T.Properties.VariableNames);
    if ~isempty(miss)
        error('read_steering_measurements:columns', '%s lacks column(s): %s.', file, strjoin(miss, ', '));
    end
    bad = ~isfinite(T.peak_pressure_mpa) | T.peak_pressure_mpa <= 0;
    if any(bad)
        warning('read_steering_measurements:dropped', '%d row(s) without a positive peak pressure dropped.', nnz(bad));
        T = T(~bad, :);
    end
end
