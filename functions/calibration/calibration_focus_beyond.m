function f = calibration_focus_beyond(I, z, skip_mm)
% CALIBRATION_FOCUS_BEYOND  FLHM centre of each axial profile, ignoring the near field
%
% A near-field peak would otherwise count as the focus, so only z >= skip_mm
% is used, as in FIT_TRANSDUCER_MODEL.
%
% Use as:
%   f = calibration_focus_beyond(I, z, skip_mm)
%
% Input:
%   I       - [N_z x N] profiles (columns)
%   z       - [N_z x 1] distance from the exit plane [mm]
%   skip_mm - near-field skip [mm]
%
% Output:
%   f - [1 x N] FLHM centre per column [mm from EP]
%
% See also: COMPUTE_FLHM_CENTER, CALIBRATION_EVALUATE

    k = z >= skip_mm;
    f = arrayfun(@(j) compute_flhm_center(I(k, j), z(k)), 1:size(I, 2));
end
