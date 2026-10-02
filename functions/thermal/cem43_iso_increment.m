function inc = cem43_iso_increment(T, T_max, dt_s)
% CEM43_ISO_INCREMENT  ISO thermal-dose (CEM43) increment for one time step
%
% Cumulative equivalent minutes at 43 °C (Sapareto & Dewey 1984), as used for
% the ISO CEM43 maps in THERMAL_SIMULATION:
%
%   CEM43 += dt/60 * R^(43 - T)
%
%   R = 0     T < 39 °C and the voxel has never reached 43 °C
%   R = 0.25  39 <= T < 43 °C and the voxel has never reached 43 °C
%   R = 0.5   T >= 43 °C, or the voxel has reached 43 °C before (T_max >= 43)
%
% R uses T_max so that a voxel that has once exceeded 43 °C keeps accumulating
% at R = 0.5 after cooling. A voxel that reaches 57 °C (irreversible damage) is
% assigned Inf. See doc_simulations-thermal.md.
%
% Use as:
%   inc = cem43_iso_increment(T, T_max, dt_s)
%
% Input:
%   T     - temperature [°C] at this step (array)
%   T_max - maximum temperature reached before this step [°C] (same size or scalar)
%   dt_s  - step duration [s]
%
% Output:
%   inc   - CEM43 increment [min], same size as T (Inf where T or T_max >= 57 °C)
%
% See also: THERMAL_SIMULATION, THERMAL_ANALYSIS

    R = 0.25 .* (T >= 39 & T < 43 & T_max < 43) + ...
        0.5  .* (T >= 43 | T_max >= 43);
    inc = dt_s ./ 60 .* R .^ (43 - T);
    inc(T >= 57 | T_max >= 57) = Inf;
end
