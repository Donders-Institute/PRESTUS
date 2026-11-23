function [cls, table, water, info] = babelbrain_medium_table(kwave_medium, water_props, freq_hz, max_materials)
% BABELBRAIN_MEDIUM_TABLE  Quantise the PRESTUS acoustic medium into BabelBrain materials
%
% For babelbrain.medium 'prestus': BabelViscoFDTD runs on the medium the
% k-Wave backend sees (smoothed properties, pseudo-CT skull) instead of
% BabelBrain's five tissue labels. Each voxel gets the class of its
% (density, sound speed, attenuation) triple, quantised on a grid that is
% coarsened until at most max_materials classes remain. A class carries the
% mean properties of its voxels. Voxels with the water properties are class 0.
%
% Use as:
%   [cls, table, water, info] = babelbrain_medium_table(kwave_medium, water_props, freq_hz)
%   [...] = babelbrain_medium_table(..., max_materials)
%
% Input:
%   kwave_medium  - struct with sound_speed [m/s], density [kg/m^3] and
%                   alpha_coeff [dB/(cm MHz^y)] (maps or scalars, same grid)
%                   and a scalar alpha_power y (medium_setup output)
%   water_props   - parameters.medium_properties.water (sound_speed, density,
%                   alpha_coeff; attenuation at freq_hz with kwave_medium.alpha_power)
%   freq_hz       - drive frequency [Hz]
%   max_materials - class limit (default 1000; BabelBrain's CT path uses up to 1023)
%
% Output:
%   cls   - uint16 class map, size of the medium maps; 0 = water, 1..K
%   table - [K x 3] class properties [density kg/m^3, sound speed m/s,
%           attenuation Np/m at freq_hz]
%   water - [1 x 3] water properties in the same columns
%   info  - struct: steps [kg/m^3, m/s, Np/m] used, n_classes, max_abs_error
%           (largest voxel deviation from its class mean, per column)
%
% See also: WRITE_BABELBRAIN_INPUT, BABELBRAIN_WRAPPER

arguments
    kwave_medium  (1,1) struct
    water_props   (1,1) struct
    freq_hz       (1,1) double {mustBePositive}
    max_materials (1,1) double {mustBeInteger, mustBePositive} = 1000
end

DB_CM_TO_NP_M = 100 / 8.685889638;   % as run_babelbrain_case.py
y     = double(kwave_medium.alpha_power);
f_MHz = freq_hz / 1e6;
c   = double(kwave_medium.sound_speed);
rho = double(kwave_medium.density);
att = double(kwave_medium.alpha_coeff) * f_MHz^y * DB_CM_TO_NP_M;
sz = size(c);
if isscalar(rho); rho = repmat(rho, sz); end
if isscalar(att); att = repmat(att, sz); end
if ~isequal(size(rho), sz) || ~isequal(size(att), sz)
    error('babelbrain_medium_table:size', 'sound_speed, density and alpha_coeff must share one grid.');
end

water_att = 0;
if isfield(water_props, 'alpha_coeff') && ~isempty(water_props.alpha_coeff)
    water_att = double(water_props.alpha_coeff) * f_MHz^y * DB_CM_TO_NP_M;
end
water = [double(water_props.density), double(water_props.sound_speed), water_att];

X = [rho(:), c(:), att(:)];
is_water = abs(X(:,1) - water(1)) < 0.5 & abs(X(:,2) - water(2)) < 0.5 & abs(X(:,3) - water(3)) < 0.01;
tissue = find(~is_water);

cls = zeros(sz, 'uint16');
steps = [1, 1, 0.1];   % kg/m^3, m/s, Np/m
table = zeros(0, 3);
max_err = [0 0 0];
if ~isempty(tissue)
    Xt = X(tissue, :);
    while true
        [~, ~, idx] = unique(round(Xt ./ steps), 'rows');
        if max(idx) <= max_materials; break; end
        steps = steps * 2;
    end
    K = max(idx);
    table = [accumarray(idx, Xt(:,1), [K 1], @mean), ...
             accumarray(idx, Xt(:,2), [K 1], @mean), ...
             accumarray(idx, Xt(:,3), [K 1], @mean)];
    cls(tissue) = uint16(idx);
    max_err = max(abs(Xt - table(idx, :)), [], 1);
end

info = struct('steps', steps, 'n_classes', size(table, 1), 'max_abs_error', max_err);
end
