function a0_fit = fitPowerLawParamsMulti(a0, y, c0, f_ref, y_ref, plot_fit)
% FITPOWERLAWPARAMSMULTI  Fit power-law absorption prefactors for the k-Wave fractional Laplacian solver
%
% k-Wave's fractional Laplacian wave equation enforces a single global
% alpha_power (y_ref). When the true tissue power law a = a0 * f^y differs
% from this reference, or when attenuation values are large, the naive
% prefactor a0 must be rescaled so that the simulated absorption at f_ref
% matches the user-specified value. The rescaling follows Eq. 40 of Treeby
% & Cox (2014) and accounts for second-order dispersive corrections. If
% spatially varying y is provided (matrix a0 and y), a separate a0_fit is
% computed for each voxel.
%
% Use as:
%   a0_fit = fitPowerLawParamsMulti(a0, y, c0, f_ref, y_ref)
%   a0_fit = fitPowerLawParamsMulti(a0, y, c0, f_ref, y_ref, plot_fit)
%
% Input:
%   a0       - desired power-law absorption prefactors [dB/(MHz^y cm)]
%   y        - same size as a0, desired power-law exponents (0 <= y <= 3)
%   c0       - same size as a0, medium sound speed [m/s]
%   f_ref    - reference frequency [Hz]
%   y_ref    - reference power-law exponent (0 <= y_ref <= 3, must not equal 1)
%   plot_fit - display fit comparison plot (optional, default: false)
%
% Output:
%   a0_fit - same size as a0, rescaled prefactor for use as
%            medium.alpha_coeff with medium.alpha_power = y_ref
%
% Reference:
%   Treeby & Cox (2014), J. Acoust. Soc. Am. 136(4):1499-1510.
%   https://doi.org/10.1121/1.4894790
%
% See also: MEDIUM_SETUP, GET_ALPHA_COEFF

arguments
    a0       {mustBeNumeric, mustBeNonnegative}
    y        {mustBeNumeric}
    c0       {mustBeNumeric, mustBeNonnegative}
    f_ref    (1,1) {mustBeNumeric, mustBeNonnegative}
    y_ref    (1,1) {mustBeNumeric}
    plot_fit (1,1) logical = false
end

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

% check for plot input
if nargin < 6
    plot_fit = false;
end

% check inputs
validateattributes(a0,       {'numeric'}, {'real', 'nonnegative'},              mfilename, 'a0',       1);
validateattributes(y,        {'numeric'}, {'real', '>=', 0, '<=', 3},           mfilename, 'y',        2);
validateattributes(c0,       {'numeric'}, {'real', 'nonnegative'},              mfilename, 'c0',       3);
validateattributes(f_ref,    {'numeric'}, {'real', 'nonnegative', 'scalar'},    mfilename, 'f_ref',    4);
validateattributes(y_ref,    {'numeric'}, {'real', 'scalar', '>=', 0, '<=', 3}, mfilename, 'y_ref',    5);
validateattributes(plot_fit, {'logical'}, {'scalar'},                           mfilename, 'plot_fit', 6);

% make sure reference value isn't 1
if y_ref == 1
    error('Input for y_ref cannot be set to 1.');
end

% define frequency in rad/s
w = 2 * pi * f_ref;

% convert user defined a0 to Nepers/((rad/s)^y m)
a0_np = db2neper(a0, y);

% define desired absorption behaviour in Nepers/m
desired_absorption = a0_np .* w.^y;

% find the corresponding values of a0 that should be used in the
% fractional Laplacian wave equation to give the desired absorption
% behaviour taking into account second order effects (see Eq. 40 in [1])
a0_fit_np = desired_absorption ./ (  w.^y_ref + desired_absorption .* (y_ref + 1) .* c0 .* tan(pi .* y_ref ./ 2) .* w.^(y_ref - 1) );

% convert absorption prefactor back to dB/(MHz^y cm)
a0_fit = neper2db(a0_fit_np, y_ref);

% plot the final fit if desired
if plot_fit
    
    % create a small frequency range around the input frequency
    f_min = f_ref / 2;
    f_max = f_ref + f_min;
    f = f_min:(f_max - f_min)/1000:f_max;
    w = 2 * pi * f;
    
    % get suitable x-axis scale factor
    [~, scale, prefix] = scaleSI(f(end));
    
    % convert from Np/m to dB/cm
    conv_factor = (0.01 * 20 * log10(exp(1)));
    desired_absorption = desired_absorption .* conv_factor;
    
    % open figure
    figure;
    hold on;
    
    % get colors
    color_order = get(gca, 'ColorOrder');
    num_colors = size(color_order, 1);
    color_ind = 1;
    
    % get the unique values of the absorption coefficient
    [a0_uniq, ind_uniq] = unique(a0);
    
    % loop through the input values
    for plot_ind = 1:length(a0_uniq)
        
        % get the index of the value for comparison
        comp_ind = ind_uniq(plot_ind);
    
        % compute absorption behaviour
        absorption_fit = a0_fit_np(comp_ind) .* w.^y_ref ./ ( 1 - (y_ref + 1) .* a0_fit_np(comp_ind) .* c0(comp_ind) .* tan(pi .* y_ref ./ 2) .* w.^(y_ref - 1) );

        % convert from Np/m to dB/cm
        absorption_fit = absorption_fit .* conv_factor;

        % plot
        plot(f_ref .* scale, desired_absorption(comp_ind), '.', 'Color', color_order(color_ind, :));
        plot(f .* scale, a0(comp_ind) .* (f .* 1e-6) .^ y(comp_ind), '--', 'Color', color_order(color_ind, :));
        plot(f .* scale, absorption_fit, '-', 'Color', color_order(color_ind, :));

        % increment color index
        color_ind = color_ind + 1;
        if color_ind > num_colors
            color_ind = 1;
        end
        
    end
    
    % annotate plot
    set(gca, 'FontSize', 12);
    xlabel(['Frequency [' prefix 'Hz]']);
    ylabel('Absorption [dB/cm]');
    legend('Desired Absorption Value', 'Original Power Law', 'Fitted Power Law', 'Location', 'NorthWest');
    axis tight;
    
end