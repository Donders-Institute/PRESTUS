function density = fit_pairwiselinear(ct_data, hounsfieldUnits, massDensity, plot_fitting)
% FIT_PAIRWISELINEAR  Piecewise linear mapping from Hounsfield Units to mass density
%
% Maps CT values in Hounsfield units to mass density using piecewise
% linear interpolation between calibration points, with linear
% extrapolation beyond the calibration range.
%
% Use as:
%   density = fit_pairwiselinear(ct_data, hounsfieldUnits, massDensity)
%   density = fit_pairwiselinear(ct_data, hounsfieldUnits, massDensity, plot_fitting)
%
% Input:
%   ct_data         - CT data in Hounsfield units (any size)
%   hounsfieldUnits - [1xP] calibration HU values, sorted ascending
%   massDensity     - [1xP] corresponding density values [kg/m^3]
%   plot_fitting    - display the piecewise fit (default: false)
%
% Output:
%   density - mass density mapped from ct_data [kg/m^3], same size as ct_data
%
% See also: PCT_SKULLMAPPING, PCT_SOFT_TISSUE_PEAK

arguments
    ct_data         {mustBeNumeric}
    hounsfieldUnits (1,:) {mustBeNumeric}
    massDensity     (1,:) {mustBeNumeric}
    plot_fitting    (1,1) logical = false
end

% Ensure inputs are column vectors and sorted
hounsfieldUnits = sort(hounsfieldUnits(:));
massDensity = massDensity(:);
assert(length(hounsfieldUnits) == length(massDensity), ...
    'hounsfieldUnits and massDensity must have the same length');

% Precompute slopes and intercepts for each segment
n_points = length(hounsfieldUnits);
slopes = diff(massDensity) ./ diff(hounsfieldUnits);
intercepts = massDensity(1:end-1) - slopes .* hounsfieldUnits(1:end-1);

% Initialize output same type and size as input
density = zeros(size(ct_data), 'like', ct_data);

% Flatten ct_data for processing
ct_vec = ct_data(:);
density_vec = zeros(size(ct_vec), 'like', ct_vec);

% Apply piecewise linear mapping with extrapolation
for i = 1:n_points-1
    mask = (ct_vec >= hounsfieldUnits(i)) & (ct_vec < hounsfieldUnits(i+1));
    density_vec(mask) = slopes(i) * ct_vec(mask) + intercepts(i);
end

% Extrapolate left of first point
mask_left = ct_vec < hounsfieldUnits(1);
density_vec(mask_left) = slopes(1) * ct_vec(mask_left) + intercepts(1);

% Extrapolate right of last point
mask_right = ct_vec >= hounsfieldUnits(end);
density_vec(mask_right) = slopes(end) * ct_vec(mask_right) + intercepts(end);

% Reshape back to original shape
density = reshape(density_vec, size(ct_data));

% Plot fitting if requested
if plot_fitting
    plot_pairwise_fit(hounsfieldUnits, massDensity);
end

end

function plot_pairwise_fit(hounsfieldUnits, massDensity)
%PLOT_PAIRWISE_FIT Plot the piecewise linear fit and extrapolation lines.
figure;
hold on;

% Plot calibration points
plot(hounsfieldUnits, massDensity, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Calibration Points');

% Compute full range for plotting
hu_range = linspace(min(hounsfieldUnits) - 200, max(hounsfieldUnits) + 200, 1000);
density_range = fit_pairwise_eval(hu_range, hounsfieldUnits, massDensity);

plot(hu_range, density_range, 'b-', 'LineWidth', 2, 'DisplayName', 'Piecewise Linear Fit');

xlabel('Hounsfield Units (HU)');
ylabel('Mass Density (kg/m^3)');
title('Piecewise Linear HU to Density Mapping');
legend('Location', 'best');
grid on;
box on;
hold off;
end

function density = fit_pairwise_eval(ct_data, hounsfieldUnits, massDensity)
% Vectorized evaluation helper for plotting (standalone)
n_points = length(hounsfieldUnits);
slopes = diff(massDensity) ./ diff(hounsfieldUnits);
intercepts = massDensity(1:end-1) - slopes .* hounsfieldUnits(1:end-1);

density = zeros(size(ct_data));
ct_vec = ct_data(:);
density_vec = zeros(size(ct_vec));

for i = 1:n_points-1
    mask = (ct_vec >= hounsfieldUnits(i)) & (ct_vec < hounsfieldUnits(i+1));
    density_vec(mask) = slopes(i) * ct_vec(mask) + intercepts(i);
end

mask_left = ct_vec < hounsfieldUnits(1);
density_vec(mask_left) = slopes(1) * ct_vec(mask_left) + intercepts(1);

mask_right = ct_vec >= hounsfieldUnits(end);
density_vec(mask_right) = slopes(end) * ct_vec(mask_right) + intercepts(end);

density = reshape(density_vec, size(ct_data));
end
