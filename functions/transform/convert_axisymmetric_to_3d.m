function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels)
% CONVERT_AXISYMMETRIC_TO_3D  Expand all axisymmetric simulation fields from 2D to 3D
%
% Applies radial expansion to every spatial field produced by a k-Wave
% axisymmetric simulation (kspaceFirstOrderAS), converting the [Nz x Nr]
% half-plane representation into a full [2*Nr x 2*Nr x Nz] Cartesian volume.
% Continuous fields (pressure, medium properties) use linear interpolation;
% discrete label/mask fields use nearest-neighbour to preserve integer values.
% Transducer and focus positions and grid.dims are updated to match.
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, ...
%    kwave_medium, kgrid, source, source_labels] = ...
%       convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, ...
%                                  medium_masks, kwave_medium, source, source_labels)
%
% Input:
%   sensor_data  - (struct) k-Wave sensor output with fields p_final, p_max_all [Pa]
%   parameters   - (struct) PRESTUS config; grid.dims = [Nz, Nr] on entry
%   segmentation - [Nz x Nr] uint8, tissue label map
%   medium_masks - [Nz x Nr] uint8, medium layer label map
%   kwave_medium - (struct) spatial medium property maps [Nz x Nr] each
%   kgrid        - (kWaveGrid) 2D axisymmetric grid
%   source       - (struct) k-Wave source with field p_mask [Nz x Nr]
%   source_labels- [Nz x Nr] numeric, transducer element label map
%
% Output:
%   All inputs returned with spatial fields expanded to [Nlateral x Nlateral x Naxial].
%   parameters.grid.dims updated to [Nlateral, Nlateral, Naxial] (original bilateral dims restored).
%   kgrid replaced with a 3D kWaveGrid.
%
% See also: RADIALEXPAND2DTO3D, GRID_AXISYMMETRY, CONVERT_AXISYMMETRIC_TO_2D

arguments
    sensor_data  (1,1) struct
    parameters   (1,1) struct
    segmentation (:,:) {mustBeNumericOrLogical}
    medium_masks (:,:) {mustBeNumericOrLogical}
    kwave_medium (1,1) struct
    kgrid        (1,1)
    source       (1,1) struct
    source_labels(:,:) {mustBeNumericOrLogical}
end

%% Resolve original bilateral dimensions
% grid_axisymmetry saved the bilateral state (post-axis-flip, pre-halve).
% Restore to those dimensions so the 3D output looks identical to what a
% full bilateral simulation would have produced.
assert(isfield(parameters.grid, 'axisym_bilateral_dims'), ...
    'convert_axisymmetric_to_3d: parameters.grid.axisym_bilateral_dims missing. Was grid_axisymmetry called?');

bilateral_dims  = parameters.grid.axisym_bilateral_dims;   % [Naxial, Nlateral]
bilateral_tp    = parameters.transducer(1).trans_pos_bilateral;   % [axial, lateral_center]
bilateral_fp    = parameters.transducer(1).focus_pos_bilateral;

Naxial   = bilateral_dims(1);
Nlateral = bilateral_dims(2);          % original full bilateral width
ax_center = bilateral_tp(2);           % axis column index in the bilateral grid

%% convert 2d axisymmetry images to 3d
% Expand to the original Nlateral × Nlateral cross-section so the output
% dimensions match a non-axisymmetric run exactly.
expand = @(A, m) radialExpand2DTo3D(A, m, Nlateral, ax_center);

% Continuous fields: linear interpolation
sensor_data.p_final   = expand(sensor_data.p_final,   'linear');
sensor_data.p_max_all = expand(sensor_data.p_max_all, 'linear');
if isfield(sensor_data, 'p_complex')
    sensor_data.p_complex = expand(real(sensor_data.p_complex), 'linear') + ...
                        1i .* expand(imag(sensor_data.p_complex), 'linear');
end
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}), 2) > 1
        kwave_medium.(fields{i}) = expand(kwave_medium.(fields{i}), 'linear');
    end
end

% Discrete label/mask fields: nearest-neighbour to preserve integer values
segmentation  = expand(segmentation,           'nearest');
medium_masks  = expand(medium_masks,           'nearest');
source.p_mask = expand(double(source.p_mask),  'nearest');
source_labels = expand(double(source_labels),  'nearest');

%% Restore original dims and positions
% Output is [Nlateral x Nlateral x Naxial], axis at (ax_center, ax_center).
parameters.grid.dims = [Nlateral, Nlateral, Naxial];

trans_pos_final = [ax_center, ax_center, bilateral_tp(1)];
focus_pos_final = [ax_center, ax_center, bilateral_fp(1)];

parameters.transducer(1).trans_pos = trans_pos_final;
parameters.transducer(1).focus_pos = focus_pos_final;

%% assert all expanded fields have consistent size
expected_size = size(sensor_data.p_max_all);
assert(isequal(size(sensor_data.p_final),  expected_size), 'convert_axisymmetric_to_3d: p_final size mismatch');
assert(isequal(size(segmentation),         expected_size), 'convert_axisymmetric_to_3d: segmentation size mismatch');
assert(isequal(size(medium_masks),         expected_size), 'convert_axisymmetric_to_3d: medium_masks size mismatch');
assert(isequal(size(source.p_mask),        expected_size), 'convert_axisymmetric_to_3d: source.p_mask size mismatch');
assert(isequal(size(source_labels),        expected_size), 'convert_axisymmetric_to_3d: source_labels size mismatch');

%% set up 3D kgrid for follow-up heating simulations
kgrid = kWaveGrid(parameters.grid.dims(1), parameters.grid.resolution_mm/1e3, ...
          parameters.grid.dims(2), parameters.grid.resolution_mm/1e3, ...
          parameters.grid.dims(3), parameters.grid.resolution_mm/1e3);

%% DEBUG: visually inspect

% % plot center of transducer
% figure; imagesc(squeeze(sensor_data.p_max_all(:,floor(parameters.grid.dims(1)/2),:)))
% 
% % plot each slice along the axial dimension
% figure;
% for index = 1:parameters.grid.dims(3)
%     imagesc(squeeze(sensor_data.p_max_all(:,:,index)))
%     pause(0.5)
% end

end