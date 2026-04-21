function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels)
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
%   All inputs returned with spatial fields expanded to [2*Nr x 2*Nr x Nz].
%   parameters.grid.dims updated to [2*Nr, 2*Nr, Nz].
%   kgrid replaced with a 3D kWaveGrid.
%
% See also: RADIALEXPAND2DTO3D, GRID_AXISYMMETRY, CONVERT_AXISYMMETRIC_TO_2D

arguments
    sensor_data  (1,1) struct
    parameters   (1,1) struct
    segmentation (:,:) {mustBeNumeric}
    medium_masks (:,:) {mustBeNumeric}
    kwave_medium (1,1) struct
    kgrid        (1,1)
    source       (1,1) struct
    source_labels(:,:) {mustBeNumeric}
end

%% convert 2d axisymmetry images to 3d

% Continuous fields: linear interpolation preserves smooth spatial variation
sensor_data.p_final   = radialExpand2DTo3D(sensor_data.p_final);
sensor_data.p_max_all = radialExpand2DTo3D(sensor_data.p_max_all);
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}), 2) > 1
        kwave_medium.(fields{i}) = radialExpand2DTo3D(kwave_medium.(fields{i}));
    end
end

% Discrete label/mask fields: nearest-neighbour interpolation prevents
% fractional blended values at tissue boundaries
segmentation  = radialExpand2DTo3D(segmentation,           'nearest');
medium_masks  = radialExpand2DTo3D(medium_masks,           'nearest');
source.p_mask = radialExpand2DTo3D(double(source.p_mask),  'nearest');
source_labels = radialExpand2DTo3D(double(source_labels),  'nearest');

%% convert positions from 2d radial to 3d 

% shift positions to diameter location
trans_pos_final = parameters.transducer(1).trans_pos;
focus_pos_final = parameters.transducer(1).focus_pos;

trans_pos_final(2) = trans_pos_final(2)+parameters.grid.dims(2);
focus_pos_final(2) = focus_pos_final(2)+parameters.grid.dims(2);

% convert: radial x axial -> diameter x diameter x axial
Nr_full = parameters.grid.dims(2) * 2;
Nz      = parameters.grid.dims(1);
parameters.grid.dims = [Nr_full, Nr_full, Nz];

% convert: transducer and focus positions
trans_pos_final = fliplr(trans_pos_final);
trans_pos_final = [trans_pos_final(1), trans_pos_final(1), trans_pos_final(2)];
focus_pos_final = fliplr(focus_pos_final);
focus_pos_final = [focus_pos_final(1), focus_pos_final(1), focus_pos_final(2)];

% encode in the transducer position
parameters.transducer.trans_pos = trans_pos_final;
parameters.transducer.focus_pos = focus_pos_final;

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