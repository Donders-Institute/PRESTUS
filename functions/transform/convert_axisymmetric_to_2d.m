function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels)
% CONVERT_AXISYMMETRIC_TO_2D  Mirror axisymmetric simulation fields into a full 2D cross-section
%
% Reflects the [Nz x Nr] axisymmetric half-plane about the axis of symmetry
% to produce a [Nz x 2*Nr] bilateral cross-section, then transposes to
% [2*Nr x Nz]. Used when a follow-up 3D thermal simulation is not requested;
% for phantom simulations convert_axisymmetric_to_3d is always used instead.
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, ...
%    kwave_medium, kgrid, source, source_labels] = ...
%       convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, ...
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
%   All inputs returned with spatial fields mirrored and transposed to [2*Nr x Nz].
%   parameters.grid.dims updated to [2*Nr, Nz].
%   kgrid replaced with a 2D kWaveGrid.
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, GRID_AXISYMMETRY

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
sensor_data.p_final = cat(2, fliplr(sensor_data.p_final), sensor_data.p_final);
sensor_data.p_max_all = cat(2, fliplr(sensor_data.p_max_all), sensor_data.p_max_all);
segmentation = cat(2, fliplr(segmentation), segmentation);
medium_masks = cat(2, fliplr(medium_masks), medium_masks);
source.p_mask = cat(2, fliplr(source.p_mask), source.p_mask);
source_labels = cat(2, fliplr(source_labels), source_labels);
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = cat(2, fliplr(kwave_medium.(fields{i})), kwave_medium.(fields{i}));
    end
end

% shift positions to diameter location
trans_pos_final = parameters.transducer(1).trans_pos;
focus_pos_final = parameters.transducer(1).focus_pos;

trans_pos_final(2) = trans_pos_final(2)+parameters.grid.dims(2);
focus_pos_final(2) = focus_pos_final(2)+parameters.grid.dims(2);

% convert radial dimension size into diameter
parameters.grid.dims(2) = parameters.grid.dims(2)*2;

% shift radial to x dim
parameters.grid.dims = fliplr(parameters.grid.dims);
trans_pos_final = fliplr(trans_pos_final);
focus_pos_final = fliplr(focus_pos_final);
segmentation = segmentation';
medium_masks = medium_masks';
sensor_data.p_max_all = sensor_data.p_max_all';
sensor_data.p_final = sensor_data.p_final';
source.p_mask = source.p_mask';
source_labels = source_labels';
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = kwave_medium.(fields{i})';
    end
end

% Retain transducer and focus positions after all grid manipulations
parameters.transducer(1).trans_pos = trans_pos_final;
parameters.transducer(1).focus_pos = focus_pos_final;

% set up kgrid again for eventual heating sim
kgrid = kWaveGrid(parameters.grid.dims(1), parameters.grid.resolution_mm/1e3, ...
          parameters.grid.dims(2), parameters.grid.resolution_mm/1e3);

end