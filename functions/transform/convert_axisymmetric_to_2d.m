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
%   All inputs returned with spatial fields mirrored and transposed to [Nlateral x Naxial].
%   parameters.grid.dims updated to [Nlateral, Naxial] (original bilateral dims restored).
%   kgrid replaced with a 2D kWaveGrid.
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, GRID_AXISYMMETRY

arguments
    sensor_data  (1,1) struct
    parameters   (1,1) struct
    segmentation (:,:) {mustBeNumericOrLogical}
    medium_masks (:,:) {mustBeNumericOrLogical}
    kwave_medium (1,1) struct
    source       (1,1) struct
    source_labels(:,:) {mustBeNumericOrLogical}
end
%% Resolve original bilateral dimensions
assert(isfield(parameters.grid, 'axisym_bilateral_dims'), ...
    'convert_axisymmetric_to_2d: parameters.grid.axisym_bilateral_dims missing. Was grid_axisymmetry called?');

bilateral_dims = parameters.grid.axisym_bilateral_dims;   % [Naxial, Nlateral]
bilateral_tp   = parameters.transducer(1).trans_pos_bilateral;
bilateral_fp   = parameters.transducer(1).focus_pos_bilateral;

Nlateral  = bilateral_dims(2);   % target bilateral width (e.g. 280)
ax_center = bilateral_tp(2);     % axis column in the bilateral grid (e.g. 140)

% Mirror the half-plane about the axis to reconstruct exactly Nlateral
% columns, with the axis landing at column ax_center.
%
% Half-grid layout (Nr columns):
%   col 1         = r = 0  (axis)
%   col 2..Nr     = r = dy .. (Nr-1)*dy
%
% Target bilateral layout (Nlateral columns):
%   col 1..(ax_center-1)   = reflected cols ax_center..2  (left of axis)
%   col ax_center          = axis (col 1 of half-grid)
%   col (ax_center+1)..Nlateral = cols 2..(Nlateral-ax_center+1)  (right)
%
% The right half takes exactly Nlateral-ax_center columns from the
% half-grid (cols 2..(Nlateral-ax_center+1)), so the half-grid must have
% at least Nlateral-ax_center+1 columns.  grid_axisymmetry guarantees this
% because Nr = dims(2) - ax_center + 1 and dims(2) = Nlateral.
right_len = Nlateral - ax_center;               % columns right of axis

mirror = @(A) cat(2, ...
    fliplr(A(:, 2:ax_center)), ...              % reflected left:  ax_center-1 cols
    A(:, 1), ...                                % axis:            1 col
    A(:, 2:right_len+1));                       % right half:      right_len cols

sensor_data.p_final   = mirror(sensor_data.p_final);
sensor_data.p_max_all = mirror(sensor_data.p_max_all);
if isfield(sensor_data, 'p_complex')
    sensor_data.p_complex = mirror(real(sensor_data.p_complex)) + ...
                        1i .* mirror(imag(sensor_data.p_complex));
end
segmentation  = mirror(segmentation);
medium_masks  = mirror(medium_masks);
source.p_mask = mirror(source.p_mask);
source_labels = mirror(source_labels);
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = mirror(kwave_medium.(fields{i}));
    end
end

%% Restore original bilateral dims and positions
% The data is currently [Naxial x Nlateral].  Transpose back to
% [Nlateral x Naxial] to match the original bilateral orientation, then
% restore the stored dims and positions so downstream code sees the same
% grid it would from a non-axisymmetric run.
parameters.grid.dims = fliplr(bilateral_dims);   % [Nlateral, Naxial]

segmentation = segmentation';
medium_masks = medium_masks';
sensor_data.p_max_all = sensor_data.p_max_all';
sensor_data.p_final   = sensor_data.p_final';
if isfield(sensor_data, 'p_complex')
    sensor_data.p_complex = sensor_data.p_complex.';
end
source.p_mask = source.p_mask';
source_labels = source_labels';
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = kwave_medium.(fields{i})';
    end
end

% Restore positions from the bilateral snapshot: [lateral, axial]
parameters.transducer(1).trans_pos = fliplr(bilateral_tp);
parameters.transducer(1).focus_pos = fliplr(bilateral_fp);

% set up kgrid again for eventual heating sim
kgrid = kWaveGrid(parameters.grid.dims(1), parameters.grid.resolution_mm/1e3, ...
          parameters.grid.dims(2), parameters.grid.resolution_mm/1e3);

end