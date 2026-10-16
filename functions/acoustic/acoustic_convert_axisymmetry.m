function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            acoustic_convert_axisymmetry(...
            parameters, sensor_data, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels, mode)
% ACOUSTIC_CONVERT_AXISYMMETRY  Expand axisymmetric simulation fields to 3D or the 2D plane
%
% Dispatches to convert_axisymmetric_to_3d or convert_axisymmetric_to_2d in
% the form chosen by AXISYMMETRIC_OUTPUT_MODE. Used after a fresh k-Wave run
% (ACOUSTIC_WRAPPER) and when an acoustic cache is reloaded (ACOUSTIC_STAGE),
% so both leave grid.dims, transducer positions and the medium in the same
% expanded state. Returns the inputs unchanged when the grid is not an
% axisymmetric half-plane.
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, ...
%    kwave_medium, kgrid, source, source_labels] = ...
%       acoustic_convert_axisymmetry(parameters, sensor_data, segmentation, ...
%                                    medium_masks, kwave_medium, kgrid, source, source_labels)
%   [...] = acoustic_convert_axisymmetry(..., mode)
%
% Input:
%   parameters   - PRESTUS config; grid.dims = [Naxial, Nr] on entry
%   sensor_data  - k-Wave sensor output with p_final, p_max_all [Pa]
%   segmentation - tissue label map [Naxial x Nr]
%   medium_masks - medium layer label map [Naxial x Nr]
%   kwave_medium - spatial medium property maps [Naxial x Nr] each
%   kgrid        - 2D axisymmetric kWaveGrid (replaced on output)
%   source       - k-Wave source with field p_mask [Naxial x Nr]
%   source_labels- transducer element label map [Naxial x Nr]
%   mode         - (optional) '3d' or '2d'; default AXISYMMETRIC_OUTPUT_MODE(parameters)
%
% Output:
%   All inputs returned with spatial fields expanded to [Nlateral x Nlateral x Naxial]
%   ('3d') or mirrored to [Nlateral x Naxial] ('2d'). Original bilateral dims
%   and positions restored from the axisym_bilateral_dims snapshot.
%
% See also: AXISYMMETRIC_OUTPUT_MODE, CONVERT_AXISYMMETRIC_TO_3D,
%           CONVERT_AXISYMMETRIC_TO_2D, ACOUSTIC_WRAPPER

arguments
    parameters   (1,1) struct
    sensor_data  (1,1) struct
    segmentation (:,:) {mustBeNumericOrLogical}
    medium_masks (:,:) {mustBeNumericOrLogical}
    kwave_medium (1,1) struct
    kgrid
    source       (1,1) struct
    source_labels(:,:) {mustBeNumericOrLogical}
    mode         (1,:) char = axisymmetric_output_mode(parameters)
end

    switch mode
        case '3d'
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_3d(sensor_data, parameters, segmentation, medium_masks, ...
                                           kwave_medium, kgrid, source, source_labels);
        case '2d'
            [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, medium_masks, ...
                                           kwave_medium, source, source_labels);
    end
end
