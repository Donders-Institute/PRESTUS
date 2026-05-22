function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            acoustic_convert_axisymmetry(...
            parameters, sensor_data, segmentation, medium_masks, kwave_medium, source, source_labels)
% ACOUSTIC_CONVERT_AXISYMMETRY  Route axisymmetric simulation fields to 2D or 3D expansion
%
% Selects between radial 3D expansion and 2D bilateral mirroring based on
% whether a follow-up heating simulation is requested. When
% parameters.modules.run_heating_sims == 1, calls convert_axisymmetric_to_3d;
% otherwise calls convert_axisymmetric_to_2d. This is an older entry-point
% that predates the inline dispatch in acoustic_wrapper; prefer acoustic_wrapper
% for new code.
%
% Use as:
%   [sensor_data, parameters, segmentation, medium_masks, ...
%    kwave_medium, kgrid, source, source_labels] = ...
%       acoustic_convert_axisymmetry(parameters, sensor_data, segmentation, ...
%                                    medium_masks, kwave_medium, source, source_labels)
%
% Input:
%   parameters   - PRESTUS config; modules.run_heating_sims controls dispatch
%   sensor_data  - k-Wave sensor output with p_final, p_max_all [Pa]
%   segmentation - tissue label map [Nz x Nr]
%   medium_masks - medium layer label map [Nz x Nr]
%   kwave_medium - spatial medium property maps [Nz x Nr] each
%   source       - k-Wave source with field p_mask [Nz x Nr]
%   source_labels- transducer element label map [Nz x Nr]
%
% Output:
%   All inputs returned with spatial fields expanded to [Nlateral x Nlateral x Naxial]
%   (heating requested) or mirrored to [Nlateral x Naxial] (no heating).
%   Original bilateral dims restored from axisym_bilateral_dims snapshot.
%   parameters.grid.dims and kgrid updated to match.
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, CONVERT_AXISYMMETRIC_TO_2D, ACOUSTIC_WRAPPER

arguments
    parameters   (1,1) struct
    sensor_data  (1,1) struct
    segmentation (:,:) {mustBeNumericOrLogical}
    medium_masks (:,:) {mustBeNumericOrLogical}
    kwave_medium (1,1) struct
    source       (1,1) struct
    source_labels(:,:) {mustBeNumericOrLogical}
end

    % if using axisymmetric settings and requesting heating simulations, reshape output to 3D
    if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims==1
        % when follow-up thermal simulation is requested, expand axisymmetric to 3D
        [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            convert_axisymmetric_to_3d(...
            sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels);
    else
        % otherwise retain 2D output (axisymmetric expanded to 2D)
        [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            convert_axisymmetric_to_2d(...
            sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels);
    end
    
end