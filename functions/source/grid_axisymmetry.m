function [parameters, segmentation, bone, medium_masks] = ...
    grid_axisymmetry(parameters, segmentation, bone, medium_masks)
% GRID_AXISYMMETRY  Adapt grid and tissue masks for k-Wave axisymmetric simulation
%
% When parameters.grid.axisymmetric == 1 and a 2D focus position is specified,
% the full bilateral grid is halved along the radial axis so that column 1
% of the returned arrays corresponds to r = 0 (the axis of symmetry), as
% required by kspaceFirstOrderAS. Transducer and focus positions are updated
% accordingly. Has no effect when axisymmetric mode is not requested.
%
% Use as:
%   [parameters, segmentation, bone, medium_masks] = ...
%       grid_axisymmetry(parameters, segmentation, bone, medium_masks)
%
% Input:
%   parameters   - (struct) PRESTUS config; must contain grid.dims [1x2],
%                    grid.axisymmetric, and transducer(1).trans_pos / focus_pos
%   segmentation - [Nx x Ny] numeric, tissue label map
%   bone         - [Nx x Ny] numeric, binary skull mask
%   medium_masks - [Nx x Ny] numeric, medium layer label map
%
% Output:
%   parameters   - updated: grid.dims halved along radial axis;
%                    trans_pos and focus_pos set to radial midline (r = 1)
%   segmentation - [Naxial x Nr] numeric, halved along radial axis
%   bone         - [Naxial x Nr] numeric, halved along radial axis
%   medium_masks - [Naxial x Nr] numeric, halved along radial axis
%
% See also: CONVERT_AXISYMMETRIC_TO_3D, CONVERT_AXISYMMETRIC_TO_2D

arguments
    parameters   (1,1) struct
    segmentation (:,:) {mustBeNumeric}
    bone         (:,:) {mustBeNumeric}
    medium_masks (:,:) {mustBeNumeric}
end
    if numel(parameters.transducer(1).focus_pos) == 2 && ...
            isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1
        if numel(parameters.transducer) > 1
            error('Axisymmetric simulations with multiple transducers are not supported (only a single transducer is allowed when axisymmetric == 1)');
        end
        trans_pos = parameters.transducer(1).trans_pos;
        focus_pos = parameters.transducer(1).focus_pos;
        % ensure that radial(y) dim is shorter than axial (x) dim
        if parameters.grid.dims(2) > parameters.grid.dims(1)
            parameters.grid.dims = fliplr(parameters.grid.dims);
            trans_pos = fliplr(trans_pos);
            focus_pos = fliplr(focus_pos);
            segmentation = segmentation';
            bone = bone';
            medium_masks = medium_masks';
        end
        % halve the grid along the radial axis
        % see http://www.k-wave.org/documentation/kspaceFirstOrderAS.php
        Ny_half = floor(parameters.grid.dims(2)/2);
        parameters.grid.dims(2) = Ny_half;
        segmentation = segmentation(:,Ny_half+1:end);
        bone = bone(:,Ny_half+1:end);
        medium_masks = medium_masks(:,Ny_half+1:end);
        % set transducer and focus position to the radial midline
        trans_pos(2) = 1; 
        focus_pos(2) = 1; 
        % Retain transducer and focus positions after grid manipulations
        parameters.transducer(1).trans_pos = trans_pos;
        parameters.transducer(1).focus_pos = focus_pos;
    end