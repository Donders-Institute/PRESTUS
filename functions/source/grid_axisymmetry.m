function [parameters, segmentation, bone, medium_masks] = ...
    grid_axisymmetry(parameters, segmentation, bone, medium_masks)
    % adapt grid dimensions to axisymmetry (if requested)
    % grid should be specified as [axial, radial x 2]
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