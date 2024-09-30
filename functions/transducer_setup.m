function [transducer_mask, source_label, transducer_pars] = transducer_setup(transducer_pars, trans_pos, focus_pos, ... 
                                                             grid_dims, grid_step_mm)

% Function that creates a transducer mask based on the transducer parameters, its coordinates
% in the computational grid, the geometric focus coordinates and the computational 
% grid sizes. It also creates the label matrix that we can use to set the signal source.
% 
% 'transducer_pars' is a structure with fields defining the transducer parameters.
% All types of transducers have n_elements as a parameter for the number of
% transducer elements. The fields 'Elements_OD_mm' and 'Elements_ID_mm'
% define the outer and inner diameters of the elements

if ~isequal(size(focus_pos), size(trans_pos))
    error('Transducer and focus positions should be arrays of equal size')
end

if isequal(size(focus_pos), [3 1])
    focus_pos = focus_pos';
    trans_pos = trans_pos';
elseif ~isequal(size(focus_pos), [1 3])
	error('Transducer and focus positions should have the size [3 1] or [1 3]')
end

% Convert from [mm] to [grid points]
% and round the diameter to the nearest odd integer

transducer_pars.Elements_OD = 2*floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % [grid points]
transducer_pars.Elements_ID = 2*floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % [grid points]

transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm==0) = 0;

% Radius of curvature (ROC) indirect from the focal length or direct
% transducer_pars.radius_mm = sqrt((transducer_pars.Elements_OD_mm(2)/2)^2 + (focal_dist_mm)^2);   % [mm] % indirect
% transducer_pars.focal_dist_mm = focal_dist_mm;
% transducer_pars.focal_dist_based_on_focal_pos_mm = norm(trans_pos - focus_pos)*grid_step_mm;

% Convert from [mm] to [grid points]
transducer_pars.radius_grid = round(transducer_pars.curv_radius_mm / grid_step_mm);  % [grid points]

transducer_mask = zeros(grid_dims);
source_label = zeros(grid_dims);

% create element bowls one by one
% for the elements after the first, the inner bowl is subtracted
for el_i = 1:transducer_pars.n_elements
    bowl = makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos);
    if transducer_pars.Elements_ID(el_i) > 0
        bowl = bowl - makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos);
    end

    % Define the binary source mask
    transducer_mask = transducer_mask + bowl;

    % Label the elements for setting the source signal and to identify the elements in a figure
    source_label = source_label + el_i*bowl;
end

end

