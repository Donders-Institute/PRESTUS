function aim_pos_m = matrix_element_aim(tr_matrix, elem_pos_m, trans_pos_m)
% MATRIX_ELEMENT_AIM  Point each matrix element faces (its normal), curved or flat
%
% Elements of a curved array face the natural focus, curv_radius_mm along
% the beam axis (+z) from the apex. Elements of a flat array (is_curved
% false, curv_radius_mm inf) face along +z; their aim point is set 1 m in
% front of each element, so that the vector element -> aim is the element
% normal for CREATE_MATRIX_KARRAY, kWaveArray.addDiscElement and the mask
% in TRANSDUCER_SETUP.
%
% Use as:
%   aim_pos_m = matrix_element_aim(tr_matrix, elem_pos_m, trans_pos_m)
%
% Input:
%   tr_matrix   - tr.matrix: is_curved, curv_radius_mm [mm]
%   elem_pos_m  - [3xN] element centres [m]
%   trans_pos_m - [3x1] transducer reference point (apex / array centre) [m]
%
% Output:
%   aim_pos_m - [3xN] point each element faces [m]
%
% See also: CREATE_MATRIX_KARRAY, TRANSDUCER_SETUP, IS_FLAT_MATRIX

    trans_pos_m = trans_pos_m(:);
    if is_flat_matrix(tr_matrix)
        aim_pos_m = elem_pos_m + [0; 0; 1];
    else
        natural_focus_pos_m = trans_pos_m + [0; 0; tr_matrix.curv_radius_mm / 1000];
        aim_pos_m = repmat(natural_focus_pos_m, 1, size(elem_pos_m, 2));
    end
end
