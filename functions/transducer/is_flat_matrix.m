function tf = is_flat_matrix(tr_matrix)
% IS_FLAT_MATRIX  True for a flat matrix array (is_curved false or infinite radius)
%
% Use as:
%   tf = is_flat_matrix(tr_matrix)
%
% Input:
%   tr_matrix - tr.matrix: is_curved, curv_radius_mm [mm]
%
% Output:
%   tf - logical
%
% See also: MATRIX_ELEMENT_AIM, VALIDATE_MATRIX_TRANSDUCER

    tf = (isfield(tr_matrix, 'is_curved') && ~isempty(tr_matrix.is_curved) && ~tr_matrix.is_curved) || ...
         (isfield(tr_matrix, 'curv_radius_mm') && ~isempty(tr_matrix.curv_radius_mm) && isinf(tr_matrix.curv_radius_mm));
end
