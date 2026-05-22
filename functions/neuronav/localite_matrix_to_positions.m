function [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
    localite_matrix_to_positions(coord_matrix, t1_header, parameters)
% LOCALITE_MATRIX_TO_POSITIONS  Convert a Localite 4×4 matrix to voxel positions.
%
% Shared position computation used by position_transducer_localite and
% neuronav_convert_trigger_to_voxels. Accepts a single 4×4 transformation
% matrix (typically the mean across a trigger series) and returns transducer
% and focus positions in both RAS (mm) and voxel space.
%
% Use as:
%   [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
%       localite_matrix_to_positions(coord_matrix, t1_header, parameters)
%
% Input:
%   coord_matrix - [4×4] Localite transformation matrix (RAS mm):
%                    column 4 = reference point position (translation)
%                    column 1 = direction vector (coil centre → head)
%   t1_header    - NIfTI header with Transform field (from niftiinfo)
%   parameters   - Simulation config struct. The offset from the Localite matrix
%                  origin to the bowl rear centre is derived in priority order from:
%                  (1) parameters.placement.localite.tracker_to_bowl_mm (explicit),
%                  (2) transducer geometry (-(curv_radius_mm − dist_geom_ep_mm)),
%                  (3) parameters.placement.localite.reference_distance_mm (fallback).
%
% Output:
%   trans_pos      - [1×3] bowl rear-centre position in voxel coordinates (ijk)
%                    (this is bowl_pos as passed to k-Wave's makeBowl)
%   focus_pos      - [1×3] acoustic focus position in voxel coordinates (ijk)
%   trans_pos_ras  - [1×3] bowl rear-centre position in RAS space (mm)
%   focus_pos_ras  - [1×3] acoustic focus position in RAS space (mm)
%
% See also: POSITION_TRANSDUCER_LOCALITE, NEURONAV_CONVERT_TRIGGER_TO_VOXELS

    arguments
        coord_matrix (4,4) double
        t1_header    (1,1) struct
        parameters   (1,1) struct
    end

    ref_pos = coord_matrix(:, 4);   % origin / reference position in RAS mm [4×1]
    ref_vec = coord_matrix(:, 1);   % direction: coil centre → head [4×1]

    % --- Offset: Localite matrix origin → bowl rear centre ---
    %
    % PRESTUS models the transducer as a bowl whose anchor point (trans_pos /
    % bowl_pos in k-Wave's makeBowl) is the centre of the rear surface of the
    % bowl — the geometric axis point at the back of the housing. The offset
    % below maps the Localite matrix origin (column 4) to that point.
    %
    % Priority order:
    %   1. parameters.placement.localite.tracker_to_bowl_mm  — explicit measured value
    %   2. Geometry-derived: -(curv_radius_mm - dist_geom_ep_mm)  — assumes the
    %      Localite matrix origin coincides with the exit plane
    %   3. parameters.placement.localite.reference_distance_mm  — legacy fallback
    %
    % A negative value means the bowl rear centre is on the housing/tracker side
    % (away from the head) relative to the matrix origin; positive means it is
    % further into the head.
    %
    % WARNING: The physical meaning of the Localite matrix origin is
    % instrument-registration-dependent and NOT validated in PRESTUS. It may
    % represent the IR tracker cluster centroid, the exit plane, or a housing
    % reference. Set tracker_to_bowl_mm to a directly measured value if spatial
    % accuracy is critical. See doc_placement_neuronav.md §"Localite matrix origin".

    lc = struct();
    if isfield(parameters, 'placement') && isfield(parameters.placement, 'localite')
        lc = parameters.placement.localite;
    end

    if isfield(lc, 'tracker_to_bowl_mm') && ~isempty(lc.tracker_to_bowl_mm)
        % Explicit calibrated value — highest priority
        reference_dist = lc.tracker_to_bowl_mm;
    else
        td      = parameters.transducer(1);
        td_type = '';
        if isfield(td, 'type'), td_type = td.type; end

        if ~isempty(td_type) && isfield(td, td_type) && ...
           isfield(td.(td_type), 'curv_radius_mm') && isfield(td.(td_type), 'dist_geom_ep_mm')
            % Assumes matrix origin = exit plane; bowl rear centre is behind it
            reference_dist = -(td.(td_type).curv_radius_mm - td.(td_type).dist_geom_ep_mm);
        elseif isfield(lc, 'reference_distance_mm') && ~isempty(lc.reference_distance_mm)
            reference_dist = lc.reference_distance_mm;
        else
            error('PRESTUS:localite:noReferenceDistance', ...
                ['Cannot determine offset from Localite matrix origin to bowl rear centre. ' ...
                 'Provide transducer.(type).curv_radius_mm and dist_geom_ep_mm, set ' ...
                 'parameters.placement.localite.tracker_to_bowl_mm (recommended), or set ' ...
                 'parameters.placement.localite.reference_distance_mm.']);
        end
    end

    % --- Focal distance (bowl rear centre → acoustic focus, along beam axis) ---
    % focal_distance_bowl is consistent with the bowl rear centre convention.
    parameters = focal_distance_calculation(parameters);
    focal_dist = parameters.transducer(1).focal_distance_bowl;

    % --- Compute positions in RAS space ---
    trans_pos_ras = (ref_pos + reference_dist * ref_vec)';   % [1×4] homogeneous
    focus_pos_ras = (ref_pos + focal_dist     * ref_vec)';

    % --- Convert RAS to voxel space ---
    trans_pos = ras_to_grid(trans_pos_ras(1:3)', t1_header);  % [1×3]
    focus_pos = ras_to_grid(focus_pos_ras(1:3)', t1_header);

    % Strip homogeneous coordinate from RAS outputs
    trans_pos_ras = trans_pos_ras(1:3);
    focus_pos_ras = focus_pos_ras(1:3);
end
