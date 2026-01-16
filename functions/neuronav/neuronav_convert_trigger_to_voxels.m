function [transducer_ras, transducer_pos, target_ras, target_pos, t1_image] = neuronav_convert_trigger_to_voxels(sub_id, sides, outputStruct, parameters, pn)
% NEURONAV_CONVERT_TRIGGER_TO_VOXELS - Convert Localite trigger positions to voxel (image) coordinates.
%
% USAGE:
%   [transducers, targets, t1_image] = neuronav_convert_trigger_to_voxels(subject_id, sides, outputStruct, parameters, pn)
%
% INPUT:
%   sub_id       - String subject if (e.g., sub-001)
%   sides        - Array of stimulation sides (1:2 for left/right)
%   outputStruct - Struct containing averaged Localite (trigger) matrices for each side
%   parameters   - Simulation configuration struct, including hardware geometry
%   pn           - Struct containing path configuration for all data locations
%
% OUTPUT:
%   transducers  - Cell array with transducer positions (in voxel indices, xyz) for each side
%   targets      - Cell array with acoustic focus positions (in voxel indices, xyz) for each side
%   t1_image     - The loaded 3D T1 MRI volume as a matrix
%
% DESCRIPTION:
%   This function extracts averaged transformation matrices for each stimulation site,
%   determines the RAS (scanner/world) coordinates of (A) the transducer surface and (B) the ultrasound acoustic focus,
%   and converts these into index space (voxel coordinates) of the subject's MRI image.
%   The RAS coordinates are determined relative to the reference point and direction encoded in the trigger matrices.
%   Voxel space positions are indispensable for visualization, ROI extraction, and volumetric analysis.
%
% COORDINATE SYSTEM NOTE:
%   The initial Localite/trigger coordinates and all geometric calculations are in scanner RAS (Right-Anterior-Superior),
%   which reflects the native space of the T1 MRI scan. These are **not** in canonical MNI (Montreal Neurological Institute) space,
%   unless the subject's T1 is already warped to the MNI template.
%   - If you want MNI coordinates, you must apply a spatial normalization step (e.g., using SPM, FSL, or ANTs) to the T1 image and map
%     these points using the appropriate deformation fields or affine transformation to MNI template space.
%   - Output coordinates from this function are thus in **native scanner space / native MRI voxel space**, not MNI.
%
%   For group analysis or reporting in standard space, do spatial normalization post hoc.

    % --- STEP 1: Locate and load T1 MRI and header/info
    t1_file = dir(fullfile(pn.data_prelocalite, sprintf('%s_*T1*.nii.gz',sub_id)));
    t1_header = niftiinfo(fullfile(t1_file.folder, t1_file.name));
    t1_image = niftiread(fullfile(t1_file.folder, t1_file.name));

    % --- STEP 2: Compute the physical offset from the matrix origin to the transducer face
    reference_dist = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

    % --- STEP 3: For each stimulation side, extract the averaged transformation matrix and compute positions
    for i = sides
        % Extract 4x4 transformation from Matrix4D fields
        matrix_flat = struct2cell(outputStruct.TriggerMarker(i).Matrix4D);
        coord_matrix = reshape(cell2mat(matrix_flat)', [4, 4])';  % Transpose to match MATLAB layout

        % -- Extract and interpret Localite coordinate system:
        %   - coord_matrix(:,4): the position (origin) in RAS mm
        %   - coord_matrix(:,1): coil local X axis, typically pointing from coil center to head
        ref_pos = coord_matrix(:, 4);
        ref_vec = coord_matrix(:, 1);

        % -- Compute the RAS mm position of the transducer surface
        transducer_ras(i,:) = ref_pos + reference_dist * ref_vec;
        % -- Compute the RAS mm position of the acoustic focal point (forward along vector)
        target_ras(i,:) = ref_pos + parameters.expected_focal_distance_mm * ref_vec;

        % -- Convert these world (RAS mm) positions into MRI voxel index space
        transducer_pos(i,:) = ras_to_grid(transducer_ras(i,1:3)', t1_header);
        target_pos(i,:) = ras_to_grid(target_ras(i,1:3)', t1_header);
    end
    transducer_ras = round(transducer_ras(:,1:3));
    target_ras = round(target_ras(:,1:3));
end
