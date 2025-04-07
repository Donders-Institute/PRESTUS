function [trans_pos, focus_pos, focus_pos_ras, trans_pos_ras] = position_transducer_localite(localite_instr_file, mri_hdr, parameters)

% POSITION_TRANSDUCER_LOCALITE Determines transducer and focus positions in voxel space using Localite data.
%
% This function uses instrument marker data from a Localite XML file and an MRI header to calculate 
% the transducer position and focus point in voxel (ijk) space. The Localite coordinates are provided 
% in RAS format and are converted to voxel space using the transformation matrix from the MRI header.
%
% Input:
%   localite_instr_file - Path to the Localite instrument markers XML file.
%   mri_hdr             - Header of the MRI file containing the world-to-voxel transformation matrix.
%   parameters          - Struct containing additional parameters:
%                         * reference_transducer_distance_mm: Distance between reference point and transducer (in mm).
%                         * expected_focal_distance_mm: Distance between transducer and focus point (in mm).
%
% Output:
%   trans_pos           - [1x3] array specifying the transducer position in voxel coordinates (ijk).
%   focus_pos           - [1x3] array specifying the focus position in voxel coordinates (ijk).
%   trans_pos_ras       - [1x4] array specifying the transducer position in RAS coordinates.
%   focus_pos_ras       - [1x4] array specifying the focus position in RAS coordinates.

    %% Notes on Inputs and Coordinate Systems
    % - 'localite_instr_file' contains the path to the Localite instrument markers XML file.
    % - 'mri_hdr' is the header of the MRI file that provides world-to-voxel transformation matrix.
    % - 'reference_transducer_distance_mm' is the distance between neuronavigation reference point and transducer (mm).
    % - 'expected_focal_distance_mm' is the distance between transducer and focus point (mm).
    %
    % The Localite coordinates are provided in RAS format. For details on NIFTI coordinate systems, refer to:
    % "Orientation information" at https://brainder.org/2012/09/23/the-nifti-file-format/.
    
    %% Parse Localite XML File
    % Read instrument marker data from Localite XML file using `xml2struct`.
    xml = xml2struct(localite_instr_file);

    % Extract 4x4 coordinate matrix from XML data.
    coord_xml = xml.Children(4).Children(2).Children(2);  % Navigate XML hierarchy to marker data
    coord_matrix = {coord_xml.Attributes.Value};          % Extract coordinate values as strings

    %% Convert Coordinate Matrix
    % Convert string values to numeric array and reshape into a 4x4 matrix.
    coord_matrix = str2double(coord_matrix);
    coord_matrix = reshape(coord_matrix', [4, 4])';

    %% Extract Reference Position and Direction Vectors
    % Reference position: translation vector from neuronavigation reference point.
    reference_pos = coord_matrix(:, 4); 

    % Direction vector pointing from coil center towards head.
    reference_center_to_head = coord_matrix(:, 1);

    %% Compute Transducer and Focus Positions in RAS Space
    % Calculate transducer position relative to reference point.
    trans_pos_ras = reference_pos + parameters.reference_transducer_distance_mm * reference_center_to_head;

    % Calculate focus position relative to transducer position.
    focus_pos_ras = trans_pos_ras + parameters.expected_focal_distance_mm * reference_center_to_head;

    %% Convert Positions from RAS Space to Voxel Space
    % Use MRI header's transformation matrix to convert RAS coordinates to voxel coordinates.
    trans_pos = round(mri_hdr.Transform.T' \ trans_pos_ras);  % Transducer position in voxel space.
    focus_pos = round(mri_hdr.Transform.T' \ focus_pos_ras);  % Focus position in voxel space.

    % Extract x, y, z components of positions (ignore homogeneous coordinate).
    trans_pos = trans_pos(1:3);
    focus_pos = focus_pos(1:3);

    %% Display Results
    disp('The transducer position is taken from the Localite files, and the expected focus position is set based on the expected focal distance.');
    disp('The transducer position & the focus position:');
    disp([trans_pos, focus_pos]);

end