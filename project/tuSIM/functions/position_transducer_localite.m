function [trans_pos, focus_pos, focus_pos_ras, trans_pos_ras] = position_transducer_localite(localite_instr_file, mri_hdr, ...
    parameters)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                         Position transducer                       %
    %                                                                   %
    % Position_transducer uses the data from the Localite instrument    %
    % markers and MRI scan header to determine the transducer coordinate%
    % in voxel (ijk) space.                                             %
    %                                                                   %
    % Some notes:                                                       %
    % - 'localite_instr_file' is the path to the Localite instrument    %
    % markers XML file.                                                 %
    % - 'mri_hdr' is the header of the MRI file that has the desired    %
    % voxel space (used to extract world to voxel transformation matrix)%
    % - 'reference_transducer_distance' is the distance between the     %
    % reference used for the neuronavigation and the transducer, mm.    %
    % - 'focal_distance' is the distance between the transducer and the %
    % focus point, mm.                                                  %
    %                                                                   %
    % The Localite coordinates are in RAS format. For the helpful info  %
    % on the NIFTI coordinates and headers, see the section             %
    % "Orientation information" here:                                   %
    % https://brainder.org/2012/09/23/the-nifti-file-format/            %
    %                                                                   %
    % The code for reading Localite data is written by Judith Lefkes    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Read in the instrument marker data to get the transducer coordinates
    % Example file of one of the participants from LocaliteData folder
    % localite_instr_file = "../data/example_transducer_Localite.xml";
    xml = xml2struct(localite_instr_file);
    
    % Instrument markers are saved in the 4x4 matrix.
    % Read in the 4D matrix of Instrument marker position
    coord_xml = xml.Children(4).Children(2).Children(2);
    coord_matrix = {coord_xml.Attributes.Value};
    
    % First convert into array then reshape
    coord_matrix = str2double(coord_matrix);
    coord_matrix = reshape(coord_matrix',[4,4])';
     
    % The x axis goes down vertically from the center of the coil housing to the
    % head; y axis goes in the direction of the handle and the z axis passes from left to
    % right, starting at the center of the coil housing (calibrated ?hotspot?) at a right 
    % angle to the x axis.
    
    % Converting into 3D vectors
    reference_pos = coord_matrix(:,4); % Position of the reference
    reference_center_to_head = coord_matrix(:,1); % From the center of the coil towards head
    % y = coord_matrix(:,2); % From the center of the coil to the left coil handle
    % z = coord_matrix(:,3); % From the center of the coil to the right coil handle
    
    trans_pos_ras = reference_pos + parameters.reference_transducer_distance_mm*reference_center_to_head; 
    focus_pos_ras = trans_pos_ras + parameters.expected_focal_distance_mm*reference_center_to_head;
    
    trans_pos = round(mri_hdr.Transform.T' \ trans_pos_ras);
    focus_pos = round(mri_hdr.Transform.T' \ focus_pos_ras);
    
    trans_pos = trans_pos(1:3);
    focus_pos = focus_pos(1:3);
    disp('The transducer position is taken from the localite files, the expected focus position is set based on the expected focal distance')
    disp('The transducer position & the focus position:')
    [trans_pos focus_pos]
    
    % actual_focus = t+35.1*x; 

end