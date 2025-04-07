function [img_mni, final_to_mni_affine, mni_header] = convert_final_to_MNI_matlab(final_img, headreco_folder, inv_final_transformation_matrix, parameters, options)

% CONVERT_FINAL_TO_MNI_MATLAB Converts an image from subject space to MNI space.
%
% This function transforms a 3D image (`final_img`) from subject-specific space 
% to MNI space using affine transformations derived from the `headreco_folder`. 
% The transformed image can be saved as a NIfTI file or loaded if it already exists.
%
% Input:
%   final_img                     - [Nx x Ny x Nz] matrix representing the 3D image in subject space.
%   headreco_folder               - String specifying the folder containing headreco outputs.
%   inv_final_transformation_matrix - [4x3] matrix representing the inverse transformation from final space to subject space.
%   parameters                    - Struct containing pipeline configuration parameters (e.g., overwrite settings).
%
% Options:
%   check_nifti_on_disk           - Boolean flag to check if the output NIfTI file exists on disk (default: 1).
%   nifti_filename                - String specifying the filename for the output NIfTI file (default: '').
%   nifti_data_type               - Datatype for the output NIfTI file (default: 'single').
%   BitsPerPixel                  - Bits per pixel for the output NIfTI file (default: taken from MNI header).
%   fill_value                    - Value used to fill empty regions during transformation (default: 0).
%
% Output:
%   img_mni                       - Transformed image in MNI space as a matrix.
%   final_to_mni_affine           - [4x4] affine transformation matrix mapping final space to MNI space.
%   mni_header                    - Header information for the MNI template image.

    arguments
        final_img(:,:,:)
        headreco_folder string
        inv_final_transformation_matrix (4, 3)
        parameters struct
        options.check_nifti_on_disk = 1
        options.nifti_filename = ''
        options.nifti_data_type = 'single'
        options.BitsPerPixel = []
        options.fill_value = 0
    end

    % Load header information for the MNI template image
    mni_header = niftiinfo(fullfile(headreco_folder, 'toMNI', 'T1fs_nu_12DOF_MNI.nii.gz'));

    % Set BitsPerPixel for the output image if not explicitly specified
    if isempty(options.BitsPerPixel)
        options.BitsPerPixel = mni_header.BitsPerPixel;
    end

    % Load header information for the original T1 conform image
    t1_orig_hdr = niftiinfo(fullfile(headreco_folder, 'T1fs_conform.nii.gz'));

    % Get dimensions of the MNI template image
    mni_template_size = mni_header.ImageSize;

    % Load transformation matrix from subject space to MNI space
    dlmopts = delimitedTextImportOptions('NumVariables', 4);
    dlmopts.Delimiter = ' ';
    dlmopts.ConsecutiveDelimitersRule = 'join';
    dlmopts.LeadingDelimitersRule = 'ignore';
    mni2subj = readmatrix(fullfile(headreco_folder, 'toMNI', 'MNI2conform_12DOF.txt'), dlmopts, 'OutputType', 'double');
    subj2mni = inv(mni2subj); % Compute inverse transformation matrix

    % Fix bug with transformation matrix length (ensure it's valid)
    if length(inv_final_transformation_matrix) > 1
        inv_final_transformation_matrix = inv_final_transformation_matrix(1);
    end

    % Compute affine transformation matrix from final space to MNI space
    final_to_mni_affine = inv(mni_header.Transform.T') * subj2mni * t1_orig_hdr.Transform.T' * inv_final_transformation_matrix.tdata.T';

    % Create affine transformation object for MATLAB's `tformarray`
    final_to_mni_tform = maketform('affine', final_to_mni_affine');

    % Check if output NIfTI file exists and handle accordingly
    if options.check_nifti_on_disk
        img_mni_hdr = mni_header; % Use MNI header as template for output image
        img_mni_hdr.Datatype = options.nifti_data_type; % Set datatype for output image
        img_mni_hdr.BitsPerPixel = options.BitsPerPixel; % Set bits per pixel for output image

        % Confirm overwriting or load existing file
        if confirm_overwriting(options.nifti_filename, parameters)
            % Transform input image to MNI space and save as NIfTI file
            img_mni = tformarray(final_img, final_to_mni_tform, ...
                        makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], mni_template_size, [], options.fill_value);

            niftiwrite(img_mni, regexprep(options.nifti_filename, '.nii.gz$', ''), 'Compressed', true);
        else
            % Load existing NIfTI file if overwriting is not allowed
            img_mni = niftiread(options.nifti_filename);
        end
    end
    
end
