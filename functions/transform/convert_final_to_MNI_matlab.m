function [img_mni, final_to_mni_affine, mni_header] = convert_final_to_MNI_matlab(final_img, headreco_folder, inv_final_transformation_matrix, parameters, options)
% CONVERT_FINAL_TO_MNI_MATLAB  Transform a 3D image from subject simulation space to MNI space
%
% Chains three affine transforms: (1) simulation-space to original T1 conform
% space via inv_final_transformation_matrix, (2) T1 conform to MNI via the
% 12-DOF registration matrix in headreco_folder/toMNI/MNI2conform_12DOF.txt,
% and (3) MNI template voxel space via the MNI header affine. The result is
% resampled to the MNI template grid with nearest-neighbour interpolation and
% optionally saved as a compressed NIfTI. If the output file already exists and
% overwriting is disabled, the existing file is loaded instead.
%
% Use as:
%   [img_mni, final_to_mni_affine, mni_header] = ...
%       convert_final_to_MNI_matlab(final_img, headreco_folder, ...
%                                   inv_final_transformation_matrix, parameters)
%   [img_mni, final_to_mni_affine, mni_header] = ...
%       convert_final_to_MNI_matlab(..., 'nifti_filename', fname, 'nifti_data_type', 'uint8')
%
% Input:
%   final_img                       - 3D image in simulation (subject) space
%   headreco_folder                 - path to headreco m2m folder containing
%                                     toMNI/T1fs_nu_12DOF_MNI.nii.gz and T1fs_conform.nii.gz
%   inv_final_transformation_matrix - [4x3] inverse affine tform from simulation grid
%                                     to T1 conform space
%   parameters                      - PRESTUS config (overwrite settings)
%   options.check_nifti_on_disk     - whether to check/write disk file (default: 1)
%   options.nifti_filename          - output NIfTI path (default: '')
%   options.nifti_data_type         - output datatype (default: 'single')
%   options.BitsPerPixel            - bits per pixel (default: from MNI header)
%   options.fill_value              - fill value outside FOV (default: 0)
%
% Output:
%   img_mni            - image resampled to MNI template grid
%   final_to_mni_affine- [4x4] composite affine from simulation space to MNI voxels
%   mni_header         - niftiinfo header of the MNI template
%
% See also: SIMULATION_NIFTI, CONVERT_FINAL_TO_MNI_SIMNIBS, RAS_TO_GRID

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

    % Check if output NIfTI file exists and handle accordingly
    if options.check_nifti_on_disk
        img_mni_hdr = mni_header; % Use MNI header as template for output image
        img_mni_hdr.Datatype = options.nifti_data_type; % Set datatype for output image
        img_mni_hdr.BitsPerPixel = options.BitsPerPixel; % Set bits per pixel for output image

        % Confirm overwriting or load existing file
        if confirm_overwriting(options.nifti_filename, parameters)
            % Transform input image to MNI space and save as NIfTI file
            img_mni = affine_resample_3d(final_img, final_to_mni_affine', mni_template_size, 'nearest', options.fill_value);

            niftiwrite(img_mni, regexprep(options.nifti_filename, '.nii.gz$', ''), 'Compressed', true);
        else
            % Load existing NIfTI file if overwriting is not allowed
            img_mni = niftiread(options.nifti_filename);
        end
    end
    
end
