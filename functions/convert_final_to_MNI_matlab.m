function [img_mni, final_to_mni_affine, mni_header] = convert_final_to_MNI_matlab(final_img, headreco_folder, inv_final_transformation_matrix, parameters, options)
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

    mni_header = niftiinfo(fullfile(headreco_folder, 'toMNI','T1fs_nu_12DOF_MNI.nii.gz'));

    if isempty(options.BitsPerPixel)
        options.BitsPerPixel = mni_header.BitsPerPixel;
    end

    t1_orig_hdr = niftiinfo(fullfile(headreco_folder, 'T1fs_conform.nii.gz'));

    mni_template_size = mni_header.ImageSize;

    dlmopts = delimitedTextImportOptions('NumVariables',4);
    dlmopts.Delimiter = ' ';
    dlmopts.ConsecutiveDelimitersRule = 'join';
    dlmopts.LeadingDelimitersRule = 'ignore';
    mni2subj = readmatrix(fullfile(headreco_folder, 'toMNI','MNI2conform_12DOF.txt'),dlmopts, 'OutputType','double');
    subj2mni = inv(mni2subj);
    
    % a weird bug fix
    if length(inv_final_transformation_matrix)>1
        inv_final_transformation_matrix = inv_final_transformation_matrix(1);
    end
    final_to_mni_affine = inv(mni_header.Transform.T')*subj2mni*t1_orig_hdr.Transform.T'*inv_final_transformation_matrix.tdata.T';

    final_to_mni_tform = maketform('affine', final_to_mni_affine' );

    if options.check_nifti_on_disk
        img_mni_hdr = mni_header;
        img_mni_hdr.Datatype = options.nifti_data_type;
        img_mni_hdr.BitsPerPixel = options.BitsPerPixel;
        if confirm_overwriting(options.nifti_filename, parameters)
            img_mni = tformarray(final_img, final_to_mni_tform, ...
                        makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], mni_template_size, [], options.fill_value);

            niftiwrite(img_mni, regexprep(options.nifti_filename, '.nii.gz$', ''), 'Compressed', true);
        else
            img_mni = niftiread(options.nifti_filename);
        end
    end
    
end