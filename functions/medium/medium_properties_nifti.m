function medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, property)
    arguments
        parameters struct
        kwave_medium struct
        inv_final_transformation_matrix
        t1_header
        property % string of the tissue property field 
    end

    orig_hdr = t1_header; % header based on original T1w
    orig_hdr.Datatype = 'single';

    file_name = fullfile(char(parameters.io.cache_dir), char(property));
    
    if ~isfield(kwave_medium, property)
        warning('Missing field: %s', property);
    else
        % Transform and save if file doesn't exist
        % niftiwrite with 'Compressed',true produces file_name.nii.gz;
        % check for that explicitly to avoid repeated writes on re-runs.
        if ~isfile([file_name '.nii.gz']) && ~isfile([file_name '.nii'])
            transformed_data = single(tformarray(...
                kwave_medium.(property), ...
                inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), ...
                [1 2 3], [1 2 3], ...
                orig_hdr.ImageSize, [], 0));

            % Use two-step write: niftiwrite then manual gzip.
            % Compressed=true is unreliable on network filesystems because
            % MATLAB writes a temp .nii file first; if that write stalls,
            % gzip fails with "file does not exist".
            niftiwrite(transformed_data, file_name, orig_hdr);
            nii_file = [file_name '.nii'];
            if isfile(nii_file)
                gzip(nii_file);
                delete(nii_file);
            end
        end
    end
end
