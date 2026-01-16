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

    file_name = fullfile(parameters.debug_dir, property);
    
    if ~isfield(kwave_medium, property)
        warning('Missing field: %s', property);
    else
        % Transform and save if file doesn't exist
        if ~isfile(file_name)
            transformed_data = single(tformarray(...
                kwave_medium.(property), ...
                inv_final_transformation_matrix, ...
                makeresampler('nearest', 'fill'), ...
                [1 2 3], [1 2 3], ...
                orig_hdr.ImageSize, [], 0));
            
            niftiwrite(transformed_data, file_name, orig_hdr, 'Compressed', true);
        end
    end
end
