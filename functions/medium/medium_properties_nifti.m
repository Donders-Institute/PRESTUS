function medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, property)
% MEDIUM_PROPERTIES_NIFTI  Transform a simulation-grid medium property map to T1 space and save as NIfTI
%
% Applies the inverse affine transformation (inv_final_transformation_matrix)
% to the requested field of kwave_medium using nearest-neighbour resampling,
% then writes the result as a gzip-compressed NIfTI file in
% parameters.io.debug_dir_medium. The output filename is <debug_dir_medium>/<property>_t1.nii.gz.
% To avoid repeated writes on re-runs the function skips the write if the
% file already exists. Uses a two-step write (niftiwrite then gzip) rather
% than Compressed=true to work reliably on network filesystems.
%
% Use as:
%   medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, property)
%
% Input:
%   parameters                      - PRESTUS config; must contain io.debug_dir_medium
%   kwave_medium                    - medium property maps; must contain a field named property
%   inv_final_transformation_matrix - affine inverse transform from simulation grid
%                                     back to T1 space (as expected by tformarray)
%   t1_header                       - NIfTI header of the original T1 image
%                                     (ImageSize used as output dimensions)
%   property                        - field name of the property to export
%                                     (e.g., 'sound_speed', 'density')
%
% See also: MEDIUM_SETUP, NIFTIWRITE

arguments
    parameters                     (1,1) struct
    kwave_medium                   (1,1) struct
    inv_final_transformation_matrix
    t1_header
    property                       (1,:) char
end

    orig_hdr = t1_header; % header based on original T1w
    orig_hdr.Datatype = 'single';

    file_name = fullfile(char(parameters.io.debug_dir_medium), [char(property) '_t1']);
    
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
