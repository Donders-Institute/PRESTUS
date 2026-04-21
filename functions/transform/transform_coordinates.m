function output_pos = transform_coordinates(parameters, input_pos, input_cs, output_cs, nii_hdr)
% TRANSFORM_COORDINATES  Convert a 3D position between MNI, RAS+, and simulation grid spaces
%
% Supports three coordinate systems: 'mni' (MNI152 world space, mm),
% 'ras_plus' (subject T1 RAS+ world space, mm), and 'grid' (PRESTUS
% simulation voxel indices). MNI→subject conversion calls
% mni2subject_coords_LDfix (SimNIBS CLI); RAS+↔grid conversion uses
% the NIfTI affine via transformPointsInverse/Forward.
%
% Use as:
%   output_pos = transform_coordinates(parameters, input_pos, input_cs, output_cs)
%   output_pos = transform_coordinates(parameters, input_pos, input_cs, output_cs, nii_hdr)
%
% Input:
%   parameters - PRESTUS config (needs path.seg, subject_id,
%                startup.simnibs_bin_path, hpc.ld_library_path)
%   input_pos  - [1x3] position in the input coordinate system [mm or voxels]
%   input_cs   - input coordinate system: 'mni' | 'ras_plus' | 'grid'
%   output_cs  - target coordinate system: 'ras_plus' | 'grid'
%   nii_hdr    - niftiinfo header for RAS+↔grid conversions
%                (optional; mandatory when input_cs or output_cs is 'ras_plus'/'grid')
%
% Output:
%   output_pos - [1x3] position in the output coordinate system [mm or voxels]
%
% See also: MNI2SUBJECT_COORDS_LDFIX, RAS_TO_GRID, NIFTIINFO

    arguments
        parameters struct
        input_pos double
        input_cs string
        output_cs string
        nii_hdr struct = struct() % [Mandatory] for RAS+ and grid input
    end

    assert(numel(input_pos) == 3, 'input_pos must be [1x3]');
    
    switch input_cs
        case 'mni'
            m2m_path = fullfile(parameters.path.seg, sprintf('m2m_sub-%03i', parameters.subject_id));
            
            % MNI -> subject RAS+
            disp("Mapping MNI to subject RAS+...")
            ras_pos = mni2subject_coords_LDfix(input_pos, m2m_path, parameters);
            
            if strcmp(output_cs, 'ras_plus')
                output_pos = ras_pos;
            elseif strcmp(output_cs, 'grid')
                output_pos = transform_coordinates(parameters, ras_pos, 'ras_plus', 'grid', nii_hdr);
            else
                error('MNI output_cs: ''ras_plus'' or ''grid'' only');
            end
            
        case 'ras_plus'

            if strcmp(output_cs, 'grid')
                disp("Mapping RAS+ to subject grid...")
                output_pos = round(transformPointsInverse(nii_hdr.Transform, input_pos));
                % alternative:
                % output_pos = ras_to_grid(input_pos, nii_hdr);
            else
                error('RAS+ input, output_cs: ''grid'' only');
            end
            
        case 'grid'

            if strcmp(output_cs, 'ras_plus')
                % Voxel -> RAS+ (forward affine)
                disp("Mapping subject grid to subject RAS+...")
                output_pos = transformPointsForward(nii_hdr.Transform, input_pos);
                % alternative:
                % tmp = nii_hdr.Transform.T * [input_pos(:); 1];
                % output_pos = tmp(1:3)';
            else
                error('grid voxel input, output_cs: ''ras_plus'' only');
            end
            
        otherwise
            error('input_cs: ''mni'' | ''ras_plus'' | ''grid''');
    end
    
    fprintf('%s [%.1f %.1f %.1f] -> %s [%.1f %.1f %.1f]\n', ...
            input_cs, input_pos, output_cs, output_pos);
end