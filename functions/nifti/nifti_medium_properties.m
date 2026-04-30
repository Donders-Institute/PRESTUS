function nifti_medium_properties(parameters, planimg, kwave_medium, m2m_folder)
% NIFTI_MEDIUM_PROPERTIES  Write medium acoustic/thermal property maps to NIfTI
%
% Writes per-property NIfTIs back-transformed to T1 space into
% dir_output/ap/. Filenames inherit the run's output_affix. MNI
% conversions follow parameters.io.save_MNI and are written to the same
% folder with 'MNI' encoded in the filename.
%
% Only called when should_save_output(parameters.io, 'save_property_maps')
% returns true (caller's responsibility).
%
% Use as:
%   nifti_medium_properties(parameters, planimg, kwave_medium, m2m_folder)
%
% Input:
%   parameters   - PRESTUS config struct (must be layered medium to reach here)
%   planimg      - struct with fields: t1_image_orig, inv_transf, t1_header
%   kwave_medium - struct with medium property fields
%   m2m_folder   - path to SimNIBS m2m folder
%
% See also: NIFTI_MEDIUM, NIFTI_TO_T1W, NIFTI_TO_MNI, SHOULD_SAVE_OUTPUT

properties = {'sound_speed', 'density', 'alpha_coeff', 'alpha_power', ...
              'thermal_conductivity', 'specific_heat', 'perfusion_coeff', 'absorption_fraction'};

try
    out_dir = fullfile(char(parameters.io.dir_output), 'ap');
    if ~isfolder(out_dir); mkdir(out_dir); end

    affix = parameters.io.output_affix;

    for prop = properties
        prop = char(prop); %#ok<FXSET>
        if ~isfield(kwave_medium, prop)
            continue
        end
        t1w_file = fullfile(out_dir, sprintf('%s_T1w%s', prop, affix));
        mni_file = fullfile(out_dir, sprintf('%s_MNI%s.nii.gz', prop, affix));

        nifti_to_t1w(single(kwave_medium.(prop)), t1w_file, parameters, planimg, ...
            'Resampler', 'nearest');
        nifti_to_mni(strcat(t1w_file, '.nii.gz'), mni_file, parameters, true, m2m_folder);
    end
catch ME
    prev = warning('off', 'backtrace');
    warn('nifti_medium_properties:writeError', ...
        'Could not write medium property NIfTIs: %s', ME.message);
    warning(prev);
end
end
