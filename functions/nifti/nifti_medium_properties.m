function nifti_medium_properties(parameters, planimg, kwave_medium, m2m_folder)
% NIFTI_MEDIUM_PROPERTIES  Write medium acoustic/thermal property maps to NIfTI
%
% Writes per-property NIfTIs back-transformed to T1 space into
% dir_output/nii/properties/. Filenames inherit the run's output_affix. MNI
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

% Water baseline fill values for out-of-FOV voxels (values outside the
% simulation grid should reflect the water medium, not zero).
w = parameters.medium_properties.water;
water_fill = struct( ...
    'sound_speed',         w.sound_speed, ...
    'density',             w.density, ...
    'alpha_coeff',         w.alpha_coeff, ...
    'alpha_power',         w.alpha_power, ...
    'thermal_conductivity', w.thermal_conductivity, ...
    'specific_heat',       w.specific_heat_capacity, ...
    'perfusion_coeff',     (w.perfusion / 60) * w.density * 1e-6, ...
    'absorption_fraction', w.absorption_fraction);

properties = fieldnames(water_fill);

try
    out_dir = fullfile(char(parameters.io.dir_output), 'nii', 'properties');
    if ~isfolder(out_dir); mkdir(out_dir); end

    affix = parameters.io.output_affix;

    % Write all T1w property maps
    for i = 1:numel(properties)
        prop = properties{i};
        if ~isfield(kwave_medium, prop)
            continue
        end
        t1w_file = fullfile(out_dir, sprintf('%s_T1w%s', prop, affix));

        prop_data = single(kwave_medium.(prop));
        if isscalar(prop_data)
            prop_data = repmat(prop_data, size(kwave_medium.sound_speed));
        end

        nifti_to_t1w(prop_data, t1w_file, parameters, planimg, ...
            'Resampler', 'nearest', ...
            'FillValue', water_fill.(prop));
    end

    % Batch all MNI conversions into a single parallel subject2mni call
    if should_save_output(parameters.io, 'save_MNI')
        mni_in    = {};
        mni_out   = {};
        mni_fills = [];
        for i = 1:numel(properties)
            prop = properties{i};
            if ~isfield(kwave_medium, prop)
                continue
            end
            t1w_gz = fullfile(out_dir, sprintf('%s_T1w%s.nii.gz', prop, affix));
            mni_gz = fullfile(out_dir, sprintf('%s_MNI%s.nii.gz', prop, affix));
            if ~confirm_overwriting(mni_gz, parameters) || ~isfile(t1w_gz)
                continue
            end
            mni_in{end+1}    = t1w_gz;              %#ok<AGROW>
            mni_out{end+1}   = mni_gz;              %#ok<AGROW>
            mni_fills(end+1) = water_fill.(prop);   %#ok<AGROW>
        end
        if ~isempty(mni_in)
            mni_warp_method = 'simnibs';
            if isfield(parameters, 'analysis') && isfield(parameters.analysis, 'mni_warp_method')
                mni_warp_method = parameters.analysis.mni_warp_method;
            end
            convert_final_to_MNI_simnibs(mni_in, m2m_folder, mni_out, parameters, ...
                'interpolation_order', 0, 'FillValues', mni_fills, 'method', mni_warp_method);
        end
    end
catch ME
    prev = warning('off', 'backtrace');
    warn('nifti_medium_properties:writeError', ...
        'Could not write medium property NIfTIs: %s', ME.message);
    warning(prev);
end
end
