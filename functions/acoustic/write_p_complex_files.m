function sensor_data = write_p_complex_files(parameters, sensor_data, filename_sensor_data, meta)
% WRITE_P_COMPLEX_FILES  Save complex pressure fields to their own cache files
%
% With io.save_p_complex, writes each complex field of sensor_data
% (meta.field_names) to <name>_p_complex[_<f>Hz]<affix>.mat next to the
% acoustic results file, with p_complex_info (meta without field_names; the
% frequency of that field; its grid size), and removes the fields from
% sensor_data. Shared by every acoustic backend, so SUPERPOSE_COMPLEX_FIELDS
% reads the same file format whatever produced it.
%
% Use as:
%   sensor_data = write_p_complex_files(parameters, sensor_data, filename_sensor_data, meta)
%
% Input:
%   parameters          - PRESTUS config (io.save_p_complex)
%   sensor_data         - acoustic results holding the fields in meta.field_names
%   filename_sensor_data - path of the acoustic results file (..._results<affix>.mat)
%   meta                - struct: field_names (cellstr), freq_hz (one per field),
%                         phase_reference, plus any further p_complex_info fields
%
% Output:
%   sensor_data - without the complex fields
%
% See also: ACOUSTIC_WRAPPER, SUPERPOSE_COMPLEX_FIELDS

    save_it = isfield(parameters.io, 'save_p_complex') && parameters.io.save_p_complex;
    [fdir, fname, fext] = fileparts(filename_sensor_data);
    for fi = 1:numel(meta.field_names)
        fld = meta.field_names{fi};
        if ~isfield(sensor_data, fld)
            continue
        end
        if save_it
            % the single-frequency file keeps its original name
            freq_suffix = fld(length('p_complex')+1:end);   % '' or e.g. '_500000Hz'
            filename_p_complex = fullfile(fdir, [strrep(fname, '_results', ['_p_complex' freq_suffix]), fext]);
            p_complex = sensor_data.(fld);
            p_complex_info = rmfield(meta, 'field_names');
            p_complex_info.freq_hz   = meta.freq_hz(fi);
            p_complex_info.grid_dims = size(p_complex);
            save(filename_p_complex, 'p_complex', 'p_complex_info', '-v7.3');
            fprintf('Saved complex pressure field (%s) to %s\n', fld, filename_p_complex);
        end
        sensor_data = rmfield(sensor_data, fld);
    end
end
