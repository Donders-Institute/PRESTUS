function sidecar_file = write_acoustic_sidecar(parameters, acoustic_provenance, sensor_data)
% WRITE_ACOUSTIC_SIDECAR  Describe the T1w pressure NIfTI so a later run can start from it
%
% Writes <...>_T1w<affix>_pressure.json next to the pressure map NIFTI_ACOUSTIC
% wrote. A later run points path.acoustic_input at the NIfTI to run acoustic
% analysis and thermal simulation without simulating (LOAD_ACOUSTIC_INPUT);
% the sidecar tells it the free-water ISPPA the map corresponds to (so it can
% be rescaled to another target) and the settings it was simulated with
% (compared on reuse, CACHE_SIGNATURE 'acoustic').
%
% Nothing is written when the pressure NIfTI does not exist (NIfTI creation
% off, or a medium without T1 space). Never errors.
%
% Use as:
%   sidecar_file = write_acoustic_sidecar(parameters, acoustic_provenance, sensor_data)
%
% Input:
%   parameters          - PRESTUS parameters after the acoustic stage
%   acoustic_provenance - provenance from the acoustic stage (baseline, signature)
%   sensor_data         - acoustic results (isppa_scale_p when scaled)
%
% Output:
%   sidecar_file - path written, or ''
%
% See also: LOAD_ACOUSTIC_INPUT, NIFTI_ACOUSTIC, WRITE_SIDECAR

    sidecar_file = '';
    try
        if ~isfield(parameters.io, 'dir_nii_T1w')
            return
        end
        nii = fullfile(parameters.io.dir_nii_T1w, sprintf('sub-%03d_%s_T1w%s_pressure.nii.gz', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        if ~isfile(nii)
            return
        end
        [~, json_file] = find_sidecar(nii);

        [baseline, baseline_source] = isppa_baseline(acoustic_provenance);
        target = [];
        t = parameters.transducer(1);
        if isfield(t, 'target_isppa_wcm2') && isnumeric(t.target_isppa_wcm2) && ...
                isscalar(t.target_isppa_wcm2) && isfinite(t.target_isppa_wcm2)
            target = double(t.target_isppa_wcm2);
        end
        scale_p = 1;
        if isstruct(sensor_data) && isfield(sensor_data, 'isppa_scale_p')
            scale_p = double(sensor_data.isppa_scale_p);
        end
        combined = isfield(acoustic_provenance, 'combination_mode') && ...
            ~isempty(acoustic_provenance.combination_mode);
        % Free-water ISPPA the stored field corresponds to: the target it was
        % scaled to, else the free-water ISPPA of the simulated drive.
        level = [];
        if ~combined
            if ~isempty(target) && baseline > 0
                level = target;
            elseif isempty(target) && baseline > 0
                level = baseline;
            end
        end
        isppa = struct( ...
            'level_wcm2',              level, ...
            'target_wcm2',             target, ...
            'freefield_baseline_wcm2', nan_to_empty(baseline), ...
            'baseline_source',         baseline_source, ...
            'scale_p',                 scale_p);
        if combined
            isppa.combination_mode = char(acoustic_provenance.combination_mode);
        end

        tr = cell(1, numel(parameters.transducer));
        for ti = 1:numel(parameters.transducer)
            ti_s = parameters.transducer(ti);
            e = struct('index', ti, 'serial', '', 'freq_hz', []);
            if isfield(ti_s, 'serial') && ~isempty(ti_s.serial); e.serial = char(ti_s.serial); end
            if isfield(ti_s, 'freq_hz'); e.freq_hz = double(ti_s.freq_hz); end
            if isfield(parameters, 'state') && isfield(parameters.state, 'planned_trans_t1_vox') && ...
                    size(parameters.state.planned_trans_t1_vox, 1) >= ti
                e.trans_pos_t1_vox = double(parameters.state.planned_trans_t1_vox(ti, :));
                e.focus_pos_t1_vox = double(parameters.state.planned_focus_t1_vox(ti, :));
            end
            tr{ti} = e;
        end

        backend = 'kwave';
        if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
            backend = char(parameters.acoustic_backend);
        end
        [~, nii_name, nii_ext] = fileparts(nii);
        content = struct( ...
            'files',       struct('pressure', [nii_name nii_ext]), ...
            'quantity',    'peak pressure amplitude (temporal maximum, p_max_all)', ...
            'units',       'Pa', ...
            'space',       'T1w', ...
            'isppa',       isppa, ...
            'transducers', {tr}, ...
            'simulation',  struct('medium', char(parameters.simulation.medium), ...
                                  'acoustic_backend', backend, ...
                                  'grid_resolution_mm', parameters.grid.resolution_mm), ...
            'source_run',  struct('dir_output', char(parameters.io.dir_output), ...
                                  'output_affix', char(parameters.io.output_affix)), ...
            'signature',   []);
        if isfield(acoustic_provenance, 'input_signature')
            content.signature = acoustic_provenance.input_signature;
        end
        sidecar_file = write_sidecar(json_file, 'acoustic_pressure', content, parameters);
    catch ME
        warning('prestus:sidecar:notWritten', ...
            'Could not write the acoustic pressure sidecar: %s', ME.message);
        sidecar_file = '';
    end
end

function v = nan_to_empty(v)
    if ~isfinite(v); v = []; end
end
