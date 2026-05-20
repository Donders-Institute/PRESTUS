function nifti_thermal(parameters, planimg, results_heating, kwave_medium)
% NIFTI_THERMAL  Export thermal-stage data to NIfTI in T1 and MNI space
%
% Writes heating maps, heat-rise maps, and CEM43 maps back-transformed to T1
% space, with edge-artefact correction applied to the temperature map.
% Optionally converts all outputs to MNI space.
%
% The T1w back-transform (tformarray) is batched by fill-value group, and the
% MNI transform (subject2mni) is run in parallel across all volumes.
%
% Use as:
%   nifti_thermal(parameters, planimg, results_heating, kwave_medium)
%
% Input:
%   parameters      - PRESTUS config struct
%   planimg         - struct with fields: t1_image_orig, inv_transf, t1_header, transf
%   results_heating - struct with fields: maxT, heating_endT, CEM43, CEM43_end,
%                     CEM43_iso, CEM43_iso_end
%   kwave_medium    - struct with field: temp_0
%
% See also: NIFTI_MEDIUM, NIFTI_ACOUSTIC, CONVERT_FINAL_TO_MNI_SIMNIBS

arguments
    parameters      (1,1) struct
    planimg         (1,1) struct
    results_heating (1,1) struct
    kwave_medium    (1,1) struct
end

    if isfield(parameters.modules, 'run_nifti_creation') && ...
            parameters.modules.run_nifti_creation == 0
        disp('No nifti creation requested...')
        return
    end

    if ~contains(parameters.simulation.medium, {'layered'; 'phantom'})
        return
    end

    is_layered = strcmp(parameters.simulation.medium, 'layered');
    m2m_folder = fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));

    data_types = ["heating", "heating_end", "heatrise", "heatrise_end", ...
                  "CEM43", "CEM43_end", "CEM43_iso", "CEM43_iso_end"];
    n = numel(data_types);

    % Build file paths up front
    orig_files    = cell(1, n);
    orig_files_gz = cell(1, n);
    mni_files     = cell(1, n);
    for k = 1:n
        dt = data_types(k);
        orig_files{k}    = fullfile(parameters.io.dir_nii_T1w, ...
            sprintf('sub-%03d_%s_T1w%s_%s', parameters.subject_id, ...
                    parameters.simulation.medium, parameters.io.output_affix, dt));
        orig_files_gz{k} = strcat(orig_files{k}, '.nii.gz');
        mni_files{k}     = fullfile(parameters.io.dir_nii_MNI, ...
            sprintf('sub-%03d_%s_MNI%s_%s.nii.gz', parameters.subject_id, ...
                    parameters.simulation.medium, parameters.io.output_affix, dt));
    end

    % Check overwrite flags so we skip files that don't need writing
    needs_t1w = false(1, n);
    for k = 1:n
        needs_t1w(k) = logical(confirm_overwriting(orig_files_gz{k}, parameters));
    end

    % -------------------------------------------------------------------------
    % T1w back-transform — batch by fill-value group to call tformarray once
    % per group instead of once per volume.
    %
    % Group A (fill = water baseline): heating, heating_end  → indices 1–2
    % Group B (fill = 0):              everything else        → indices 3–8
    % -------------------------------------------------------------------------
    if is_layered
        orig_hdr          = planimg.t1_header;
        orig_hdr.Datatype = 'single';
        t1_size           = size(planimg.t1_image_orig);
        fill_water        = double(parameters.thermal.temp_0.water);

        % Group A
        grp_a = [1 2];
        if any(needs_t1w(grp_a))
            stack_a = cat(4, single(results_heating.maxT), single(results_heating.heating_endT));
            backtransf_a = tformarray(stack_a, planimg.inv_transf, ...
                makeresampler('linear', 'fill'), [1 2 3], [1 2 3], t1_size, [], fill_water);
        end

        % Group B — temperature rises use linear interpolation; CEM43 volumes use
        % nearest-neighbour to preserve spatial peaks across sequential-run adoption
        % (linear smoothing would reduce peak values causing a visible drop at the
        % run boundary when the file is re-imported by thermal_simulation.m).
        grp_b = 3:8;
        cem_local_idx = [3 4 5 6]; % indices within grp_b that are CEM43 types
        if any(needs_t1w(grp_b))
            stack_hrise = cat(4, ...
                single(results_heating.maxT - kwave_medium.temp_0), ...
                single(results_heating.heating_endT - kwave_medium.temp_0));
            backtransf_hrise = tformarray(stack_hrise, planimg.inv_transf, ...
                makeresampler('linear', 'fill'), [1 2 3], [1 2 3], t1_size, [], 0);

            stack_cem = cat(4, ...
                single(results_heating.CEM43), ...
                single(results_heating.CEM43_end), ...
                single(results_heating.CEM43_iso), ...
                single(results_heating.CEM43_iso_end));
            backtransf_cem = tformarray(stack_cem, planimg.inv_transf, ...
                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], t1_size, [], 0);
        end

        % Write Group A volumes
        for ki = 1:numel(grp_a)
            k = grp_a(ki);
            if ~needs_t1w(k); continue; end
            niftiwrite(backtransf_a(:,:,:,ki), orig_files{k}, orig_hdr, 'Compressed', true);
        end

        % Write Group B volumes (with CEM43 floor correction)
        for ki = 1:numel(grp_b)
            k = grp_b(ki);
            if ~needs_t1w(k); continue; end
            if ismember(ki, cem_local_idx)
                v = backtransf_cem(:,:,:,ki - 2);
                v(v <= 0) = 0.0000001;
            else
                v = backtransf_hrise(:,:,:,ki);
            end
            niftiwrite(v, orig_files{k}, orig_hdr, 'Compressed', true);
        end

    else
        % Phantom medium: no back-transform needed
        for k = 1:n
            if ~needs_t1w(k); continue; end
            switch data_types(k)
                case "heating",      data = single(results_heating.maxT);
                case "heating_end",  data = single(results_heating.heating_endT);
                case "heatrise",     data = single(results_heating.maxT - kwave_medium.temp_0);
                case "heatrise_end", data = single(results_heating.heating_endT - kwave_medium.temp_0);
                case "CEM43",        data = single(results_heating.CEM43);
                case "CEM43_end",    data = single(results_heating.CEM43_end);
                case "CEM43_iso",    data = single(results_heating.CEM43_iso);
                case "CEM43_iso_end",data = single(results_heating.CEM43_iso_end);
            end
            niftiwrite(data, orig_files{k}, 'Compressed', true);
        end
    end

    % -------------------------------------------------------------------------
    % MNI transform — batch all volumes into a single parallel subject2mni run
    % -------------------------------------------------------------------------
    if ~is_layered || ~should_save_output(parameters.io, 'save_MNI')
        return
    end

    % Collect volumes that need MNI conversion (T1w file must exist)
    mni_in    = {};
    mni_out   = {};
    mni_fills = [];
    for k = 1:n
        if ~confirm_overwriting(mni_files{k}, parameters); continue; end
        if ~isfile(orig_files_gz{k}); continue; end
        mni_in{end+1}    = orig_files_gz{k}; %#ok<AGROW>
        mni_out{end+1}   = mni_files{k};     %#ok<AGROW>
        % heating/heating_end: fill out-of-FOV with water baseline; others: 0
        if any(k == [1 2])
            mni_fills(end+1) = fill_water; %#ok<AGROW>
        else
            mni_fills(end+1) = 0;          %#ok<AGROW>
        end
    end

    if ~isempty(mni_in)
        convert_final_to_MNI_simnibs(mni_in, m2m_folder, mni_out, parameters, ...
            'interpolation_order', 0, 'FillValues', mni_fills);
    end
end
