function [sensor_data, acoustic_provenance] = load_acoustic_input(file, parameters, planimg, acoustic_provenance)
% LOAD_ACOUSTIC_INPUT  Start a run from a prior pressure map instead of simulating
%
% Reads a peak-pressure NIfTI in T1w space [Pa] (path.acoustic_input), e.g.
% the *_T1w*_pressure.nii.gz of an earlier PRESTUS run, and resamples it onto
% this run's acoustic grid as sensor_data.p_max_all, so acoustic analysis and
% the thermal stage run exactly as after a simulation.
%
% The sidecar next to the NIfTI (WRITE_ACOUSTIC_SIDECAR) is used when present:
%   - its signature (CACHE_SIGNATURE 'acoustic': subject, medium, grid,
%     transducer) is compared with this run's; a mismatch is handled per
%     io.cache_signature_check (error | warn | off);
%   - its isppa.level_wcm2 (free-water ISPPA the map corresponds to) lets the
%     map be rescaled to this run's transducer.target_isppa_wcm2.
% Without a sidecar (e.g. a BabelBrain or other external map) the map is used
% as it is, with a warning.
%
% The map holds no peak negative pressure: the mechanical index falls back to
% the peak pressure (ACOUSTIC_ANALYSIS warns). The T1w map is at the T1's
% resolution, so a map from a finer simulation grid has lost some peak
% sharpness in the export.
%
% Use as:
%   [sensor_data, acoustic_provenance] = ...
%       load_acoustic_input(file, parameters, planimg, acoustic_provenance)
%
% Input:
%   file                - pressure NIfTI [Pa] in the T1w space of this subject
%   parameters          - PRESTUS parameters at the acoustic stage
%   planimg             - planning image struct (transf: T1 → grid; t1_image_orig;
%                         fov_offset_ac)
%   acoustic_provenance - provenance (full_ac_dims when an acoustic FOV is active;
%                         input_signature)
%
% Output:
%   sensor_data         - struct with p_max_all [Pa] on the acoustic grid
%                         (and isppa_scale_p relative to the simulated drive)
%   acoustic_provenance - with freefield baseline and acoustic_input recorded
%
% See also: WRITE_ACOUSTIC_SIDECAR, ACOUSTIC_INPUT_FILE, ACOUSTIC_STAGE

arguments
    file                {mustBeTextScalar}
    parameters          (1,1) struct
    planimg             (1,1) struct
    acoustic_provenance (1,1) struct
end

    file = char(file);
    if ~isfile(file)
        error('prestus:acousticInput:notFound', 'path.acoustic_input not found: %s', file);
    end
    if ~isfield(planimg, 'transf') || isempty(planimg.transf) || ...
            ~isfield(planimg, 't1_image_orig') || isempty(planimg.t1_image_orig)
        error('prestus:acousticInput:noT1Space', ...
            ['path.acoustic_input needs the T1 → grid transform of the grid stage ' ...
             '(simulation.medium layered or phantom with modules.run_grid_setup: 1).']);
    end
    if numel(parameters.grid.dims) ~= 3 || ...
            (isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric == 1)
        error('prestus:acousticInput:not3D', 'path.acoustic_input needs a 3D (not axisymmetric) grid.');
    end

    fprintf('[pipeline] Starting from the pressure map in path.acoustic_input:\n  %s\n', file);
    vol = single(niftiread(file));
    t1_dims = size(planimg.t1_image_orig);
    if ~isequal(size(vol), t1_dims)
        error('prestus:acousticInput:spaceMismatch', ...
            ['path.acoustic_input has dimensions [%s], the T1 of this run [%s]: it must be ' ...
             'a map in this subject''s T1w space (e.g. a *_T1w*_pressure.nii.gz).'], ...
            num2str(size(vol)), num2str(t1_dims));
    end

    % --- Sidecar: settings check and the ISPPA level of the map ---------------
    level = [];
    stored_scale = 1;
    sidecar = find_sidecar(file);
    remedy = ['Point path.acoustic_input at a map simulated with the current settings, or ' ...
              'clear it to simulate'];
    current_sig = struct();
    if isfield(acoustic_provenance, 'input_signature')
        current_sig = acoustic_provenance.input_signature;
    end
    if isempty(sidecar)
        check_sidecar_signature(file, [], current_sig, parameters, remedy);   % warns: unverified
    else
        sc = read_sidecar(sidecar);
        if ~isfield(sc, 'kind') || ~strcmp(char(sc.kind), 'acoustic_pressure')
            warning('prestus:acousticInput:kind', ...
                'Sidecar %s is not an acoustic_pressure sidecar; its settings are not checked.', sidecar);
        else
            stored_sig = [];
            if isfield(sc, 'signature'); stored_sig = sc.signature; end
            check_sidecar_signature(sidecar, stored_sig, current_sig, parameters, remedy);
            if isfield(sc, 'isppa') && isstruct(sc.isppa)
                level = positive_or_empty(sc.isppa, 'level_wcm2');
                s = positive_or_empty(sc.isppa, 'scale_p');
                if ~isempty(s); stored_scale = s; end
                b = positive_or_empty(sc.isppa, 'freefield_baseline_wcm2');
                if ~isempty(b)
                    acoustic_provenance.freefield_isppa_wcm2 = b;
                    if isfield(sc.isppa, 'baseline_source') && ischar(sc.isppa.baseline_source)
                        acoustic_provenance.freefield_source = sc.isppa.baseline_source;
                    end
                end
            end
        end
    end

    % --- Rescale to this run's target ISPPA -----------------------------------
    target = [];
    t = parameters.transducer(1);
    if numel(parameters.transducer) > 1 && any(arrayfun(@(x) isfield(x, 'target_isppa_wcm2') && ...
            isnumeric(x.target_isppa_wcm2) && any(isfinite(x.target_isppa_wcm2)), parameters.transducer))
        error('prestus:isppa:multiTransducerAcousticInput', ...
            ['path.acoustic_input with several transducers: a combined pressure map cannot be ' ...
             'rescaled to per-transducer ISPPA targets. Set target_isppa_wcm2: [] for every ' ...
             'transducer to use the map as it is.']);
    end
    if isfield(t, 'target_isppa_wcm2') && isnumeric(t.target_isppa_wcm2) && ...
            isscalar(t.target_isppa_wcm2) && isfinite(t.target_isppa_wcm2)
        target = double(t.target_isppa_wcm2);
    end
    scale_p = 1;
    if ~isempty(target) && ~isempty(level)
        scale_p = isppa_scale_factor(target, ...
            struct('freefield_isppa_wcm2', level, 'freefield_source', 'acoustic_input'), ...
            'AllowCalibrationBaseline', true, 'Context', 'path.acoustic_input');
        if scale_p ~= 1
            fprintf('Scaling the input pressure map: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
                level, target, scale_p);
        end
    elseif ~isempty(target)
        warning('prestus:acousticInput:unscaled', ...
            ['The free-water ISPPA of %s is unknown (no PRESTUS sidecar), so it cannot be ' ...
             'scaled to transducer.target_isppa_wcm2 = %.2f W/cm²; it is used as it is.'], file, target);
    end

    % --- Resample T1w → acoustic grid (full grid, then the acoustic FOV) ------
    full_dims = parameters.grid.dims;
    if isfield(acoustic_provenance, 'full_ac_dims') && ~isempty(acoustic_provenance.full_ac_dims)
        full_dims = acoustic_provenance.full_ac_dims;
    end
    p = max(single(affine_resample_3d(vol, planimg.transf, full_dims, 'cubic', 0)), 0);
    clear vol
    if ~isequal(full_dims(:)', parameters.grid.dims(:)')
        o = planimg.fov_offset_ac;
        d = parameters.grid.dims;
        p = p(o(1):o(1)+d(1)-1, o(2):o(2)+d(2)-1, o(3):o(3)+d(3)-1);
    end

    sensor_data = struct('p_max_all', p * scale_p);
    if stored_scale * scale_p ~= 1
        sensor_data.isppa_scale_p = stored_scale * scale_p;   % relative to the simulated drive
    end
    acoustic_provenance.acoustic_input = file;
end

function v = positive_or_empty(s, field)
    v = [];
    if isfield(s, field) && isnumeric(s.(field)) && isscalar(s.(field)) && ...
            isfinite(s.(field)) && s.(field) > 0
        v = double(s.(field));
    end
end
