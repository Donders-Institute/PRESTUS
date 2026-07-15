function sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters)
% APPLY_ISPPA_SCALING  Scale the simulated pressure field so the free-water
%                      ISPPA matches parameters.transducer(1).target_isppa_wcm2.
%
% Rescaling is only applied when a target is set AND the baseline ISPPA came
% from an actually-run free-water simulation (isppa_source is one of the
% *_water_sim / *_simulation tags) and the two differ by more than 1%.  This
% covers the two legitimate cases:
%
%   * a multi-ISPPA sweep, whose Stage-1 always runs a free-water baseline and
%     scales each requested level from that measured baseline; and
%   * any run with modules.run_water_baseline=1 whose simulated free-water peak
%     does not correspond to the requested target.
%
% A single/coherent/sequential run is NOT force-rescaled off the stored library
% or hydrophone-measured calibration value — without a run free-water baseline
% the field is reported at its raw simulated drive level.  Set
% modules.run_water_baseline=1 to calibrate such a run to target_isppa_wcm2.
%
% The target is read from the first transducer entry (single-transducer and
% coherent multi-transducer runs share one combined field driven by
% transducer(1)).  For async multi-transducer runs each transducer is scaled
% independently inside COMBINE_ASYNC_INTENSITY before the fields are summed —
% this function is not called in that code path.
%
% All pressure fields relevant to downstream analysis are scaled by the same
% factor: p_max_all (Isppa/Psptp), p_min_all (peak rarefaction → MI/MI_tc), and
% p_final.  Scaling p_max_all alone would leave the mechanical index computed at
% the unscaled drive level, inconsistent with the reported intensities.
%
% Called after loading a cached acoustic result in a thermal-only job
% spawned by MULTI_ISPPA_PIPELINE or PRESTUS_PIPELINE.

    % Combined async/FOV caches have p_max_all (thermal) and p_max_async (safety)
    % already scaled to the per-transducer targets inside combine_async_intensity /
    % assemble_limited_fov_fields.  Re-applying scaling here would double-scale
    % p_max_all while leaving p_max_async untouched, producing an inconsistency
    % between thermal and acoustic safety fields.  Skip for those caches.
    if isfield(acoustic_provenance, 'combination_mode') && ...
            any(strcmp(acoustic_provenance.combination_mode, ...
                       {'async_dual_field', 'limited_fov_incoherent_sum'}))
        return;
    end

    target = get_target(parameters);
    if isempty(target) || ~isfinite(target)
        return;
    end

    [baseline, baseline_source] = get_baseline(acoustic_provenance);
    if isempty(baseline) || ~(baseline > 0)
        warning('apply_isppa_scaling:noBaseline', ...
            ['apply_isppa_scaling: neither freefield_isppa_wcm2 nor tx_calibration.isppa_wcm2 ' ...
             'found in provenance — pressure scaling skipped. Re-run the acoustic stage ' ...
             'or ensure the transducer library sets tx_calibration.isppa_wcm2.']);
        return;
    end

    % Only rescale off an actually-run free-water simulation.  A library or
    % hydrophone-measured calibration value is not a per-run baseline, so the
    % field is left at its raw simulated drive level (see header).
    SIM_SOURCES = {'kwave_simulation', 'jwave_water_sim', 'babelbrain_water_sim'};
    if ~any(strcmp(baseline_source, SIM_SOURCES))
        fprintf(['apply_isppa_scaling: baseline source ''%s'' is not a run free-water ' ...
                 'simulation — reporting raw simulated drive level (no target rescale). ' ...
                 'Set modules.run_water_baseline=1 to calibrate to target_isppa_wcm2.\n'], ...
                baseline_source);
        return;
    end

    scale_p  = sqrt(target / baseline);

    if abs(scale_p - 1) <= 0.01
        return;
    end

    if scale_p > 4 || scale_p < 0.25
        warning('apply_isppa_scaling:largeScale', ...
            ['apply_isppa_scaling: large pressure scale factor (%.2fx) from %.1f to ' ...
             '%.1f W/cm². Verify that the target is in the linear acoustic regime ' ...
             'relative to the baseline simulation amplitude.'], ...
            scale_p, baseline, target);
    end

    fprintf('Scaling cached pressure field: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
            baseline, target, scale_p);
    sensor_data.p_max_all = sensor_data.p_max_all * scale_p;
    % Scale peak-rarefaction and final fields with the same factor so that the
    % mechanical index (computed from p_min_all) stays consistent with the
    % scaled intensities (computed from p_max_all).
    if isfield(sensor_data, 'p_min_all')
        sensor_data.p_min_all = sensor_data.p_min_all * scale_p;
    end
    if isfield(sensor_data, 'p_final')
        sensor_data.p_final = sensor_data.p_final * scale_p;
    end
end

function [baseline, source] = get_baseline(acoustic_provenance)
% Resolve the free-water baseline ISPPA and the provenance source tag that
% records how it was obtained.  Prefers the free-water field if present, else
% the transducer-calibration value.  NOTE: freefield_isppa_wcm2 is populated
% with the calibration value in the baseline-skip path (no sim run), so the
% source tag — not the presence of freefield_isppa_wcm2 — determines whether an
% actual simulation produced it.
    baseline = [];
    source   = '';
    if isfield(acoustic_provenance, 'tx_calibration') && ...
            isfield(acoustic_provenance.tx_calibration, 'isppa_source')
        source = acoustic_provenance.tx_calibration.isppa_source;
    end
    if isfield(acoustic_provenance, 'freefield_isppa_wcm2') && ...
            ~isempty(acoustic_provenance.freefield_isppa_wcm2) && ...
            acoustic_provenance.freefield_isppa_wcm2 > 0
        baseline = acoustic_provenance.freefield_isppa_wcm2;
    elseif isfield(acoustic_provenance, 'tx_calibration') && ...
            isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2') && ...
            ~isempty(acoustic_provenance.tx_calibration.isppa_wcm2) && ...
            acoustic_provenance.tx_calibration.isppa_wcm2 > 0
        baseline = acoustic_provenance.tx_calibration.isppa_wcm2;
    end
end

function target = get_target(parameters)
% Read target_isppa_wcm2 from transducer(1), with fallback to the deprecated
% calibration field so that callers that bypass load_parameters still work.
    target = [];
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer) && ...
            isfield(parameters.transducer(1), 'target_isppa_wcm2')
        target = parameters.transducer(1).target_isppa_wcm2;
    elseif isfield(parameters, 'calibration') && ...
            isfield(parameters.calibration, 'target_isppa_wcm2')
        target = parameters.calibration.target_isppa_wcm2;
    end
    % Discard if it is a vector — scaling a vector target makes no sense here;
    % multi-ISPPA pipeline sets a scalar before calling this path.
    if numel(target) > 1
        target = [];
    end
end
