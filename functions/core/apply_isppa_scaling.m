function sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters)
% APPLY_ISPPA_SCALING  Scale sensor_data.p_max_all so the free-water ISPPA
%                      matches parameters.transducer(1).target_isppa_wcm2.
%
% Only applies when both the target and the free-water baseline ISPPA are
% available and differ by more than 1%.  The target is read from the first
% transducer entry (single-transducer and coherent multi-transducer runs
% share one combined acoustic field driven by transducer(1) settings).
%
% For async multi-transducer runs each transducer is scaled independently
% inside COMBINE_ASYNC_INTENSITY before the fields are summed — this
% function is not called in that code path.
%
% Called after loading a cached acoustic result in a thermal-only job
% spawned by MULTI_ISPPA_PIPELINE or PRESTUS_PIPELINE.

    target = get_target(parameters);
    if isempty(target) || ~isfinite(target)
        return;
    end

    if ~isfield(acoustic_provenance, 'freefield_isppa_wcm2') || ...
            isempty(acoustic_provenance.freefield_isppa_wcm2) || ...
            acoustic_provenance.freefield_isppa_wcm2 <= 0
        warn(['apply_isppa_scaling: acoustic_provenance.freefield_isppa_wcm2 is missing ' ...
                 'or invalid — pressure scaling skipped. Re-run the acoustic stage to ' ...
                 'generate provenance.']);
        return;
    end

    baseline = acoustic_provenance.freefield_isppa_wcm2;
    scale_p  = sqrt(target / baseline);

    if abs(scale_p - 1) <= 0.01
        return;
    end

    if scale_p > 4 || scale_p < 0.25
        warn(['apply_isppa_scaling: large pressure scale factor (%.2fx) from %.1f to ' ...
                 '%.1f W/cm². Verify that the target is in the linear acoustic regime ' ...
                 'relative to the baseline simulation amplitude.'], ...
                scale_p, baseline, target);
    end

    fprintf('Scaling cached pressure field: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
            baseline, target, scale_p);
    sensor_data.p_max_all = sensor_data.p_max_all * scale_p;
    if isfield(sensor_data, 'p_final')
        sensor_data.p_final = sensor_data.p_final * scale_p;
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
