function sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters)
% Scale sensor_data.p_max_all (and p_final if present) so that the
% free-water ISPPA matches parameters.calibration.target_isppa_wcm2.
% Only applies when both the target and the baseline ISPPA are available
% and differ by more than 1%.
%
% Called after loading a cached acoustic result inside a thermal-only job
% spawned by multi_isppa_pipeline. The scale factor is derived entirely
% from the user-facing ISPPA values — elem_amp is never inspected.
    if ~isfield(parameters, 'calibration') || ...
            ~isfield(parameters.calibration, 'target_isppa_wcm2') || ...
            isempty(parameters.calibration.target_isppa_wcm2) || ...
            ~isfinite(parameters.calibration.target_isppa_wcm2)
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

    target   = parameters.calibration.target_isppa_wcm2;
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
