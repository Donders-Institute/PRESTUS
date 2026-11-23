function sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters)
% APPLY_ISPPA_SCALING  Scale the simulated pressure field so its free-water ISPPA matches the target
%
% Rescaling is applied when a target (parameters.transducer(1).target_isppa_wcm2)
% is set and the free-water baseline (ISPPA_BASELINE) differs from it by more
% than 1 %. The baseline must come from an actually-run free-water simulation
% (k-Wave, jwave or BabelBrain):
%
%   * a multi-ISPPA sweep, whose Stage-1 always runs a free-water baseline and
%     scales each requested level from that measured baseline; and
%   * any run whose free-water baseline ran (modules.run_water_baseline = 'auto'
%     runs it whenever a target is set without a matching calibration).
%
% A run is never force-rescaled off the stored library or hydrophone-measured
% calibration value, and a requested target is never silently ignored: when
% the target needs a rescale but no simulated baseline exists, this errors
% (prestus:isppa:noBaseline / prestus:isppa:unsimulatedBaseline) rather than
% reporting the field at its raw simulated drive level. See ISPPA_SCALE_FACTOR.
%
% The target is read from the first transducer entry (single-transducer
% runs). Multi-transducer runs are not rescaled here: coherent runs set each
% transducer's drive to its own target before the simulation
% (COHERENT_DRIVE_CALIBRATION, recorded in acoustic_provenance.coherent_drive),
% and async runs scale each transducer inside COMBINE_ASYNC_INTENSITY.
%
% All pressure fields relevant to downstream analysis are scaled by the same
% factor: p_max_all (Isppa/Psptp), p_min_all (peak rarefaction → MI/MI_tc), and
% p_final.  Scaling p_max_all alone would leave the mechanical index computed at
% the unscaled drive level, inconsistent with the reported intensities.
%
% Called by ACOUSTIC_STAGE after a simulation and by CACHE_ACOUSTIC when a
% cached acoustic result is loaded.

    % Combined async/FOV caches have p_max_all (thermal) and p_max_async (safety)
    % already scaled to the per-transducer targets inside combine_async_intensity
    % (and, in older caches, the removed partial-FOV assembly, 'limited_fov_incoherent_sum').
    % Re-applying scaling here would double-scale
    % p_max_all while leaving p_max_async untouched, producing an inconsistency
    % between thermal and acoustic safety fields.  Skip for those caches.
    if isfield(acoustic_provenance, 'combination_mode') && ...
            (startsWith(acoustic_provenance.combination_mode, 'async_dual_field') || ...   % incl. _fov_backprojected
             strcmp(acoustic_provenance.combination_mode, 'limited_fov_incoherent_sum'))
        return;
    end

    % Coherent multi-transducer: drives were set per transducer before the
    % simulation; the combined field is final.
    if isfield(acoustic_provenance, 'coherent_drive') && ~isempty(acoustic_provenance.coherent_drive)
        return;
    end
    n_tx = 0;
    if isfield(parameters, 'transducer'), n_tx = numel(parameters.transducer); end
    if n_tx > 1 && ~is_async_mode(parameters) && isfinite_scalar(get_target(parameters))
        error('prestus:isppa:coherentWithoutDriveCalibration', ...
            ['This coherent %d-transducer result has no per-transducer drive calibration ' ...
             '(acoustic_provenance.coherent_drive): it was simulated before ISPPA targets were ' ...
             'applied per transducer, so it cannot be scaled to them. Re-run the acoustic ' ...
             'simulation (io.overwrite_files).'], n_tx);
    end

    target = get_target(parameters);
    [scale_p, baseline] = isppa_scale_factor(target, acoustic_provenance);
    if scale_p == 1
        return;
    end

    fprintf('Scaling cached pressure field: %.2f → %.2f W/cm² (scale_p = %.4f)\n', ...
            baseline, target, scale_p);
    sensor_data.p_max_all = sensor_data.p_max_all * scale_p;
    sensor_data.isppa_scale_p = scale_p;  % recorded for the free-water references in PLOT_AXIAL_PROFILE
    if isfield(sensor_data, 'p_min_all')
        sensor_data.p_min_all = sensor_data.p_min_all * scale_p;
    end
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

function tf = isfinite_scalar(x)
    tf = isnumeric(x) && isscalar(x) && isfinite(x);
end
