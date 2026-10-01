function acoustic_provenance = acoustic_water_baseline(parameters, acoustic_provenance, kgrid, source, sensor, source_labels)
% ACOUSTIC_WATER_BASELINE  Decide and run the free-water baseline simulation.
%
% Determines whether a water-baseline acoustic simulation is needed, runs it
% if so, and stores the result in acoustic_provenance.freefield_isppa_wcm2
% and the free-water beam-axis profile in acoustic_provenance.freefield_axial.
%
% A baseline is skipped when the calibration ISPPA already matches the target
% (common case: simulate at the drive level from the transducer library).
% It is forced when modules.run_water_baseline == 1 or when the target ISPPA
% differs from the calibration level by more than 1%.
%
% Use as:
%   acoustic_provenance = acoustic_water_baseline(parameters, acoustic_provenance, kgrid, source, sensor)
%   acoustic_provenance = acoustic_water_baseline(parameters, acoustic_provenance, kgrid, source, sensor, source_labels)

arguments
    parameters          (1,1) struct
    acoustic_provenance (1,1) struct
    kgrid
    source
    sensor
    source_labels = []
end

has_per_transducer_target = isfield(parameters, 'transducer') && ...
    any(arrayfun(@(t) isfield(t, 'target_isppa_wcm2') && ...
                      ~isempty(t.target_isppa_wcm2) && ...
                      any(isfinite(t.target_isppa_wcm2)), parameters.transducer));

cal_isppa = NaN;
if isfield(acoustic_provenance, 'tx_calibration') && ...
        isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2')
    cal_isppa = acoustic_provenance.tx_calibration.isppa_wcm2;
end

target_differs_from_cal = has_per_transducer_target && isfinite(cal_isppa) && ...
    abs(parameters.transducer(1).target_isppa_wcm2(1) - cal_isppa) / cal_isppa > 0.01;

run_baseline = (isfield(parameters.modules, 'run_water_baseline') && ...
                parameters.modules.run_water_baseline == 1) || ...
               target_differs_from_cal;

if has_per_transducer_target && ~run_baseline && isfinite(cal_isppa)
    acoustic_provenance.freefield_isppa_wcm2 = cal_isppa;
    fprintf('Skipping free-water baseline: target ISPPA matches calibration level (%.2f W/cm²).\n', cal_isppa);
end

if run_baseline && ...
        (~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup == 1)
    log_timer('start', 'freefield_baseline', parameters.io.dir_output);
    [acoustic_provenance.freefield_isppa_wcm2, acoustic_provenance.freefield_axial] = ...
        water_baseline(parameters, kgrid, source, sensor, source_labels);
    if isfield(acoustic_provenance, 'tx_calibration')
        acoustic_provenance.tx_calibration.isppa_wcm2   = acoustic_provenance.freefield_isppa_wcm2;
        if isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'jwave')
            acoustic_provenance.tx_calibration.isppa_source = 'jwave_water_sim';
        else
            acoustic_provenance.tx_calibration.isppa_source = 'kwave_simulation';
        end
    end
    log_timer('stop', 'freefield_baseline');
end

end
