function acoustic_provenance = acoustic_water_baseline(parameters, acoustic_provenance, kgrid, source, sensor, source_labels)
% ACOUSTIC_WATER_BASELINE  Decide and run the free-water baseline simulation.
%
% Determines whether a water-baseline acoustic simulation is needed, runs it
% if so, and stores the result in acoustic_provenance.freefield_isppa_wcm2
% and the free-water beam-axis profile in acoustic_provenance.freefield_axial.
%
% Whether it runs is decided by NEEDS_WATER_BASELINE (modules.run_water_baseline
% = 1, 'auto' or 0). It is skipped when the calibration ISPPA already matches
% the target (common case: simulate at the drive level from the transducer
% library); the calibration value is then recorded as the baseline. The origin
% of a simulated baseline is recorded in acoustic_provenance.freefield_source.
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

if is_acoustic_dry_run(parameters)
    % No field will be simulated, so there is nothing to scale.
    fprintf('Skipping free-water baseline: acoustic dry run.\n');
    return;
end

[run_baseline, reason] = needs_water_baseline(parameters, acoustic_provenance);

cal_isppa = NaN;
if isfield(acoustic_provenance, 'tx_calibration') && ...
        isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2')
    cal_isppa = acoustic_provenance.tx_calibration.isppa_wcm2;
end
has_target = isfield(parameters, 'transducer') && ...
    any(arrayfun(@(t) isfield(t, 'target_isppa_wcm2') && ...
                      ~isempty(t.target_isppa_wcm2) && ...
                      any(isfinite(t.target_isppa_wcm2)), parameters.transducer));

if has_target && ~run_baseline && isfinite(cal_isppa)
    acoustic_provenance.freefield_isppa_wcm2 = cal_isppa;
    fprintf('Skipping free-water baseline: target ISPPA matches calibration level (%.2f W/cm²).\n', cal_isppa);
end

source_ready = (~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup == 1) && ...
    ~isempty(kgrid) && ~isempty(source);

if run_baseline && ~source_ready
    % BabelBrain runs its own water simulation; otherwise the cached
    % provenance supplies the baseline (ISPPA_SCALE_FACTOR errors if not).
    fprintf('Free-water baseline needed (%s) but no k-Wave source is set up in this run; ', reason);
    fprintf('the baseline must come from the cached acoustic result or the backend.\n');
elseif run_baseline
    if strcmp(reason, 'no_calibration')
        fprintf(['Running free-water baseline: transducer.target_isppa_wcm2 is set and the ' ...
                 'transducer has no calibration ISPPA (modules.run_water_baseline = ''auto'').\n']);
    end
    log_timer('start', 'freefield_baseline', parameters.io.dir_output);
    [acoustic_provenance.freefield_isppa_wcm2, acoustic_provenance.freefield_axial] = ...
        water_baseline(parameters, kgrid, source, sensor, source_labels);
    if isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'jwave')
        sim_source = 'jwave_water_sim';
    else
        sim_source = 'kwave_simulation';
    end
    acoustic_provenance.freefield_source = sim_source;
    if isfield(acoustic_provenance, 'tx_calibration')
        acoustic_provenance.tx_calibration.isppa_wcm2   = acoustic_provenance.freefield_isppa_wcm2;
        acoustic_provenance.tx_calibration.isppa_source = sim_source;
    end
    log_timer('stop', 'freefield_baseline');
end

end
