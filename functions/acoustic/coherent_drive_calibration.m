function [parameters, drive_cal] = coherent_drive_calibration(parameters, simulate_baseline)

% COHERENT_DRIVE_CALIBRATION  Set each transducer's drive to its own ISPPA target (coherent runs)
%
% In a coherent multi-transducer run all transducers fire together in one
% simulation, so the combined field cannot be scaled per transducer
% afterwards. Instead, each transducer's drive is set before the simulation
% so that the transducer alone gives its own target free-water ISPPA:
%
%   elem_p_amp(k) <- elem_p_amp(k) * sqrt(target(k) / baseline(k))
%
% where baseline(k) is the free-water ISPPA of transducer k alone at its
% configured drive. As for a single transducer (NEEDS_WATER_BASELINE), the
% baseline is the transducer's calibration ISPPA (tx_calibration.isppa_wcm2,
% e.g. from the equipment library) when the target matches it within 1 %
% (no rescale), and otherwise a free-water simulation of transducer k alone
% (the other transducers' amplitudes set to zero, same grid and source
% geometry). The combined field is then simulated at these drives and is not
% rescaled (APPLY_ISPPA_SCALING skips runs carrying this calibration).
%
% A transducer without a finite target keeps its configured drive. Vector
% targets (sweeps) are not possible for coherent multi-transducer runs: each
% sweep point needs its own simulation (MULTI_ISPPA_PIPELINE errors).
%
% Use as:
%   [parameters, drive_cal] = coherent_drive_calibration(parameters, simulate_baseline)
%
% Input:
%   parameters        - PRESTUS config, coherent coupling, >= 2 transducers
%   simulate_baseline - function handle @(p) -> free-water ISPPA [W/cm²] of the
%                       source described by p (the pipeline builds the source
%                       and runs WATER_BASELINE); called with all transducers
%                       but k at zero amplitude
%
% Output:
%   parameters - elem_p_amp of each transducer with a target rescaled
%   drive_cal  - [1 x N] struct: target_isppa_wcm2, baseline_isppa_wcm2,
%                baseline_source ('calibration' | 'simulation' | ''),
%                amplitude_scale, reason (as NEEDS_WATER_BASELINE)
%
% See also: APPLY_ISPPA_SCALING, NEEDS_WATER_BASELINE, WATER_BASELINE, ISPPA_SCALE_FACTOR

    arguments
        parameters        (1,1) struct
        simulate_baseline (1,1) function_handle
    end

    n_tx = numel(parameters.transducer);
    drive_cal = repmat(struct('target_isppa_wcm2', NaN, 'baseline_isppa_wcm2', NaN, ...
        'baseline_source', '', 'amplitude_scale', 1, 'reason', ''), 1, n_tx);

    for k = 1:n_tx
        tr = parameters.transducer(k);
        target = NaN;
        if isfield(tr, 'target_isppa_wcm2') && isnumeric(tr.target_isppa_wcm2) && ~isempty(tr.target_isppa_wcm2)
            if numel(tr.target_isppa_wcm2) > 1
                error('prestus:isppa:coherentSweep', ...
                    ['Transducer %d has %d ISPPA targets. A coherent multi-transducer run takes one ' ...
                     'target per transducer (each sweep point needs its own simulation): run one ' ...
                     'config per point, or use simulation.transducer_coupling: async.'], ...
                    k, numel(tr.target_isppa_wcm2));
            end
            target = double(tr.target_isppa_wcm2);
        end
        drive_cal(k).target_isppa_wcm2 = target;
        if ~isfinite(target)
            continue   % no target: configured drive
        end

        % Same decision as for one transducer: its own target vs its own calibration
        p_k = parameters;
        p_k.transducer = tr;
        prov_k = struct();
        if isfield(tr, 'tx_calibration'), prov_k.tx_calibration = tr.tx_calibration; end
        [run_sim, reason] = needs_water_baseline(p_k, prov_k);
        drive_cal(k).reason = reason;

        if run_sim
            p_sim = parameters;
            for j = setdiff(1:n_tx, k)
                g = p_sim.transducer(j).type;
                p_sim.transducer(j).(g).elem_p_amp = 0 * p_sim.transducer(j).(g).elem_p_amp;
            end
            fprintf('[coherent_drive_calibration] Transducer %d: free-water baseline of transducer %d alone (%s).\n', ...
                k, k, reason);
            baseline = simulate_baseline(p_sim);
            drive_cal(k).baseline_source = 'simulation';
        else
            % Target matches the calibration (or none is needed): calibration ISPPA
            baseline = NaN;
            if isfield(prov_k, 'tx_calibration') && isfield(prov_k.tx_calibration, 'isppa_wcm2')
                baseline = double(prov_k.tx_calibration.isppa_wcm2);
            end
            if ~(isscalar(baseline) && isfinite(baseline) && baseline > 0)
                error('prestus:isppa:noBaseline', ...
                    ['Transducer %d has target_isppa_wcm2 = %g W/cm² but no calibration ISPPA, and ' ...
                     'modules.run_water_baseline = 0 prevents a free-water baseline. Set it to ''auto''.'], ...
                    k, target);
            end
            drive_cal(k).baseline_source = 'calibration';
        end
        if ~(isscalar(baseline) && isfinite(baseline) && baseline > 0)
            error('prestus:isppa:badBaseline', ...
                'Transducer %d: free-water baseline %g W/cm² is not usable.', k, baseline);
        end
        drive_cal(k).baseline_isppa_wcm2 = baseline;

        s = sqrt(target / baseline);
        if abs(target - baseline) / baseline <= 0.01
            s = 1;   % within 1 %: as APPLY_ISPPA_SCALING, no rescale
        end
        drive_cal(k).amplitude_scale = s;
        g = tr.type;
        parameters.transducer(k).(g).elem_p_amp = s * tr.(g).elem_p_amp;
        fprintf(['[coherent_drive_calibration] Transducer %d: target %.2f W/cm², baseline %.2f W/cm² (%s) ' ...
                 '-> drive x %.4f.\n'], k, target, baseline, drive_cal(k).baseline_source, s);
    end
end
