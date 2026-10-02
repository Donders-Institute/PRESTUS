function [run_baseline, reason] = needs_water_baseline(parameters, acoustic_provenance)
% NEEDS_WATER_BASELINE  Decide whether a free-water baseline simulation must run
%
% The free-water baseline (WATER_BASELINE) measures the ISPPA the current drive
% produces in water, so that APPLY_ISPPA_SCALING can scale the field to
% transducer.target_isppa_wcm2. modules.run_water_baseline selects when it runs:
%
%   1 / true  always
%   'auto'    (default) when a target is set and either no calibration ISPPA
%             is recorded for the transducer or the target differs from it by
%             more than 1 %
%   0 / false only when the target differs from a recorded calibration ISPPA
%             (an uncalibrated transducer with a target then fails in
%             ISPPA_SCALE_FACTOR instead of being reported unscaled)
%
% Use as:
%   [run_baseline, reason] = needs_water_baseline(parameters, acoustic_provenance)
%
% Input:
%   parameters          - PRESTUS parameter struct
%   acoustic_provenance - provenance struct; tx_calibration.isppa_wcm2 is the
%                         calibration ISPPA, when one is recorded
%
% Output:
%   run_baseline - logical
%   reason       - 'requested', 'no_calibration', 'target_differs', or '' (no run)
%
% See also: ACOUSTIC_WATER_BASELINE, ISPPA_SCALE_FACTOR, WATER_BASELINE

    reason = '';

    mode = 'auto';
    if isfield(parameters, 'modules') && isfield(parameters.modules, 'run_water_baseline') && ...
            ~isempty(parameters.modules.run_water_baseline)
        mode = parameters.modules.run_water_baseline;
    end
    if ischar(mode) || isstring(mode)
        mode = lower(char(mode));
        if ~strcmp(mode, 'auto')
            error('prestus:invalidRunWaterBaseline', ...
                'modules.run_water_baseline must be 0, 1 or ''auto'' (got ''%s'').', mode);
        end
    elseif isscalar(mode) && (islogical(mode) || isnumeric(mode))
        mode = double(mode) ~= 0;
    else
        error('prestus:invalidRunWaterBaseline', ...
            'modules.run_water_baseline must be 0, 1 or ''auto''.');
    end

    if isequal(mode, true)
        run_baseline = true;
        reason = 'requested';
        return
    end

    target = NaN;
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer) && ...
            isfield(parameters.transducer, 'target_isppa_wcm2')
        t = parameters.transducer(1).target_isppa_wcm2;
        if isnumeric(t) && ~isempty(t) && isfinite(t(1))
            target = double(t(1));
        end
    end

    cal_isppa = NaN;
    if isfield(acoustic_provenance, 'tx_calibration') && ...
            isfield(acoustic_provenance.tx_calibration, 'isppa_wcm2')
        c = acoustic_provenance.tx_calibration.isppa_wcm2;
        if isnumeric(c) && isscalar(c) && isfinite(c) && c > 0
            cal_isppa = double(c);
        end
    end

    run_baseline = false;
    if ~isfinite(target)
        return
    end
    if isfinite(cal_isppa)
        if abs(target - cal_isppa) / cal_isppa > 0.01
            run_baseline = true;
            reason = 'target_differs';
        end
    elseif strcmp(mode, 'auto')
        run_baseline = true;
        reason = 'no_calibration';
    end
end
