function [msg, is_mismatch] = isppa_intent(parameters, t_i)
% ISPPA_INTENT  One-line statement of the free-water ISPPA a run will model
%
% Compares transducer(t_i).target_isppa_wcm2 with the ISPPA the configured
% drive was calibrated for (transducer.tx_calibration.isppa_wcm2, written by
% the calibration pipeline and the transducer library) and says which of the
% two the results will report:
%
%   * no target                → the drive level as configured (no scaling)
%   * target = calibration     → the calibrated intensity (no scaling)
%   * target ≠ calibration     → rescaled to the target (mismatch)
%   * target, no calibration   → rescaled to the target from a free-water
%                                baseline (mismatch when the target is the
%                                toolbox default rather than set by the user)
%
% A target counts as the toolbox default when no config layer or override set
% it and it still holds the default value (parameters.state.target_isppa_default_wcm2,
% recorded by LOAD_PARAMETERS).
%
% Use as:
%   [msg, is_mismatch] = isppa_intent(parameters, t_i)
%
% Input:
%   parameters - PRESTUS config
%   t_i        - transducer index (default 1)
%
% Output:
%   msg         - char, the statement (without trailing newline)
%   is_mismatch - logical, true when the modelled ISPPA differs from the
%                 calibrated one, or follows the default without a calibration
%
% See also: NEEDS_WATER_BASELINE, APPLY_ISPPA_SCALING, PRINT_PARAMETER_SUMMARY

    if nargin < 2; t_i = 1; end
    msg = '';
    is_mismatch = false;
    if ~isfield(parameters, 'transducer') || numel(parameters.transducer) < t_i
        return
    end
    tr = parameters.transducer(t_i);

    target = [];
    if isfield(tr, 'target_isppa_wcm2') && isnumeric(tr.target_isppa_wcm2)
        target = double(tr.target_isppa_wcm2(:)');
        target = target(isfinite(target));
    end
    cal = NaN;
    if isfield(tr, 'tx_calibration') && isstruct(tr.tx_calibration) && ...
            isfield(tr.tx_calibration, 'isppa_wcm2')
        c = tr.tx_calibration.isppa_wcm2;
        if isnumeric(c) && isscalar(c) && isfinite(c) && c > 0
            cal = double(c);
        end
    end
    from_default = false;
    if ~isempty(target) && isfield(parameters, 'state') && ...
            isfield(parameters.state, 'target_isppa_default_wcm2')
        from_default = isequal(target, double(parameters.state.target_isppa_default_wcm2(:)'));
    end
    src = '';
    if from_default; src = ' (toolbox default; not set in your config)'; end

    if isempty(target)
        msg = 'drive level as configured, no scaling (target_isppa_wcm2: [])';
        if isfinite(cal)
            msg = sprintf('%s; drive calibrated for %g W/cm²', msg, cal);
        end
    elseif numel(target) > 1
        msg = sprintf('sweep of %s W/cm² free-water ISPPA%s', mat2str(target), src);
        if isfinite(cal)
            msg = sprintf('%s; drive calibrated for %g W/cm²', msg, cal);
        end
    elseif isfinite(cal) && abs(target - cal) / cal <= 0.01
        msg = sprintf('%g W/cm² free-water ISPPA, the calibrated intensity of the drive (no scaling)', target);
    elseif isfinite(cal)
        is_mismatch = true;
        msg = sprintf(['%g W/cm² free-water ISPPA%s, but the drive is calibrated for %g W/cm²: ' ...
            'results are rescaled to %g W/cm². Set transducer.target_isppa_wcm2: %g to model ' ...
            'the calibrated intensity'], target, src, cal, target, cal);
    else
        is_mismatch = from_default;
        msg = sprintf(['%g W/cm² free-water ISPPA%s, scaled from a free-water baseline ' ...
            '(the drive carries no calibration ISPPA, transducer.tx_calibration.isppa_wcm2)'], target, src);
    end
end
