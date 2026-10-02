function [baseline, source] = isppa_baseline(acoustic_provenance)
% ISPPA_BASELINE  Free-water ISPPA an acoustic field was produced at, and its origin
%
% Reads the free-water reference ISPPA from acoustic provenance. A measured
% free-water baseline (freefield_isppa_wcm2, from ACOUSTIC_WATER_BASELINE or the
% BabelBrain water run) takes precedence over the library calibration value
% (tx_calibration.isppa_wcm2).
%
% The origin is read from freefield_source (written since the baseline fix),
% then tx_calibration.isppa_source. Caches written before freefield_source
% existed carry a freefield value without a calibration only when a free-water
% simulation produced it, so that case is reported as 'kwave_simulation'.
%
% Use as:
%   [baseline, source] = isppa_baseline(acoustic_provenance)
%
% Input:
%   acoustic_provenance - struct; may hold freefield_isppa_wcm2, freefield_source,
%                         tx_calibration.isppa_wcm2 / .isppa_source
%
% Output:
%   baseline - free-water ISPPA [W/cm²]; NaN when none is recorded
%   source   - 'kwave_simulation', 'jwave_water_sim', 'babelbrain_water_sim',
%              a library calibration source, or '' when unknown
%
% See also: ISPPA_SCALE_FACTOR, APPLY_ISPPA_SCALING, ACOUSTIC_WATER_BASELINE

    baseline = NaN;
    source   = '';
    if ~isstruct(acoustic_provenance)
        return
    end

    has_cal = isfield(acoustic_provenance, 'tx_calibration') && ...
        isstruct(acoustic_provenance.tx_calibration);
    cal_source = '';
    if has_cal && isfield(acoustic_provenance.tx_calibration, 'isppa_source')
        cal_source = char(acoustic_provenance.tx_calibration.isppa_source);
    end

    if positive_scalar(acoustic_provenance, 'freefield_isppa_wcm2')
        baseline = double(acoustic_provenance.freefield_isppa_wcm2);
        if isfield(acoustic_provenance, 'freefield_source') && ...
                ~isempty(acoustic_provenance.freefield_source)
            source = char(acoustic_provenance.freefield_source);
        elseif ~isempty(cal_source)
            source = cal_source;
        elseif ~has_cal
            source = 'kwave_simulation';   % legacy cache: only a simulation writes this
        end
    elseif has_cal && positive_scalar(acoustic_provenance.tx_calibration, 'isppa_wcm2')
        baseline = double(acoustic_provenance.tx_calibration.isppa_wcm2);
        source   = cal_source;
    end
end

function ok = positive_scalar(s, field)
    ok = isfield(s, field) && isnumeric(s.(field)) && isscalar(s.(field)) && ...
        isfinite(s.(field)) && s.(field) > 0;
end
