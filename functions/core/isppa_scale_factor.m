function [scale_p, baseline, source] = isppa_scale_factor(target, acoustic_provenance, options)
% ISPPA_SCALE_FACTOR  Pressure scale factor that brings a field to a target free-water ISPPA
%
% Shared by APPLY_ISPPA_SCALING, COMBINE_ASYNC_INTENSITY and
% ASSEMBLE_LIMITED_FOV_FIELDS. Pressure scales linearly with drive, so
% scale_p = sqrt(target / baseline).
%
% A requested target is never silently ignored: when a finite target needs a
% rescale (more than 1 % in pressure) but no usable baseline exists, this errors
% instead of returning the field at the simulated drive level.
%
% Use as:
%   scale_p = isppa_scale_factor(target, acoustic_provenance)
%   [scale_p, baseline, source] = isppa_scale_factor(..., Context='transducer 2', ...)
%
% Input:
%   target              - target free-water ISPPA [W/cm²]; NaN/[] = no target (scale 1)
%   acoustic_provenance - provenance struct (see ISPPA_BASELINE)
%   options.AllowCalibrationBaseline - accept a library calibration value as the
%                         baseline (default false: only a free-water simulation)
%   options.Context     - label for messages, e.g. 'transducer 2' (default '')
%
% Output:
%   scale_p  - pressure scale factor (1 when no rescale is needed)
%   baseline - baseline ISPPA [W/cm²] used (NaN if none)
%   source   - origin of the baseline (see ISPPA_BASELINE)
%
% See also: ISPPA_BASELINE, APPLY_ISPPA_SCALING, ACOUSTIC_WATER_BASELINE

    arguments
        target
        acoustic_provenance
        options.AllowCalibrationBaseline (1,1) logical = false
        options.Context char = ''
    end

    SIM_SOURCES = {'kwave_simulation', 'jwave_water_sim', 'babelbrain_water_sim'};

    scale_p = 1;
    [baseline, source] = isppa_baseline(acoustic_provenance);

    if isempty(target) || ~isscalar(target) || ~isfinite(target)
        return
    end

    ctx = '';
    if ~isempty(options.Context)
        ctx = sprintf(' (%s)', options.Context);
    end

    if ~(baseline > 0)
        error('prestus:isppa:noBaseline', ...
            ['target_isppa_wcm2 = %.2f W/cm²%s, but no free-water baseline is recorded, so the ' ...
             'field cannot be scaled to it. Set modules.run_water_baseline to ''auto'' or 1 ' ...
             '(with modules.run_source_setup = 1) and re-run the acoustic stage, or remove ' ...
             'target_isppa_wcm2 to report the simulated drive level.'], target, ctx);
    end

    s = sqrt(target / baseline);
    if abs(s - 1) <= 0.01
        return
    end

    if ~options.AllowCalibrationBaseline && ~any(strcmp(source, SIM_SOURCES))
        error('prestus:isppa:unsimulatedBaseline', ...
            ['target_isppa_wcm2 = %.2f W/cm²%s differs from the recorded baseline %.2f W/cm² ' ...
             '(source ''%s''), which is not a free-water simulation. Set ' ...
             'modules.run_water_baseline to ''auto'' or 1 (with modules.run_source_setup = 1) ' ...
             'and re-run the acoustic stage.'], target, ctx, baseline, source);
    end

    if s > 4 || s < 0.25
        warning('prestus:isppa:largeScale', ...
            ['Large pressure scale factor%s (%.2fx) from %.1f to %.1f W/cm². Verify that ' ...
             'the target is in the linear acoustic regime relative to the baseline amplitude.'], ...
            ctx, s, baseline, target);
    end
    scale_p = s;
end
