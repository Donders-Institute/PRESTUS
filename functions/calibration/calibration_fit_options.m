function fo = calibration_fit_options(preset, overrides)
% CALIBRATION_FIT_OPTIONS  FIT_TRANSDUCER_MODEL options for a per-depth calibration preset
%
% Single source of the per-depth fit settings shared by BUILD_TRANSDUCER_LIBRARY,
% REFINE_TRANSDUCER_LIBRARY, CALIBRATE_LIBRARY_ENTRY and CALIBRATION_TRANSDUCER.
%
% Presets:
%   'rayleigh_phase_gain' (default) - per-ring phase offset and gain, Rayleigh forward model,
%               intensity-domain pseudo-Huber loss (delta 0.8), weak pull to the
%               prior (Lambda 1e-3), near field below 10 mm ignored. The axial
%               shift is not fitted: it is not applied at run time.
%   'oneil_phase_only' - the behaviour of the former PERFORM_GLOBAL_SEARCH path: O'Neil
%               forward model, FEX-minimize multi-start search, phases only, no
%               pull to the prior, (practically) least-squares loss on
%               intensity, no near-field penalty. Gaussian focal weighting
%               (calibration.opt_weights) maps to FocalWeight.
%
% Use as:
%   fo = calibration_fit_options()                    % 'rayleigh_phase_gain'
%   fo = calibration_fit_options(preset)
%   fo = calibration_fit_options(preset, overrides)   % overrides win
%
% Input:
%   preset    - 'rayleigh_phase_gain' (default) | 'oneil_phase_only'
%   overrides - struct of FIT_TRANSDUCER_MODEL options applied on top
%
% Output:
%   fo - options struct for FIT_TRANSDUCER_MODEL
%
% See also: FIT_TRANSDUCER_MODEL, CALIBRATION_FIT_DEPTHS, CALIBRATE_LIBRARY_ENTRY

    if nargin < 1 || isempty(preset); preset = 'rayleigh_phase_gain'; end
    if nargin < 2 || isempty(overrides); overrides = struct(); end

    switch char(preset)
        case 'rayleigh_phase_gain'
            fo = struct('FreeParams', {{'phase', 'amp'}}, 'ForwardModel', 'rayleigh', ...
                'Domain', 'intensity', 'HuberDelta', 0.8, 'Lambda', 1e-3, 'NearFieldSkipMm', 10);
        case 'oneil_phase_only'
            fo = struct('FreeParams', {{'phase'}}, 'ForwardModel', 'oneil', ...
                'Optimizer', 'fexminimize', 'Domain', 'intensity', 'HuberDelta', 10, ...
                'Lambda', 0, 'NearFieldPenalty', 0, 'NearFieldSkipMm', 10);
        otherwise
            error('calibration_fit_options:preset', ...
                'Unknown calibration preset ''%s''; use ''rayleigh_phase_gain'' or ''oneil_phase_only''.', char(preset));
    end

    fn = fieldnames(overrides);
    for k = 1:numel(fn); fo.(fn{k}) = overrides.(fn{k}); end
end
