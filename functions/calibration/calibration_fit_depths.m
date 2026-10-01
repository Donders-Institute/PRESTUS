function fits = calibration_fit_depths(tran, prof, fit_options, progress)
% CALIBRATION_FIT_DEPTHS  One FIT_TRANSDUCER_MODEL fit per depth of a profile set
%
% The per-depth calibration step shared by BUILD_TRANSDUCER_LIBRARY,
% CALIBRATE_LIBRARY_ENTRY (REFINE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER)
% and CALIBRATION_GUI. Each column of prof is fitted on its own; when prof
% carries profile_type (INTERPOLATE_AXIAL_PROFILES), each fit is tagged with
% it so LIBRARY_REFINEMENT_FROM_FITS can store it.
%
% Use as:
%   fits = calibration_fit_depths(tran, prof, fit_options)
%   fits = calibration_fit_depths(tran, prof, fit_options, progress)
%
% Input:
%   tran        - transducer struct from load_equipment_config (annular)
%   prof        - .depths_ep_mm [1 x N_d], .z_ep_mm, .I [N_z x N_d],
%                 optional .profile_type {1 x N_d}
%   fit_options - FIT_TRANSDUCER_MODEL options (e.g. CALIBRATION_FIT_OPTIONS)
%   progress    - optional @(k, n, depth_ep_mm, profile_type) called before each fit
%
% Output:
%   fits - {1 x N_d} FIT_TRANSDUCER_MODEL results
%
% See also: FIT_TRANSDUCER_MODEL, CALIBRATION_FIT_OPTIONS, LIBRARY_REFINEMENT_FROM_FITS

    if nargin < 4; progress = []; end
    n     = numel(prof.depths_ep_mm);
    typed = isfield(prof, 'profile_type');
    fits  = cell(1, n);
    for k = 1:n
        pt = 'measured';
        if typed; pt = prof.profile_type{k}; end
        if ~isempty(progress); progress(k, n, prof.depths_ep_mm(k), pt); end
        pk = struct('depths_ep_mm', prof.depths_ep_mm(k), 'z_ep_mm', prof.z_ep_mm, 'I', prof.I(:, k));
        fits{k} = fit_transducer_model(tran, pk, fit_options);
        if typed; fits{k}.profile_type = pt; end
    end
end
