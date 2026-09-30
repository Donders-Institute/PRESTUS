function lib = library_entry_from_model(lib, model, isppa_wcm2, c_water, rho_water)
% LIBRARY_ENTRY_FROM_MODEL  Write a FIT_TRANSDUCER_MODEL result into a library struct
%
% Fills the static model (global_model) of a transducer library entry from a
% multi-depth fit. Existing refinements are kept; legacy per-depth phase
% tables are removed because the new static model supersedes them.
%
% Use as:
%   lib = library_entry_from_model(lib, model, isppa_wcm2)
%   lib = library_entry_from_model(lib, model, isppa_wcm2, c_water, rho_water)
%
% Input:
%   lib        - library struct (from yaml.loadFile) or struct() for a new entry
%   model      - output of FIT_TRANSDUCER_MODEL
%   isppa_wcm2 - reference ISPPA for amplitude_scaling / amplitude_by_depth [W/cm²]
%   c_water, rho_water - medium for velocity → source amplitude (default 1500, 994)
%
% Output:
%   lib - updated struct; write with yaml.dumpFile
%
% See also: FIT_TRANSDUCER_MODEL, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 4 || isempty(c_water);   c_water   = 1500; end
    if nargin < 5 || isempty(rho_water); rho_water = 994;  end
    if ~isfield(lib, 'global_model') || ~isstruct(lib.global_model); lib.global_model = struct(); end
    gm = lib.global_model;
    for legacy = {'depths_ep_mm', 'elem_phase_deg'}
        if isfield(gm, legacy{1}); gm = rmfield(gm, legacy{1}); end
    end

    % Source amplitude per depth that gives isppa_wcm2 at that depth in water:
    % drive_velocity reproduces the measured peak; scale as sqrt(I)
    amp_d = model.drive_velocity * rho_water * c_water .* sqrt(isppa_wcm2 ./ model.measured_peak_wcm2);

    gm.elem_phase_correction_deg = model.elem_phase_correction_deg;
    gm.elem_amp_weight           = model.elem_amp_weight;
    gm.steering_sound_speed      = model.steering_sound_speed;
    gm.forward_model             = model.forward_model;
    gm.free_params               = model.free_params;
    gm.fit_depths_ep_mm          = model.fit_depths_ep_mm;
    gm.fit_z_shift_mm            = model.z_shift_mm;
    gm.fit_r2                    = model.metrics.r2(:)';
    key = sprintf('i%s', strrep(num2str(isppa_wcm2), '.', 'p'));
    gm.amplitude_scaling = struct(key, round(median(amp_d)));
    gm.amplitude_by_depth = struct('depths_ep_mm', model.fit_depths_ep_mm, ...
        'elem_p_amp', round(amp_d), 'isppa_wcm2', isppa_wcm2);
    lib.global_model = gm;
end
