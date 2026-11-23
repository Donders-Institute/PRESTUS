function lib = library_refinement_from_fits(lib, fits, isppa_wcm2, blend_mm, source, c_water, rho_water)
% LIBRARY_REFINEMENT_FROM_FITS  Write per-depth FIT_TRANSDUCER_MODEL results as library refinements
%
% Each fit (one depth) contributes its phase offset, its per-ring gains (if
% 'amp' was a free parameter), its fitted z shift (diagnostic) and the drive that gives
% isppa_wcm2 in water at that depth with the refined phases — the same
% conversion LIBRARY_ENTRY_FROM_MODEL uses for the static model:
%   elem_p_amp = drive_velocity * rho * c * sqrt(isppa / measured_peak)
%
% Use as:
%   lib = library_refinement_from_fits(lib, fits, isppa_wcm2)
%   lib = library_refinement_from_fits(lib, fits, isppa_wcm2, blend_mm, source)
%
% Input:
%   lib        - library struct (from yaml.loadFile) or struct()
%   fits       - {1 x K} single-depth FIT_TRANSDUCER_MODEL results; if every
%                fit carries .profile_type ('measured' | 'interpolated' |
%                'extrapolated', set by the caller), it is stored per depth
%                as refinement.profile_type
%   isppa_wcm2 - reference ISPPA for the stored drives [W/cm²]
%   blend_mm   - refinement.blend_mm ([] = keep existing / 5)
%   source     - provenance string, e.g. the profile CSV
%   c_water, rho_water - medium for velocity → pressure (default 1500, 994)
%
% See also: LIBRARY_ADD_REFINEMENTS, REFINE_TRANSDUCER_LIBRARY, KWAVE_CALIBRATE_DRIVE

arguments
    lib
    fits cell
    isppa_wcm2 {mustBeNumeric}
    blend_mm = []
    source {mustBeTextScalar} = ''
    c_water {mustBeNumeric} = 1500
    rho_water {mustBeNumeric} = 994
end

    if nargin < 4; blend_mm = []; end
    if nargin < 5; source = ''; end
    if nargin < 6 || isempty(c_water);   c_water   = 1500; end
    if nargin < 7 || isempty(rho_water); rho_water = 994;  end

    depths = cellfun(@(f) f.fit_depths_ep_mm, fits);
    offs   = cellfun(@(f) f.elem_phase_correction_deg, fits, 'UniformOutput', false);
    amp_d  = cellfun(@(f) f.drive_velocity * rho_water * c_water * ...
                     sqrt(isppa_wcm2 / f.measured_peak_wcm2), fits);
    amp = struct('elem_p_amp', amp_d, 'isppa_wcm2', isppa_wcm2, ...
        'z_shift_mm', cellfun(@(f) f.z_shift_mm, fits));
    if all(cellfun(@(f) isfield(f, 'profile_type'), fits))          % set by callers fitting derived profiles
        amp.profile_type = cellfun(@(f) char(f.profile_type), fits, 'UniformOutput', false);
    end
    if any(cellfun(@(f) any(strcmp(f.free_params, 'amp')), fits))   % per-depth gains were fitted
        amp.elem_amp_weight = cellfun(@(f) f.elem_amp_weight, fits, 'UniformOutput', false);
    end
    lib = library_add_refinements(lib, depths, offs, blend_mm, fits{1}.forward_model, source, amp);
end
