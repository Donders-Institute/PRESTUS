function lib = library_set_ring_profile(lib, rp)
% LIBRARY_SET_RING_PROFILE  Set or remove the within-ring velocity profile of a library entry
%
% A ring profile (FIT_RING_PROFILE) splits every ring into sub-rings, so the
% per-element arrays of an entry refer to a different geometry once the
% profile is added, removed or changed. This function stores the profile in
% global_model.ring_profile ([] removes it) and drops what referred to the
% old geometry: the refinements (with a warning naming their depths) and the
% static per-element offsets and gains of the global model. An unchanged
% profile leaves the entry as it is.
%
% Use as:
%   lib = library_set_ring_profile(lib, rp)
%
% Input:
%   lib - library entry (struct)
%   rp  - ring profile struct, or [] for uniform rings
%
% Output:
%   lib - entry with the new profile and without stale per-element arrays
%
% See also: FIT_RING_PROFILE, LIBRARY_TRANSDUCER, CALIBRATE_LIBRARY_ENTRY, CALIBRATION_FIT_ENTRY

    if isempty(lib); lib = struct(); end
    old = [];
    if isfield(lib, 'global_model') && isfield(lib.global_model, 'ring_profile'); old = lib.global_model.ring_profile; end
    if isequal(old, rp); return; end
    if isfield(lib, 'refinement') && isfield(lib.refinement, 'depths_ep_mm') && ~isempty(lib.refinement.depths_ep_mm)
        warning('library_set_ring_profile:reset', ...
            ['The ring profile changes the element geometry: the refinements at %s mm are dropped ' ...
             'and only the depths calibrated now are stored.'], mat2str(double(lib.refinement.depths_ep_mm(:)')));
        lib = rmfield(lib, 'refinement');
    end
    if ~isfield(lib, 'global_model'); lib.global_model = struct(); end
    for f = {'elem_phase_correction_deg', 'elem_amp_weight'}
        if isfield(lib.global_model, f{1}); lib.global_model = rmfield(lib.global_model, f{1}); end
    end
    if isempty(rp)
        if isfield(lib.global_model, 'ring_profile'); lib.global_model = rmfield(lib.global_model, 'ring_profile'); end
    else
        lib.global_model.ring_profile = rp;
    end
end
