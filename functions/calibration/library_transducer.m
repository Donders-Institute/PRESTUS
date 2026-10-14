function tran = library_transducer(tran, lib)
% LIBRARY_TRANSDUCER  Transducer geometry a library entry was calibrated with
%
% Library entries are calibrated on the equipment geometry, unless the entry
% carries a within-ring velocity profile (global_model.ring_profile,
% FIT_RING_PROFILE): then every ring is split into sub-rings with fixed
% weights (RING_PROFILE_TRANSDUCER) and the entry's per-element arrays refer
% to those sub-rings. Every function that combines a library entry with the
% equipment geometry resolves the geometry through this function, so an
% entry with a profile is modelled, fitted and simulated consistently.
% Idempotent: a transducer that already carries a profile is returned as is.
%
% Use as:
%   tran = library_transducer(tran, lib)
%
% Input:
%   tran - transducer struct from load_equipment_config
%   lib  - library entry (yaml.loadFile(..., 'ConvertToArray', true)), or []
%
% Output:
%   tran - tran, with the entry's ring profile applied if it has one
%
% See also: RING_PROFILE_TRANSDUCER, FIT_RING_PROFILE, LIBRARY_MODEL_AT_DEPTH

    rp = library_ring_profile(lib);
    if ~isempty(rp)
        tran = ring_profile_transducer(tran, rp);
    end
end

function rp = library_ring_profile(lib)
    rp = [];
    if isstruct(lib) && isfield(lib, 'global_model') && isstruct(lib.global_model) && ...
            isfield(lib.global_model, 'ring_profile') && isstruct(lib.global_model.ring_profile)
        rp = lib.global_model.ring_profile;
    end
end
