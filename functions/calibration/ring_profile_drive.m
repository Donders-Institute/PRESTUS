function ann = ring_profile_drive(ann, elem_phase_deg, elem_p_amp)
% RING_PROFILE_DRIVE  Simulation drive of an annular transducer with a ring profile
%
% Library phases and drives are per channel (equal for the sub-rings of a
% ring); the ring profile's weights (annular.elem_weight, set by
% RING_PROFILE_TRANSDUCER) are applied by the analytical models. A simulated
% source has no such weights, so each sub-ring is driven with its channel's
% phase plus the weight's phase and its channel's drive times |weight|. The
% analytical fields (elem_weight, ring_profile) are removed, so nothing
% applies the weights twice. Without a profile, the phases and drive are set
% unchanged.
%
% Use as:
%   ann = ring_profile_drive(ann, elem_phase_deg, elem_p_amp)
%
% Input:
%   ann            - transducer.annular (possibly from RING_PROFILE_TRANSDUCER)
%   elem_phase_deg - [1 x N] library phases [deg] (LIBRARY_MODEL_AT_DEPTH)
%   elem_p_amp     - [1 x N] library drive [Pa], or [] for none
%
% Output:
%   ann - with elem_phase_deg, elem_phase_rad and elem_p_amp (if given) set
%
% See also: RING_PROFILE_TRANSDUCER, LOAD_TRANSDUCER_FROM_LIBRARY

    ph = reshape(double(elem_phase_deg), 1, []);
    pa = reshape(double(elem_p_amp), 1, []);
    if isfield(ann, 'elem_weight') && ~isempty(ann.elem_weight)
        u  = reshape(ann.elem_weight, 1, []);
        ph = mod(ph + angle(u) * 180 / pi, 360);
        if ~isempty(pa); pa = round(pa .* abs(u)); end
    end
    for f = {'elem_weight', 'ring_profile'}
        if isfield(ann, f{1}); ann = rmfield(ann, f{1}); end
    end
    ann.elem_phase_deg = ph;
    ann.elem_phase_rad = ph * pi / 180;
    if ~isempty(pa); ann.elem_p_amp = pa; end
end
