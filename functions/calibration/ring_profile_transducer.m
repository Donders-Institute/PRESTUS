function tran = ring_profile_transducer(tran, rp)
% RING_PROFILE_TRANSDUCER  Annular transducer with a within-ring velocity profile
%
% A ring profile describes how each ring of an annular transducer vibrates
% across its width: every ring is split into rp.n_sub sub-rings of equal
% width, each with a fixed complex weight (rp.weight_abs, rp.weight_phase_deg;
% the same at every focal depth, i.e. a property of the hardware). The driver
% still sets one phase and gain per channel; the weights shape the field of
% each ring (FIT_RING_PROFILE). Uniform rings are rp.weight_abs = 1 and
% rp.weight_phase_deg = 0.
%
% The returned transducer has the split geometry (annular.elem_n = rings x
% n_sub, elem_id_mm / elem_od_mm per sub-ring), the complex weights in
% annular.elem_weight (applied by AXIAL_TRANSFER_MATRIX), the steering
% radius of the parent ring in annular.steer_r_mid_mm (COMPUTE_PHASES,
% STEERING_PATH_LENGTHS: a ring is steered as a whole) and the channel count
% in n_elem, so FIT_TRANSDUCER_MODEL groups the sub-rings of a channel.
% A transducer that already carries a profile is returned unchanged.
%
% Use as:
%   tran = ring_profile_transducer(tran, rp)
%
% Input:
%   tran - transducer struct from load_equipment_config (annular)
%   rp   - ring profile: .n_sub, .weight_abs, .weight_phase_deg ([1 x rings*n_sub],
%          or [] for uniform); [] or n_sub <= 1 returns tran unchanged
%
% Output:
%   tran - transducer with split geometry, .transducer.annular.elem_weight,
%          .transducer.annular.ring_profile and .ring_profile_applied = true
%
% See also: LIBRARY_TRANSDUCER, FIT_RING_PROFILE, AXIAL_TRANSFER_MATRIX,
%           LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 2 || isempty(rp) || ~isstruct(rp) || ~isfield(rp, 'n_sub') || double(rp.n_sub) <= 1
        return;
    end
    if isfield(tran, 'ring_profile_applied') && tran.ring_profile_applied
        return;
    end
    a  = tran.transducer.annular;
    M  = double(rp.n_sub);
    nr = double(a.elem_n);
    if ~isfield(tran, 'n_elem') || isempty(tran.n_elem); tran.n_elem = nr; end   % channels

    ids = zeros(1, nr * M); ods = ids;
    for i = 1:nr
        e = linspace(double(a.elem_id_mm(i)), double(a.elem_od_mm(i)), M + 1);
        ids((i - 1) * M + (1:M)) = e(1:M);
        ods((i - 1) * M + (1:M)) = e(2:M + 1);
    end
    wa = field_or(rp, 'weight_abs', ones(1, nr * M));
    wp = field_or(rp, 'weight_phase_deg', zeros(1, nr * M));
    wa = double(wa(:)'); wp = double(wp(:)');
    if numel(wa) ~= nr * M || numel(wp) ~= nr * M
        error('ring_profile_transducer:length', ...
            'Ring profile has %d / %d weights; %d rings x %d sub-rings = %d expected.', ...
            numel(wa), numel(wp), nr, M, nr * M);
    end
    a.elem_n      = nr * M;
    a.elem_id_mm  = ids;
    a.elem_od_mm  = ods;
    a.elem_weight = wa .* exp(1i * wp * pi / 180);
    % The driver steers a ring as a whole: every sub-ring takes its ring's steering phase
    a.steer_r_mid_mm = kron((double(tran.transducer.annular.elem_id_mm(:)') / 2 + ...
        double(tran.transducer.annular.elem_od_mm(:)') / 2) / 2, ones(1, M));
    a.ring_profile = rp;
    tran.transducer.annular = a;
    tran.ring_profile_applied = true;
end

function v = field_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
