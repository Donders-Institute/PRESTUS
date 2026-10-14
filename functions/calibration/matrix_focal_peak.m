function T = matrix_focal_peak(parameters, tr, targets_mm, opts)
% MATRIX_FOCAL_PEAK  Modelled free-field focal peak of a matrix array per steering target
%
% For each target, steers the array geometrically (phase k |element - target|,
% as CREATE_MATRIX_KARRAY does) at a uniform element amplitude and returns
% the largest pressure amplitude in the plane through the target normal to
% the beam axis (z = target z), from the Rayleigh integral
% (MATRIX_RAYLEIGH_PRESSURE). This is the quantity a hydrophone scan of
% that plane measures, so measured / modelled peak is the steering
% correction (FIT_STEERING_CORRECTION). When the maximum sits on the edge
% of the scan window, the window is enlarged (up to three times).
%
% Use as:
%   T = matrix_focal_peak(parameters, tr, targets_mm)
%   T = matrix_focal_peak(parameters, tr, targets_mm, opts)
%
% Input:
%   parameters - PRESTUS config: medium_properties.water.sound_speed
%   tr         - matrix transducer struct (freq_hz, matrix.*; elem_p_amp(1)
%                is the element amplitude [Pa])
%   targets_mm - [N x 3] steering targets [mm], transducer frame (origin at
%                the apex / array centre, beam along +z)
%   opts       - optional: .half_width_mm (scan half-width, default 10),
%                .step_mm (default 0.5), .patch_per_lambda (default 6)
%
% Output:
%   T - table: target_x_mm, target_y_mm, target_z_mm, peak_pressure_pa,
%       peak_x_mm, peak_y_mm, on_edge (true if still on the window edge)
%
% See also: MATRIX_RAYLEIGH_PRESSURE, FIT_STEERING_CORRECTION, BUILD_STEERING_CORRECTION

    if nargin < 4; opts = struct(); end
    hw0 = get_opt(opts, 'half_width_mm', 10);
    step = get_opt(opts, 'step_mm', 0.5);
    ropts = struct('patch_per_lambda', get_opt(opts, 'patch_per_lambda', 6));
    c = parameters.medium_properties.water.sound_speed;
    f = tr.freq_hz;
    k = 2 * pi * f / c;
    L = matrix_element_layout(parameters, tr);
    amp = tr.matrix.elem_p_amp(1);

    N = size(targets_mm, 1);
    pk = zeros(N, 1); px = pk; py = pk; edge = false(N, 1);
    for i = 1:N
        tgt = targets_mm(i, :)' * 1e-3;
        phase = mod(k * vecnorm(L.pos_m - tgt, 2, 1), 2 * pi);
        hw = hw0;
        for attempt = 1:4
            g = -hw:step:hw;
            [GX, GY] = ndgrid(g, g);
            pts = [tgt(1) + GX(:)' * 1e-3; tgt(2) + GY(:)' * 1e-3; repmat(tgt(3), 1, numel(GX))];
            a = abs(matrix_rayleigh_pressure(L, f, c, amp, phase, pts, ropts));
            [pk(i), j] = max(a);
            on_edge = abs(abs(GX(j)) - hw) < step / 2 || abs(abs(GY(j)) - hw) < step / 2;
            if ~on_edge || attempt == 4; break; end
            hw = hw * 1.5;
        end
        px(i) = targets_mm(i, 1) + GX(j); py(i) = targets_mm(i, 2) + GY(j); edge(i) = on_edge;
    end
    T = table(targets_mm(:, 1), targets_mm(:, 2), targets_mm(:, 3), pk, px, py, edge, ...
        'VariableNames', {'target_x_mm', 'target_y_mm', 'target_z_mm', 'peak_pressure_pa', ...
        'peak_x_mm', 'peak_y_mm', 'on_edge'});
    if any(edge)
        warning('matrix_focal_peak:onEdge', '%d target(s): the peak is still on the scan window edge.', nnz(edge));
    end
end

function v = get_opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
