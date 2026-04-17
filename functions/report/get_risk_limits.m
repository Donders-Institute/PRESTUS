function limits = get_risk_limits(is_layered)
% GET_RISK_LIMITS  Returns struct with ITRUSST consensus non-significant risk limits (Aubry et al., 2025).
%
%   limits = get_risk_limits(true)   — tissue-specific limits (layered/phantom simulations)
%   limits = get_risk_limits(false)  — global-only limits (water/free-field simulations)
%
% Color thresholds: green < 50% of limit <= amber < limit <= red.

    arguments
        is_layered (1,1) logical = true
    end

    if ~is_layered
        limits = struct();
        limits.Isppa = struct('label', 'ISPPA (global)', 'limit', Inf, 'unit', 'W/cm²');
        limits.Psptp = struct('label', 'Psptp',          'limit', Inf, 'unit', 'Pa');
        return
    end

    limits = struct();

    % MI (ITRUSST: MI <= 1.9)
    limits.MI_tc    = struct('label', 'MI (transcranial)', 'limit', 1.9, 'unit', '');
    limits.MI_brain = struct('label', 'MI (brain)',        'limit', 1.9, 'unit', '');
    limits.MI_skull = struct('label', 'MI (skull)',        'limit', 1.9, 'unit', '');
    limits.MI_skin  = struct('label', 'MI (skin)',         'limit', 1.9, 'unit', '');

    % Temperature rise (ITRUSST: <= 2 °C)
    limits.riseT_brain = struct('label', 'Temp rise (brain)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skull = struct('label', 'Temp rise (skull)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skin  = struct('label', 'Temp rise (skin)',  'limit', 2.0, 'unit', [char(176) 'C']);

    % CEM43 (ITRUSST: brain <= 2, skull <= 16, skin <= 21 min)
    limits.CEM43_brain = struct('label', 'CEM43 (brain)', 'limit', 2.0,  'unit', 'min');
    limits.CEM43_skull = struct('label', 'CEM43 (skull)', 'limit', 16.0, 'unit', 'min');
    limits.CEM43_skin  = struct('label', 'CEM43 (skin)',  'limit', 21.0, 'unit', 'min');

    % Absolute temperature (ITRUSST: <= 39 °C)
    limits.maxT_brain = struct('label', 'Max temp (brain)', 'limit', 39.0, 'unit', [char(176) 'C']);
    limits.maxT_skull = struct('label', 'Max temp (skull)', 'limit', 39.0, 'unit', [char(176) 'C']);
    limits.maxT_skin  = struct('label', 'Max temp (skin)',  'limit', 39.0, 'unit', [char(176) 'C']);

    % ISPPA (informational, no ITRUSST limit)
    limits.Isppa_brain = struct('label', 'ISPPA (brain)', 'limit', Inf, 'unit', 'W/cm²');
    limits.Isppa_skull = struct('label', 'ISPPA (skull)', 'limit', Inf, 'unit', 'W/cm²');
    limits.Isppa_skin  = struct('label', 'ISPPA (skin)',  'limit', Inf, 'unit', 'W/cm²');
end
