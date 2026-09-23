function limits = get_risk_limits(is_layered)
% GET_RISK_LIMITS  Return ITRUSST consensus non-significant risk limits
%
% Returns a struct of acoustic and thermal safety limits based on the
% ITRUSST consensus (Aubry et al., 2025). Tissue-specific limits are
% returned for layered/phantom simulations; global-only limits for
% water/free-field simulations.
%
% Use as:
%   limits = get_risk_limits(is_layered)
%
% Input:
%   is_layered - true for tissue-specific limits, false for global-only (default: true)
%
% Output:
%   limits - struct with fields MI_tc, MI_brain, Isppa, Ispta, etc.
%
% See also: GENERATE_SIMULATION_REPORT, RISK_COLOR

    arguments
        is_layered (1,1) logical = true
    end

    if ~is_layered
        limits = struct();
        limits.Isppa = struct('label', 'ISPPA (global)',   'limit', Inf, 'unit', 'W/cm²');
        limits.MI    = struct('label', 'MI (free water)',  'limit', 1.9, 'unit', '');
        % MI_tc is layered-only data, but the field is always defined so the
        % Safety board can show both MI tiles side by side; with no MI_tc
        % column in a water-medium CSV it naturally renders as N/A (gray).
        limits.MI_tc = struct('label', 'MI (transcranial)', 'limit', 1.9, 'unit', '');
        limits.Psptp = struct('label', 'Max pressure',     'limit', 2e6, 'unit', 'Pa');
        limits.TIC   = struct('label', 'Cranial TI (TIC)', 'limit', Inf, 'unit', '');
        return
    end

    limits = struct();

    % MI (ITRUSST: MI <= 1.9)
    limits.MI_tc    = struct('label', 'MI (transcranial)', 'limit', 1.9, 'unit', '');
    % MI (free water) is water-medium-only data, but the field is always
    % defined so the Safety board can show both MI tiles side by side; with
    % no MI column in a layered-medium CSV it naturally renders as N/A (gray)
    % — i.e. MI (free water) is effectively grayed out whenever MItc is available.
    limits.MI       = struct('label', 'MI (free water)',   'limit', 1.9, 'unit', '');
    limits.MI_brain = struct('label', 'MI (brain)',        'limit', 1.9, 'unit', '');
    limits.MI_skull = struct('label', 'MI (skull)',        'limit', 1.9, 'unit', '');
    limits.MI_skin  = struct('label', 'MI (skin)',         'limit', 1.9, 'unit', '');

    % Peak pressure (kept below 2 MPa)
    limits.Psptp = struct('label', 'Max pressure', 'limit', 2e6, 'unit', 'Pa');

    % Cranial Thermal Index (IEC 62359) — informational, no hard ITRUSST limit
    limits.TIC      = struct('label', 'Cranial TI (TIC)',  'limit', Inf, 'unit', '');

    % Temperature rise (ITRUSST: <= 2 °C)
    limits.riseT_brain = struct('label', 'Temp rise (brain)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skull = struct('label', 'Temp rise (skull)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skin  = struct('label', 'Temp rise (skin)',  'limit', 2.0, 'unit', [char(176) 'C']);

    % Absolute temperature rise from a fixed 37 °C start (ITRUSST: <= 2 °C)
    limits.riseT37_brain = struct('label', [char(916) 'T from 37' char(176) 'C (brain)'], 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT37_skull = struct('label', [char(916) 'T from 37' char(176) 'C (skull)'], 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT37_skin  = struct('label', [char(916) 'T from 37' char(176) 'C (skin)'],  'limit', 2.0, 'unit', [char(176) 'C']);

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
