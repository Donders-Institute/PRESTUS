function color = risk_color(value, limit)
% RISK_COLOR  Return a colour string based on a value relative to an ITRUSST risk limit
%
% Returns 'green', 'amber', or 'red' based on value vs the ITRUSST
% non-significant risk limit: green < 50 %, amber 50–100 %, red > 100 %.
% Returns 'gray' for NaN/empty and 'info' for unlimited metrics.
%
% Use as:
%   color = risk_color(value, limit)
%
% Input:
%   value - numeric scalar to evaluate
%   limit - numeric risk limit (Inf = no limit)
%
% Output:
%   color - 'green', 'amber', 'red', 'gray', or 'info'
%
% See also: GET_RISK_LIMITS, GENERATE_SIMULATION_REPORT
    if isnan(value) || isempty(value)
        color = 'gray';
    elseif isinf(limit)
        color = 'info'; % informational, no limit
    elseif value > limit
        color = 'red';
    elseif value >= 0.5 * limit
        color = 'amber';
    else
        color = 'green';
    end
end
