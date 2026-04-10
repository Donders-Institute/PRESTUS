function color = risk_color(value, limit)
% RISK_COLOR  Returns 'green', 'amber', or 'red' based on value vs ITRUSST non-significant risk limit.
% Green: < 50% of limit. Amber: 50-100%. Red: exceeds limit.
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
