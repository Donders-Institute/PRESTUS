function flag = neuronav_series_qc_flag(r)
% NEURONAV_SERIES_QC_FLAG  Tracking problem of one detected series, as text
%
% Use as:
%   flag = neuronav_series_qc_flag(r)
%
% Input:
%   r - series statistics (NEURONAV_DETECT_SERIES .results{k})
%
% Output:
%   flag - 'untracked (all positions 0)', '<n> pulse(s) at origin', or '' when fine
%
% See also: NEURONAV_DETECT_SERIES, NEURONAV_NAME_SERIES

arguments
    r (1,1) struct
end

    P = reshape(r.matrix4d_all(:, 1:3, 4), [], 3);
    if all(abs(P(:)) < 1e-9)
        flag = 'untracked (all positions 0)';
    elseif any(all(abs(P) < 1e-9, 2))
        flag = sprintf('%d pulse(s) at origin', sum(all(abs(P) < 1e-9, 2)));
    else
        flag = '';
    end
end
