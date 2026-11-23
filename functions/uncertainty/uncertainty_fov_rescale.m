function t = uncertainty_fov_rescale(t, mismatch)
% UNCERTAINTY_FOV_RESCALE  Map a partial-FOV metrics table onto the whole grid
%
% Applies the factors of UNCERTAINTY_FOV_MISMATCH (default medium, whole grid
% against partial FOV) to the metrics table of a variant simulated in the
% partial FOV: 'ratio' metrics are multiplied by their factor, 'rise' metrics
% (absolute temperatures) keep their baseline and have the rise above it
% multiplied. Metrics of kind 'none', metrics without a valid factor, and
% columns the mismatch table does not list are returned unchanged.
%
% Use as:
%   t = uncertainty_fov_rescale(t, mismatch)
%
% Input:
%   t        - metrics table of a variant simulated in the partial FOV
%   mismatch - table from UNCERTAINTY_FOV_MISMATCH
%
% Output:
%   t        - the table with rescaled values (all rows)
%
% See also: UNCERTAINTY_FOV_MISMATCH, GENERATE_UNCERTAINTY_REPORT

arguments
    t        table
    mismatch table
end

for m = 1:height(mismatch)
    name = mismatch.metric{m};
    f    = mismatch.factor(m);
    if ~isfinite(f) || ~ismember(name, t.Properties.VariableNames) || ~isnumeric(t.(name))
        continue
    end
    switch mismatch.kind{m}
        case 'ratio'
            t.(name) = double(t.(name)) * f;
        case 'rise'
            b = mismatch.baseline(m);
            t.(name) = b + (double(t.(name)) - b) * f;
    end
end
end
