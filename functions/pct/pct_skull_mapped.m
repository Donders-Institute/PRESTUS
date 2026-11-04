function tf = pct_skull_mapped(parameters)
% PCT_SKULL_MAPPED  True when the skull layer's properties come from the pseudo-CT
%
% The skull is mapped from the pseudo-CT when pct.enabled is 1, unless the
% run is a bound variant of the uncertainty pipeline
% (uncertainty.homogeneous_skull): those share the pseudo-CT run's head
% geometry but give the skull layer the homogeneous properties of their
% medium file.
%
% Use as:
%   tf = pct_skull_mapped(parameters)
%
% See also: MEDIUM_SETUP, UNCERTAINTY_PIPELINE

    tf = isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && ...
         ~isempty(parameters.pct.enabled) && parameters.pct.enabled == 1 && ...
         ~(isfield(parameters, 'uncertainty') && isfield(parameters.uncertainty, 'homogeneous_skull') && ...
           parameters.uncertainty.homogeneous_skull);
end
