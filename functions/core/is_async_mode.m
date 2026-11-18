function result = is_async_mode(parameters)
% IS_ASYNC_MODE  Return true when async multi-transducer mode is active.
%
% Async mode is active when ALL of the following hold:
%   1. parameters.simulation.transducer_coupling equals 'async'
%   2. parameters.transducer contains more than one entry
% or for coherent transducers at different frequencies
% (IS_SIMULTANEOUS_MULTIFREQ), which run through the same pipeline and are
% combined as firing together.
%
% Use as:
%   tf = is_async_mode(parameters)
%
% See also: ASYNC_TRANSDUCER_PIPELINE, PRESTUS_PIPELINE_START

result = isfield(parameters, 'simulation') && ...
         isfield(parameters.simulation, 'transducer_coupling') && ...
         strcmp(parameters.simulation.transducer_coupling, 'async') && ...
         isfield(parameters, 'transducer') && ...
         numel(parameters.transducer) > 1;
result = result || is_simultaneous_multifreq(parameters);
end
