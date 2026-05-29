function result = is_async_mode(parameters)
% IS_ASYNC_MODE  Return true when async multi-transducer mode is active.
%
% Async mode is active when ALL of the following hold:
%   1. parameters.simulation.transducer_coupling equals 'async'
%   2. parameters.transducer contains more than one entry
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
end
