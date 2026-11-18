function warn(varargin)
% WARN  Issue a warning without the stack-trace printout
%
% Use as:
%   warn(msg, A1, ...)
%   warn(id, msg, A1, ...)
%
% Takes the same arguments as WARNING.
%
% See also: WARNING
    prev = warning('off', 'backtrace');
    warning(varargin{:});
    warning(prev);
end
