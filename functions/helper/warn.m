function warn(varargin)
% Wrapper around warning() that suppresses the stack-trace printout.
    prev = warning('off', 'backtrace');
    warning(varargin{:});
    warning(prev);
end
