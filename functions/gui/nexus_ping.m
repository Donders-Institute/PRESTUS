function out = nexus_ping(args)
% NEXUS_PING  Answer a NEXUS_BRIDGE health check
%
% Use as:
%   out = nexus_ping(args)
%
% Input:
%   args - struct; echoed back as out.echo
%
% Output:
%   out - .matlab (version), .prestus_root, .echo
%
% See also: NEXUS_BRIDGE

arguments
    args (1,1) struct
end

    out = struct('matlab', version, 'prestus_root', get_prestus_path(), 'echo', args);
end
