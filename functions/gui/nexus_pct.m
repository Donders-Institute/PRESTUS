function out = nexus_pct(args)
% NEXUS_PCT  Pseudo-CT steps of the PRESTUS Nexus, served by NEXUS_BRIDGE
%
% Use as:
%   out = nexus_pct(struct('action', 'predict', 'config', cfg))
%
% Input:
%   args - struct with .action and .config (path of the saved study config):
%          'predict'  PCT_PREDICT_SKULL_PROPERTIES: writes the predicted skull
%                     properties to the output's cache folder
%
% Output:
%   out - .file (written JSON), .source ('grid' | 'native'), .n_voxels
%
% See also: NEXUS_BRIDGE, PCT_PREDICT_SKULL_PROPERTIES

    p = load_parameters(char(args.config));
    switch char(args.action)
        case 'predict'
            s = pct_predict_skull_properties(p);
            out = struct('file', s.file, 'source', s.source, 'n_voxels', s.n_voxels);
        otherwise
            error('nexus_pct:action', 'Unknown action ''%s''.', char(args.action));
    end
end
