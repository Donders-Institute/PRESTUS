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
%          'curves'   PCT_MAPPING_CURVES: HU -> density, sound speed and
%                     attenuation of the configured mappings. .config may be
%                     '' (toolbox defaults); .pct and .freq_hz, when given,
%                     replace the config's (the Nexus's unsaved edits)
%
% Output:
%   out - predict: .file (written JSON), .source ('grid' | 'native'), .n_voxels
%         curves:  the PCT_MAPPING_CURVES struct
%
% See also: NEXUS_BRIDGE, PCT_PREDICT_SKULL_PROPERTIES

arguments
    args (1,1) struct
end

    if isfield(args, 'config') && ~isempty(args.config)
        p = load_parameters(char(args.config));
    else
        p = load_parameters();
    end
    switch char(args.action)
        case 'predict'
            s = pct_predict_skull_properties(p);
            out = struct('file', s.file, 'source', s.source, 'n_voxels', s.n_voxels);
        case 'curves'
            if isfield(args, 'pct') && isstruct(args.pct)
                for f = fieldnames(args.pct)'
                    p.pct.(f{1}) = args.pct.(f{1});
                end
            end
            if isfield(args, 'freq_hz') && ~isempty(args.freq_hz)
                p.transducer(1).freq_hz = args.freq_hz;
            end
            out = pct_mapping_curves(p);
        otherwise
            error('nexus_pct:action', 'Unknown action ''%s''.', char(args.action));
    end
end
