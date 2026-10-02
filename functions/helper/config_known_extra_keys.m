function keys = config_known_extra_keys(prestus_root)
% CONFIG_KNOWN_EXTRA_KEYS  Config keys PRESTUS reads that config_default.yaml omits
%
% Returns the dotted keys listed under 'extra:' and 'deprecated:' in the PRESTUS
% Nexus overlay (nexus/prestus_nexus/overlay.yaml). The overlay is the single
% list of such keys: the Nexus drift test (nexus/tests/test_schema.py) checks
% that each one is still read by functions/, so MATLAB and the Nexus agree on
% which study-config keys are valid.
%
% The overlay's keys contain dots, which MATLAB struct fields cannot hold, so
% the two sections are read as text (two-space-indented 'key:' lines).
%
% Use as:
%   keys = config_known_extra_keys()
%   keys = config_known_extra_keys(prestus_root)
%
% Output:
%   keys - cellstr of dotted keys ({} when the overlay is not found)
%
% See also: CONFIG_UNKNOWN_KEYS, LOAD_PARAMETERS

    if nargin < 1 || isempty(prestus_root)
        prestus_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    end
    keys = {};
    overlay = fullfile(prestus_root, 'nexus', 'prestus_nexus', 'overlay.yaml');
    if ~isfile(overlay)
        return
    end
    lines = splitlines(fileread(overlay));
    section = '';
    for i = 1:numel(lines)
        ln = lines{i};
        top = regexp(ln, '^([A-Za-z_]\w*):', 'tokens', 'once');
        if ~isempty(top)
            section = top{1};
            continue
        end
        if any(strcmp(section, {'extra', 'deprecated'}))
            tok = regexp(ln, '^  ([A-Za-z_][\w.]*):', 'tokens', 'once');
            if ~isempty(tok)
                keys{end+1} = tok{1}; %#ok<AGROW>
            end
        end
    end
    keys = unique(keys, 'stable');
end
