function preset = babelbrain_preset(bb_root, tx_system)
% BABELBRAIN_PRESET  Resolve a BabelBrain preset from BabelBrain's transducer registry
%
% Looks the preset up in BabelBrain/SelFiles/transducer_list.yaml (BabelBrain
% 0.8.8 and later) and returns its geometry family and default.yaml, the way
% BabelBrain itself locates a device. Mirrors resolve_preset in
% run_babelbrain_case.py. A renamed preset (REMOPD, now TheraPhas) resolves to
% its current name.
%
% Use as:
%   preset = babelbrain_preset(bb_root, tx_system)
%
% Input:
%   bb_root   - BabelBrain source root (contains TranscranialModeling/)
%   tx_system - BabelBrain preset name, e.g. 'CTX_500'
%
% Output:
%   preset - struct with fields
%       .name - current BabelBrain preset name
%       .type - geometry family (simple_focused, focused_annular_array,
%               flat_annular_array, focused_array, flat_array_2D)
%       .yaml - path to the preset default.yaml
%       .cfg  - contents of default.yaml (empty struct if unreadable)
%
% See also: BABELBRAIN_TX_GEOMETRY, BABELBRAIN_WRAPPER, EXPORT_BABELBRAIN_INPUTS

    arguments
        bb_root   (1,:) char
        tx_system (1,:) char
    end

    aliases = struct('REMOPD', 'TheraPhas');
    name = tx_system;
    if isfield(aliases, name)
        name = aliases.(name);
    end

    reg_path = fullfile(bb_root, 'BabelBrain', 'SelFiles', 'transducer_list.yaml');
    if ~isfile(reg_path)
        error('PRESTUS:babelbrain:noRegistry', ...
            ['BabelBrain transducer registry not found: %s. PRESTUS needs BabelBrain 0.8.8 ' ...
             'or later (git submodule update --init external/BabelBrain).'], reg_path);
    end
    reg = yaml.loadFile(reg_path);
    if isstruct(reg)
        reg = num2cell(reg);
    end

    entry = [];
    for k = 1:numel(reg)
        if strcmp(char(reg{k}.name), name)
            entry = reg{k};
            break
        end
    end
    if isempty(entry)
        names = cellfun(@(e) char(e.name), reg, 'UniformOutput', false);
        error('PRESTUS:babelbrain:unknownPreset', ...
            'Unknown BabelBrain preset ''%s''. Known presets: %s.', ...
            tx_system, strjoin(sort(names), ', '));
    end

    preset.name = name;
    preset.type = char(entry.transducer_type);
    preset.yaml = fullfile(bb_root, 'BabelBrain', 'babel_transducers', preset.type, ...
        char(entry.module_name), 'default.yaml');
    preset.cfg  = struct();
    if isfile(preset.yaml)
        try
            preset.cfg = yaml.loadFile(preset.yaml, 'ConvertToArray', true);
        catch ME
            warning('PRESTUS:babelbrain:presetYaml', ...
                'Could not read BabelBrain preset yaml %s: %s', preset.yaml, ME.message);
        end
    end
end
