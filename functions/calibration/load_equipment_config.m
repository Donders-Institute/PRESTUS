function equip = load_equipment_config(equipment_path)
% LOAD_EQUIPMENT_CONFIG  Load equipment configuration from individual YAML files
%
% Reads equipment_info.yaml for general settings, then loads all transducer
% YAMLs (type: transducer). Combos are embedded in each transducer file
% under a 'combos' section keyed by DS serial; the loader expands them into
% equip.combos keyed by '{tran_serial}_{ds_serial}'.
%
% Use as:
%   equip = load_equipment_config()
%   equip = load_equipment_config(equipment_path)
%
% Output:
%   equip.gen     — general settings (from equipment_info.yaml)
%   equip.trans   — struct keyed by tran_serial
%   equip.combos  — struct keyed by combo_name ({tran_serial}_{ds_serial})

    if nargin < 1 || isempty(equipment_path)
        equipment_path = fullfile(get_prestus_path(), 'config', 'equipment');
    end

    equip.gen    = yaml.loadFile(fullfile(equipment_path, 'equipment_info.yaml'), 'ConvertToArray', true);
    equip.trans  = struct();
    equip.combos = struct();

    files = dir(fullfile(equipment_path, '*.yaml'));
    for i = 1:numel(files)
        if strcmp(files(i).name, 'equipment_info.yaml')
            continue;
        end
        d = yaml.loadFile(fullfile(equipment_path, files(i).name), 'ConvertToArray', true);
        if ~isfield(d, 'type') || ~strcmp(d.type, 'transducer')
            continue;
        end

        equip.trans.(d.serial) = d;

        if ~isfield(d, 'combos') || isempty(fieldnames(d.combos))
            continue;
        end
        ds_keys = fieldnames(d.combos);
        for k = 1:numel(ds_keys)
            ds_serial  = char(ds_keys{k});
            combo_name = [char(d.serial), '_', ds_serial];
            entry = d.combos.(ds_serial);
            entry.tran_serial = d.serial;
            equip.combos.(combo_name) = entry;
        end
    end
end
