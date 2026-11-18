function equip = load_equipment_config(equipment_path, user_folders)
% LOAD_EQUIPMENT_CONFIG  Load the equipment library (transducers and driving systems)
%
% Reads equipment_info.yaml for general settings, then every YAML file in
% the folder and its subfolders:
%   type: transducer      a device model or measured unit
%                         (shipped: model-<MODEL>/tran-<SERIAL>_transducer.yaml;
%                         a flat <SERIAL>.yaml is read as well)
%   type: driving_system  a driving system (shipped: ds/ds-<DS>_driver.yaml)
% Combos are embedded in each transducer entry under 'combos', keyed by DS
% serial; the loader expands them into equip.combos keyed by
% '{tran_serial}_{ds_serial}' and fills a combo's 'ds' block from the
% driving-system entry when the combo has none.
%
% User folders (your own libraries, see TRANSDUCER_LIBRARY_FOLDERS; most
% specific first) may hold equipment entries next to calibration entries:
% files named *_transducer.yaml or *_driver.yaml, in the folder or its
% subfolders. They add serials or replace shipped entries of the same serial;
% an earlier folder in the list replaces a later one. A user unit may name a
% shipped model.
%
% A unit (model: <MODEL>, another entry) takes the fields it omits from its
% model entry: transducer geometry, min_foc/max_foc (the hardware steering
% range), n_elem and manufact.
%
% Use as:
%   equip = load_equipment_config()
%   equip = load_equipment_config(equipment_path)
%   equip = load_equipment_config(equipment_path, user_folders)   % char or cellstr
%
% Output:
%   equip.gen     — general settings (from equipment_info.yaml; empty struct without it)
%   equip.trans   — struct keyed by tran_serial; each model entry comes before
%                   its units (a unit's 'model', else its model-<MODEL> folder)
%   equip.ds      — struct keyed by driving-system serial
%   equip.combos  — struct keyed by combo_name ({tran_serial}_{ds_serial})
%   equip.files   — struct keyed by tran_serial: the file each entry came from
%   equip.user    — struct keyed by tran_serial: true for entries from a user folder

    if nargin < 1 || isempty(equipment_path)
        equipment_path = fullfile(get_prestus_path(), 'config', 'equipment');
    end

    equip.path   = equipment_path;
    info = fullfile(equipment_path, 'equipment_info.yaml');
    equip.gen    = struct();     % a folder without it (scratch folder of a single entry) has no general settings
    if isfile(info); equip.gen = yaml.loadFile(info, 'ConvertToArray', true); end
    equip.trans  = struct();
    equip.ds     = struct();
    equip.combos = struct();
    equip.files  = struct();
    equip.user   = struct();
    if nargin < 2 || isempty(user_folders); user_folders = {}; end
    user_folders = cellstr(user_folders);

    files = dir(fullfile(equipment_path, '**', '*.yaml'));
    files = files(~strcmp({files.name}, 'equipment_info.yaml'));
    is_user = false(numel(files), 1);
    for u = numel(user_folders):-1:1          % most specific last, so it wins
        if isempty(user_folders{u}) || ~isfolder(user_folders{u}) || ...
                strcmp(canonical(user_folders{u}), canonical(equipment_path))
            continue;
        end
        % only named equipment files: the folder also holds calibration entries
        uf = [dir(fullfile(user_folders{u}, '**', '*_transducer.yaml')); ...
              dir(fullfile(user_folders{u}, '**', '*_driver.yaml'))];
        files = [files(:); uf(:)];
        is_user = [is_user; true(numel(uf), 1)]; %#ok<AGROW>
    end
    for i = 1:numel(files)
        fpath = fullfile(files(i).folder, files(i).name);
        d = yaml.loadFile(fpath, 'ConvertToArray', true);
        if ~isstruct(d) || ~isfield(d, 'type') || ~isfield(d, 'serial')
            continue;
        end
        serial = char(d.serial);
        switch char(d.type)
            case 'transducer'
                if is_user(i) && isfield(equip.trans, serial)
                    fprintf('Equipment: %s from %s replaces %s.\n', serial, fpath, equip.files.(serial));
                end
                equip.trans.(serial) = d;
                equip.files.(serial) = fpath;
                equip.user.(serial)  = is_user(i);
            case 'driving_system'
                equip.ds.(serial) = d;
        end
    end

    % Units take what they omit from their model entry.
    inherited = {'transducer', 'min_foc', 'max_foc', 'n_elem', 'manufact'};
    serials = fieldnames(equip.trans);
    for i = 1:numel(serials)
        d = equip.trans.(serials{i});
        if ~isfield(d, 'model') || isempty(d.model) || strcmp(char(d.model), serials{i}) || ...
                ~isfield(equip.trans, char(d.model))
            continue;
        end
        m = equip.trans.(char(d.model));
        for f = inherited
            if (~isfield(d, f{1}) || isempty(d.(f{1}))) && isfield(m, f{1})
                d.(f{1}) = m.(f{1});
            end
        end
        equip.trans.(serials{i}) = d;
    end

    % Each model before its units, so listings start with the generic entry
    [equip.trans, perm] = orderfields(equip.trans, models_first(equip.trans, equip.files));
    serials = fieldnames(equip.trans);
    equip.files = orderfields(equip.files, perm);
    equip.user  = orderfields(equip.user, perm);

    % Combos, with the driving system described by its own entry
    for i = 1:numel(serials)
        d = equip.trans.(serials{i});
        if ~isfield(d, 'combos') || ~isstruct(d.combos) || isempty(fieldnames(d.combos))
            continue;
        end
        ds_keys = fieldnames(d.combos);
        for k = 1:numel(ds_keys)
            ds_serial  = char(ds_keys{k});
            entry = d.combos.(ds_serial);
            if ~isstruct(entry); entry = struct(); end
            if ~isfield(entry, 'ds') && isfield(equip.ds, ds_serial)
                ds = equip.ds.(ds_serial);
                entry.ds = rmfield(ds, intersect(fieldnames(ds), {'type'}));
            end
            d.combos.(ds_serial) = entry;
            entry.tran_serial = d.serial;
            equip.combos.([char(d.serial), '_', ds_serial]) = entry;
        end
        equip.trans.(serials{i}) = d;
    end
end

function order = models_first(trans, files)
% Field order listing each model entry before its units. A device belongs to
% its 'model' entry, else to its shipped model-<MODEL> folder's entry, else
% to itself; groups keep the order in which their first entry loaded.
    serials = fieldnames(trans);
    group = serials;
    for i = 1:numel(serials)
        d = trans.(serials{i});
        [~, folder] = fileparts(fileparts(files.(serials{i})));
        if isfield(d, 'model') && ~isempty(d.model) && isfield(trans, char(d.model))
            group{i} = char(d.model);
        elseif startsWith(folder, 'model-') && isfield(trans, extractAfter(folder, 'model-'))
            group{i} = extractAfter(folder, 'model-');
        end
    end
    [~, first] = unique(group, 'stable');
    [~, g] = ismember(group, group(sort(first)));
    is_unit = ~strcmp(group, serials);
    [~, order] = sortrows([g(:), double(is_unit(:)), (1:numel(serials))']);
end

function p = canonical(p)
    p = char(java.io.File(char(p)).getCanonicalPath());
end
