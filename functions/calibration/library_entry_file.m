function [file, found, index] = library_entry_file(library_path, key, tran_serial, ds_serial, model)
% LIBRARY_ENTRY_FILE  Path of a calibration-library entry in a flat or nested library
%
% A library folder holds entries in one of two layouts (both are read):
%   flat     <library_path>/<KEY>.yaml                      (user libraries)
%   nested   <library_path>/model-<MODEL>/tran-<TRAN>[_ds-<DS>]_calib.yaml
% where KEY is <TRAN>_<DS> for a transducer–driving-system entry and <TRAN>
% for a driving-system-agnostic one. The shipped library (config/calibration)
% is nested.
%
% Use as:
%   [file, found] = library_entry_file(library_path, key)
%   file = library_entry_file(library_path, key, tran_serial, ds_serial, model)
%   [~, ~, index] = library_entry_file(library_path)     % all entries
%
% Inputs:
%   library_path - library folder
%   key          - entry key (<TRAN>_<DS> or <TRAN>)
%   tran_serial, ds_serial, model - optional; place a new entry in a nested
%                  library ('' ds_serial = driving-system agnostic). Without
%                  them, or in a flat library, a new entry is written flat.
%
% Outputs:
%   file  - the existing entry; when none exists, where to write it
%   found - true when the entry exists
%   index - struct array with fields key, file (every entry in the folder)
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, LOAD_EQUIPMENT_CONFIG

    if nargin < 3; tran_serial = ''; end
    if nargin < 4; ds_serial = ''; end
    if nargin < 5; model = ''; end
    library_path = char(library_path);

    index = list_entries(library_path);
    file = ''; found = false;
    if nargin < 2 || isempty(key)
        return;
    end
    key = char(key);
    hit = find(strcmp({index.key}, key), 1);
    if ~isempty(hit)
        file = index(hit).file; found = true;
        return;
    end

    % New entry: nested when the folder is nested and the transducer is known
    file = fullfile(library_path, [key '.yaml']);
    if isempty(tran_serial)
        return;
    end
    models = dir(fullfile(library_path, 'model-*'));
    models = models([models.isdir]);
    if isempty(models)
        return;
    end
    tran_serial = char(tran_serial);
    if isempty(model)
        try     % the model named in the transducer's equipment entry
            eq = load_equipment_config();
            if isfield(eq.trans, tran_serial) && isfield(eq.trans.(tran_serial), 'model') && ...
                    ~isempty(eq.trans.(tran_serial).model)
                model = char(eq.trans.(tran_serial).model);
            end
        catch
        end
    end
    if isempty(model)
        % the longest model folder whose name prefixes the transducer serial
        names = cellfun(@(n) n(7:end), {models.name}, 'UniformOutput', false);
        match = names(cellfun(@(n) strcmp(tran_serial, n) || startsWith(tran_serial, [n '_']), names));
        if isempty(match)
            model = tran_serial;
        else
            [~, k] = max(cellfun(@numel, match));
            model = match{k};
        end
    end
    stem = ['tran-' tran_serial];
    if ~isempty(ds_serial)
        stem = [stem '_ds-' char(ds_serial)];
    end
    file = fullfile(library_path, ['model-' char(model)], [stem '_calib.yaml']);
end

function index = list_entries(library_path)
    index = struct('key', {}, 'file', {});
    flat = dir(fullfile(library_path, '*.yaml'));
    for i = 1:numel(flat)
        index(end+1) = struct('key', flat(i).name(1:end-5), ...
            'file', fullfile(library_path, flat(i).name)); %#ok<AGROW>
    end
    nested = dir(fullfile(library_path, 'model-*', 'tran-*_calib.yaml'));
    for i = 1:numel(nested)
        stem = nested(i).name(6:end-11);              % strip 'tran-' and '_calib.yaml'
        parts = strsplit(stem, '_ds-');
        key = parts{1};
        if numel(parts) > 1
            key = [key '_' strjoin(parts(2:end), '_ds-')];
        end
        index(end+1) = struct('key', key, 'file', fullfile(nested(i).folder, nested(i).name)); %#ok<AGROW>
    end
end
