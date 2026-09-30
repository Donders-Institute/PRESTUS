function [manifest, info] = study_read_manifest(manifest_path)
% STUDY_READ_MANIFEST  Read a study manifest: one row per simulation run
%
% A manifest is a tab- (.tsv) or comma-separated (.csv) table. Each row is
% one run: a subject, optionally under a named condition, with its own
% settings. Reserved columns:
%
%   subject_id  (required) integer, '003' or 'sub-003'
%   condition   (optional) run label, e.g. target or protocol name; becomes
%               the run's output affix (see STUDY_RUN_TABLE) and the
%               grouping variable in the group report
%   include     (optional) 1/0, true/false, yes/no; default 1
%   config      (optional) per-row YAML merged over the study config; path
%               relative to the manifest folder
%   notes       (optional) free text, ignored (also: note, comment)
%
% Every other column is a config key in dotted form (e.g.
% 'placement.localite.target', 'transducer.target_isppa_wcm2',
% 'transducer(2).serial'). Cells are parsed as YAML scalars or flow
% sequences ('30', '[-38, 20, 52]', 'L_MD', 'true'); an empty cell keeps the
% study value, '[]' sets an empty value. Lines starting with '#' are
% comments.
%
% Use as:
%   [manifest, info] = study_read_manifest(manifest_path)
%
% Input:
%   manifest_path - path to the .tsv/.csv manifest
%
% Output:
%   manifest - table, one row per manifest row, with columns
%              row (manifest row number), subject_id (double), condition
%              (string), include (logical), config (string, absolute or ""),
%              overrides (cell; each an N×2 {key, value} cell)
%   info     - struct: path, dir, override_keys (cellstr, column order)
%
% See also: STUDY_RUN_TABLE, STUDY_ROW_PARAMETERS, PRESTUS_STUDY_START,
%           APPLY_CONFIG_OVERRIDES

arguments
    manifest_path {mustBeTextScalar}
end

    manifest_path = char(manifest_path);
    if ~isfile(manifest_path)
        error('study_read_manifest:notFound', 'Study manifest not found: %s', manifest_path);
    end
    [mdir, ~, ext] = fileparts(manifest_path);
    if isempty(mdir); mdir = pwd; end
    switch lower(ext)
        case '.tsv', delim = '\t';
        case '.csv', delim = ',';
        otherwise
            error('study_read_manifest:badExtension', ...
                'Manifest must be a .tsv or .csv file: %s', manifest_path);
    end

    raw = read_text_table(manifest_path, delim);
    cols = raw.Properties.VariableNames;
    cols_lc = lower(strtrim(cols));
    if ~ismember('subject_id', cols_lc)
        error('study_read_manifest:noSubjectColumn', ...
            'Manifest %s has no subject_id column.', manifest_path);
    end

    reserved = {'subject_id', 'condition', 'include', 'config', 'notes', 'note', 'comment'};
    is_key   = ~ismember(cols_lc, reserved);
    key_cols = cols(is_key);
    for c = 1:numel(key_cols)
        if isempty(regexp(strtrim(key_cols{c}), '^[A-Za-z]\w*(\(\d+\))?(\.[A-Za-z]\w*(\(\d+\))?)*$', 'once'))
            error('study_read_manifest:badColumn', ...
                ['Manifest column "%s" is neither a reserved column (%s) nor a dotted ' ...
                 'config key.'], key_cols{c}, strjoin(reserved(1:4), ', '));
        end
    end

    n = height(raw);
    subject_id = nan(n, 1);
    condition  = strings(n, 1);
    include    = true(n, 1);
    config     = strings(n, 1);
    overrides  = cell(n, 1);
    get = @(r, name) strtrim(raw{r, cols{strcmp(cols_lc, name)}});

    for r = 1:n
        subject_id(r) = parse_subject_id(get(r, 'subject_id'), r);
        if ismember('condition', cols_lc); condition(r) = get(r, 'condition'); end
        if ismember('include', cols_lc);   include(r)   = parse_include(get(r, 'include'), r); end
        if ismember('config', cols_lc)
            cfg = get(r, 'config');
            if strlength(cfg) > 0
                cfg = char(cfg);
                if ~is_absolute(cfg); cfg = fullfile(mdir, cfg); end
                config(r) = string(cfg);
            end
        end
        ov = cell(0, 2);
        for c = 1:numel(key_cols)
            txt = strtrim(raw{r, key_cols{c}});
            if strlength(txt) == 0; continue; end   % empty cell: keep the study value
            ov(end+1, :) = {strtrim(key_cols{c}), parse_manifest_value(txt)}; %#ok<AGROW>
        end
        overrides{r} = ov;
    end

    row = (1:n)';
    manifest = table(row, subject_id, condition, include, config, overrides);
    info = struct('path', manifest_path, 'dir', mdir, 'override_keys', {strtrim(key_cols)});
end

% -------------------------------------------------------------------------
function T = read_text_table(file, delim)
% All columns as text, header names preserved (dots, parentheses), '#' comments.
    lines = splitlines(string(fileread(file)));
    keep = strlength(strtrim(lines)) > 0 & ~startsWith(strtrim(lines), "#");
    lines = lines(keep);
    if isempty(lines)
        error('study_read_manifest:empty', 'Manifest %s is empty.', file);
    end
    sep = sprintf(delim);
    header = strtrim(split_line(lines(1), sep));
    if any(header == "")
        error('study_read_manifest:emptyHeader', 'Manifest %s has an empty column name.', file);
    end
    [~, iu] = unique(lower(header), 'stable');
    if numel(iu) ~= numel(header)
        error('study_read_manifest:duplicateColumn', 'Manifest %s has duplicate column names.', file);
    end
    body = strings(numel(lines) - 1, numel(header));
    for r = 2:numel(lines)
        cells = split_line(lines(r), sep);
        if numel(cells) > numel(header)
            if any(strtrim(cells(numel(header)+1:end)) ~= "")
                error('study_read_manifest:tooManyCells', ...
                    'Manifest %s line %d has more cells than columns.', file, r);
            end
            cells = cells(1:numel(header));
        end
        cells(end+1:numel(header)) = "";
        body(r-1, :) = strtrim(cells);
    end
    T = array2table(body, 'VariableNames', cellstr(header));
end

function cells = split_line(line, sep)
% Split one line at sep, honouring double quotes ("a, b" stays one cell;
% "" inside quotes is a literal quote).
    line = char(line);
    cells = strings(1, 0);
    cur = ''; inq = false; k = 1;
    while k <= numel(line)
        ch = line(k);
        if inq
            if ch == '"'
                if k < numel(line) && line(k+1) == '"'
                    cur(end+1) = '"'; k = k + 1; %#ok<AGROW>
                else
                    inq = false;
                end
            else
                cur(end+1) = ch; %#ok<AGROW>
            end
        elseif ch == '"'
            inq = true;
        elseif ch == sep
            cells(end+1) = string(cur); cur = ''; %#ok<AGROW>
        else
            cur(end+1) = ch; %#ok<AGROW>
        end
        k = k + 1;
    end
    cells(end+1) = string(cur);
end

function id = parse_subject_id(txt, r)
    tok = regexp(char(txt), '^(?:sub-)?(\d+)$', 'tokens', 'once');
    if isempty(tok)
        error('study_read_manifest:badSubject', ...
            'Manifest row %d: subject_id "%s" is not an integer (or sub-NNN).', r, txt);
    end
    id = str2double(tok{1});
end

function tf = parse_include(txt, r)
    switch lower(char(txt))
        case {'', '1', 'true', 'yes', 'y'}, tf = true;
        case {'0', 'false', 'no', 'n'},     tf = false;
        otherwise
            error('study_read_manifest:badInclude', ...
                'Manifest row %d: include must be 1/0, true/false or yes/no (got "%s").', r, txt);
    end
end

function v = parse_manifest_value(txt)
% YAML scalar or flow sequence, typed like values read from a config file.
    try
        v = yaml.load(txt, "ConvertToArray", true);
    catch
        v = txt;   % not valid YAML (e.g. an unquoted colon): keep the text
    end
    if isa(v, 'yaml.Null'); v = []; end
    if isstring(v) && isscalar(v); v = char(v); end
end

function tf = is_absolute(p)
    tf = startsWith(p, filesep) || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')) || startsWith(p, '~');
end
