function debug_register(parameters, stage, file, caption, kind, group)
% DEBUG_REGISTER  Record a debug output in the run's debug manifest
%
% Appends one JSON line to <io.dir_debug>/manifest.jsonl, which
% GENERATE_DEBUG_REPORT reads to show every debug output with its caption,
% stage and the function that wrote it. Files inside io.dir_debug are stored
% relative to it, so the output folder can be moved. Never errors: a failed
% write is a warning.
%
% Use as:
%   debug_register(parameters, stage, file, caption)
%   debug_register(parameters, stage, file, caption, kind)
%   debug_register(parameters, stage, file, caption, kind, group)
%
% Input:
%   parameters - PRESTUS config (io.dir_debug)
%   stage      - debug stage (DEBUG_STAGES)
%   file       - path of the written file (an extension-less NIfTI path from
%                niftiwrite is resolved to the file actually written)
%   caption    - one-line description shown in the debug report
%   kind       - 'figure' | 'volume' | 'other' (default: from the extension)
%   group      - name shared by several files of one figure (its views,
%                DEBUG_SAVE_FIGURE); the debug report shows them side by side
%
% See also: DEBUG_SAVE_FIGURE, DEBUG_ON, GENERATE_DEBUG_REPORT

arguments
    parameters (1,1) struct
    stage      {mustBeTextScalar}
    file       {mustBeTextScalar}
    caption    {mustBeTextScalar} = ''
    kind       {mustBeTextScalar} = ''
    group      {mustBeTextScalar} = ''
end

    try
        if ~isfield(parameters, 'io') || ~isfield(parameters.io, 'dir_debug') || isempty(parameters.io.dir_debug)
            return
        end
        file = resolve_written(char(file));
        if strlength(kind) == 0
            kind = kind_from_ext(file);
        end
        root = char(parameters.io.dir_debug);
        if ~isfolder(root); mkdir(root); end

        rel = file;
        prefix = [root filesep];
        if strncmp(file, prefix, numel(prefix))
            rel = file(numel(prefix)+1:end);
            rel = strrep(rel, filesep, '/');
        end

        st = dbstack(1);
        src = '';
        if ~isempty(st)
            src = st(1).name;
            if strcmp(src, 'debug_save_figure') && numel(st) > 1
                src = st(2).name;
            end
        end
        affix = '';
        if isfield(parameters.io, 'output_affix'); affix = char(parameters.io.output_affix); end

        rec = struct('stage', char(stage), 'file', rel, 'caption', char(caption), ...
            'kind', char(kind), 'group', char(group), 'source', src, 'affix', affix, ...
            'time', char(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss')));
        fid = fopen(fullfile(root, 'manifest.jsonl'), 'a', 'n', 'UTF-8');
        if fid == -1
            warning('PRESTUS:debug:manifest', 'Cannot write the debug manifest in %s.', root);
            return
        end
        fprintf(fid, '%s\n', jsonencode(rec));
        fclose(fid);
    catch ME
        warning('PRESTUS:debug:manifest', 'Could not register debug output %s: %s', char(file), ME.message);
    end
end

function file = resolve_written(file)
% niftiwrite/saveas may add an extension to the requested name
    if isfile(file), return; end
    for ext = {'.gz', '.nii.gz', '.nii', '.png'}
        if isfile([file ext{1}])
            file = [file ext{1}];
            return
        end
    end
end

function kind = kind_from_ext(file)
    f = lower(file);
    if endsWith(f, {'.png', '.jpg', '.jpeg', '.svg'})
        kind = 'figure';
    elseif endsWith(f, {'.nii', '.nii.gz'})
        kind = 'volume';
    else
        kind = 'other';
    end
end
