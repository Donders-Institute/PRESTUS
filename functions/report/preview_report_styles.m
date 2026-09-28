function out = preview_report_styles(target, opts)
% PREVIEW_REPORT_STYLES  Re-render an existing simulation report in several
% visual styles, without re-running any simulation ("report test mode")
%
% Loads the parameter snapshot that the pipeline saved in <subject>/cache,
% re-points its io paths at the given directory (so outputs copied from the
% HPC or mounted elsewhere work), and calls generate_simulation_report once
% per style. All reports are written to a sandbox folder, never over the
% original report, together with an index.html that shows two styles side by
% side. Because the report is rebuilt from the CSV/PNG/log files on disk,
% this is also the quickest way to test edits to the report code itself.
%
% Use as:
%   preview_report_styles('/path/to/tussim/<run>/sub-602')
%   preview_report_styles(dir, 'affix', '_kwave_cpp_pfov')     % pick one run in a sweep folder
%   preview_report_styles(dir, 'styles', {'instrument','classic'}, 'open', true)
%   preview_report_styles('/path/to/cache/sub-602_layered_parameters_260925_1233.mat')
%
% Input:
%   target        - subject output directory, or a *_parameters_*.mat snapshot
%   opts.styles   - cell array of style names (default: report_styles.names())
%   opts.affix    - output_affix of the snapshot to use when a directory holds
%                   several runs; default [] = newest snapshot
%   opts.out_dir  - destination (default: <subject dir>/report_preview)
%   opts.switcher - also write one report with all styles + live dropdown (default true)
%   opts.open     - open index.html in the system browser (default false)
%
% Output:
%   out - struct with fields out_dir, index, files (cellstr), failed (cellstr)
%
% See also: REPORT_STYLES, GENERATE_SIMULATION_REPORT

arguments
    target {mustBeTextScalar}
    opts.styles cell = report_styles.names()
    opts.affix = []
    opts.out_dir {mustBeTextScalar} = ''
    opts.switcher (1,1) logical = true
    opts.open (1,1) logical = false
end

target = char(target);
[parameters, sub_dir, snap] = load_snapshot(target, opts.affix);
parameters = repoint_io(parameters, sub_dir);
fprintf('preview_report_styles: using snapshot %s\n', snap);

out_dir = char(opts.out_dir);
if isempty(out_dir); out_dir = fullfile(sub_dir, 'report_preview'); end
if ~isfolder(out_dir); mkdir(out_dir); end

% Render into a scratch dir, then rename: the report filename is derived
% from subject/medium/affix, and the affix also locates the input images.
scratch = tempname; mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch, 's'));
parameters.io.dir_reports = scratch;

variants = struct('name', {}, 'switcher', {});
for k = 1:numel(opts.styles)
    variants(end+1) = struct('name', char(opts.styles{k}), 'switcher', false); %#ok<AGROW>
end
if opts.switcher
    variants(end+1) = struct('name', 'instrument', 'switcher', true);
end

out = struct('out_dir', out_dir, 'index', '', 'files', {{}}, 'failed', {{}});
labels = {};
for k = 1:numel(variants)
    v = variants(k);
    p = parameters;
    p.modules.report_style = v.name;
    p.modules.report_style_switcher = v.switcher;
    label = v.name;
    if v.switcher; label = 'switcher'; end
    rp = generate_simulation_report(p);
    if isempty(rp) || ~isfile(rp)
        out.failed{end+1} = label;
        continue
    end
    [~, base] = fileparts(rp);
    dest = fullfile(out_dir, sprintf('%s_style-%s.html', base, label));
    movefile(rp, dest, 'f');
    out.files{end+1} = dest;
    labels{end+1} = label; %#ok<AGROW>
end

out.index = write_index(out_dir, out.files, labels, snap);
fprintf('preview_report_styles: %d report(s) in %s\n  index: %s\n', ...
    numel(out.files), out_dir, out.index);
if ~isempty(out.failed)
    warning('preview_report_styles:failed', 'Report generation failed for: %s', ...
        strjoin(out.failed, ', '));
end
if opts.open && ~isempty(out.files)
    web(out.index, '-browser');
end
end

%% ------------------------------------------------------------------------

function [parameters, sub_dir, snap] = load_snapshot(target, affix)
% Locate + load the parameters snapshot written by path_log_setup.
    if isfile(target)
        snap = target;
        sub_dir = fileparts(fileparts(snap));   % <sub>/cache/<file>.mat
    elseif isfolder(target)
        sub_dir = target;
        % Preferred: final state saved by prestus_pipeline right before the
        % report (exact re-render). Else: start-of-run snapshot, which lacks
        % a few runtime-derived fields (e.g. resolved grid size, k-Wave version).
        rep = dir(fullfile(sub_dir, 'cache', '*_report_parameters*.mat'));
        if ~isempty(rep)
            rep_affix = regexprep({rep.name}, '^.*_report_parameters(.*)\.mat$', '$1');
            if ~isempty(affix); rep = rep(strcmp(rep_affix, char(affix))); end
        end
        if ~isempty(rep)
            [~, newest] = max([rep.datenum]);
            snap = fullfile(sub_dir, 'cache', rep(newest).name);
            S = load(snap, 'parameters');
            parameters = S.parameters;
            return
        end
        files = dir(fullfile(sub_dir, 'cache', '*_parameters_*.mat'));
        if isempty(files)
            error('preview_report_styles:noSnapshot', ...
                'No *_parameters_*.mat snapshot in %s', fullfile(sub_dir, 'cache'));
        end
        names = {files.name};
        tok = regexp(names, '_parameters_(\d{6}_\d{4})(.*)\.mat$', 'tokens', 'once');
        keep = ~cellfun(@isempty, tok);
        names = names(keep); tok = tok(keep);
        stamps = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
        affixes = cellfun(@(t) t{2}, tok, 'UniformOutput', false);
        if ~isempty(affix)
            m = strcmp(affixes, char(affix));
            if ~any(m)
                error('preview_report_styles:noAffix', ...
                    'No snapshot with affix "%s". Available: %s', char(affix), ...
                    strjoin(unique(affixes), ', '));
            end
            names = names(m); stamps = stamps(m);
        end
        [~, newest] = max(datenum(stamps, 'yymmdd_HHMM')); %#ok<DATNM>
        snap = fullfile(sub_dir, 'cache', names{newest});
    else
        error('preview_report_styles:notFound', '%s is neither a folder nor a file.', target);
    end
    S = load(snap, 'parameters');
    parameters = S.parameters;
end

function p = repoint_io(p, sub_dir)
% Swap the saved output-dir prefix for sub_dir in every io path field, so
% snapshots written on the HPC (or before a folder was moved) resolve here.
    if ~isfield(p, 'io'); p.io = struct(); end
    old = '';
    if isfield(p.io, 'dir_output'); old = char(p.io.dir_output); end
    p.io.dir_output = sub_dir;
    f = fieldnames(p.io);
    for i = 1:numel(f)
        v = p.io.(f{i});
        if ~(ischar(v) || (isstring(v) && isscalar(v))) || strcmp(f{i}, 'dir_output'); continue; end
        v = char(v);
        if ~isempty(old) && startsWith(v, old)
            p.io.(f{i}) = [sub_dir v(numel(old)+1:end)];
        end
    end
    % Fill anything still missing (older snapshots)
    p = resolve_io_dirs(p);
end

function index = write_index(out_dir, files, labels, snap)
% Side-by-side comparison page: two iframes, each with its own style picker.
    index = fullfile(out_dir, 'index.html');
    opts_html = '';
    for k = 1:numel(files)
        [~, n, e] = fileparts(files{k});
        opts_html = [opts_html sprintf('<option value="%s">%s</option>', ...
            html_utils.escape([n e]), html_utils.escape(labels{k}))]; %#ok<AGROW>
    end
    right = min(2, numel(files));
    html = ['<!DOCTYPE html><html><head><meta charset="UTF-8"><title>Report styles</title><style>' ...
        'body{margin:0;font-family:system-ui,sans-serif;background:#e9edf1}' ...
        'header{padding:8px 14px;font-size:.85rem;color:#3d4956;display:flex;gap:16px;align-items:center;flex-wrap:wrap}' ...
        'header code{font-size:.75rem;color:#6b7886}' ...
        '.panes{display:grid;grid-template-columns:1fr 1fr;gap:6px;height:calc(100vh - 44px);padding:0 6px 6px}' ...
        '.pane{display:flex;flex-direction:column;min-width:0}.pane select{margin-bottom:4px;font-size:.85rem}' ...
        'iframe{flex:1;border:1px solid #cdd6df;border-radius:6px;background:#fff;width:100%}' ...
        'label{display:flex;gap:6px;align-items:center}' ...
        '</style></head><body><header><strong>PRESTUS report styles</strong>' ...
        sprintf('<code>%s</code>', html_utils.escape(snap)) ...
        '<label><input type="checkbox" id="sync" checked> sync scroll</label></header>' ...
        '<div class="panes">' ...
        sprintf('<div class="pane"><select data-pane="a">%s</select><iframe id="a"></iframe></div>', opts_html) ...
        sprintf('<div class="pane"><select data-pane="b">%s</select><iframe id="b"></iframe></div>', opts_html) ...
        '</div><script>(function(){var sel=document.querySelectorAll("select[data-pane]");' ...
        sprintf('sel[0].selectedIndex=0;sel[1].selectedIndex=%d;', right - 1) ...
        'function load(s){document.getElementById(s.dataset.pane).src=s.value;}' ...
        'sel.forEach(function(s){load(s);s.addEventListener("change",function(){load(s);});});' ...
        'var lock=false;function hook(id,other){var f=document.getElementById(id);f.addEventListener("load",function(){' ...
        'try{f.contentWindow.addEventListener("scroll",function(){if(!document.getElementById("sync").checked||lock)return;' ...
        'var w=f.contentWindow,o=document.getElementById(other).contentWindow,d=w.document.documentElement,od=o.document.documentElement;' ...
        'lock=true;o.scrollTo(0,w.scrollY/Math.max(1,d.scrollHeight-w.innerHeight)*(od.scrollHeight-o.innerHeight));' ...
        'setTimeout(function(){lock=false;},30);});}catch(e){}});}' ...
        'hook("a","b");hook("b","a");})();</script></body></html>'];
    fid = fopen(index, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s', html);
    fclose(fid);
end
