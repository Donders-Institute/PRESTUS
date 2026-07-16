function html = build_config_dump(parameters)
% BUILD_CONFIG_DUMP  Searchable HTML block of the full parameters struct.
%
% Serialises the parameters actually used by prestus_pipeline to YAML so the
% complete acoustic/thermal configuration is visible in the report. Large
% numeric arrays, cached masks and handles are pruned to keep the dump
% readable and the file small. Falls back to a plain text serialisation when
% the YAML toolbox (SnakeYAML) is unavailable.
%
% Use as:
%   html = build_config_dump(parameters)   % returns an HTML fragment
%
% See also: GENERATE_SIMULATION_REPORT, BUILD_CONFIG_SUMMARY

    try
        pruned = local_prune(parameters, 0);
        try
            txt = char(yaml.dump(pruned, "block"));
        catch
            txt = local_struct2text(pruned, 0);
        end
    catch ME
        html = sprintf('<p class="placeholder">Full configuration unavailable: %s</p>', ...
            html_utils.escape(ME.message));
        return
    end

    safe = html_utils.escape(txt);
    html = ['<div class="config-dump-bar">' ...
        '<input id="config-dump-search" type="search" placeholder="Filter the full configuration…" autocomplete="off">' ...
        '<span style="font-size:.74rem;color:var(--ink-4)">parameters as used by prestus_pipeline</span>' ...
        '</div>'];
    html = [html '<pre class="config-dump" id="config-dump">' safe '</pre>'];
    html = [html local_dump_script()];
end

% ------------------------------------------------------------------------
function out = local_prune(s, depth)
% Recursively copy a struct, replacing bulky / non-serialisable content with
% short placeholders so the YAML dump stays small and readable.
    MAX_DEPTH = 7;
    MAX_NUMEL = 64;
    if depth > MAX_DEPTH
        out = '[... pruned: max depth ...]';
        return
    end
    if isstruct(s)
        if numel(s) ~= 1
            out = struct('struct_array', sprintf('[%s struct]', local_size_str(size(s))));
            if ~isempty(s)
                try, out.first_element = local_prune(s(1), depth + 1); catch, end
            end
            return
        end
        out = struct();
        f = fieldnames(s);
        for i = 1:numel(f)
            out.(f{i}) = local_prune(s.(f{i}), depth + 1);
        end
    elseif isnumeric(s) || islogical(s)
        if numel(s) > MAX_NUMEL
            out = sprintf('[%s %s pruned]', class(s), local_size_str(size(s)));
        elseif isempty(s)
            out = '[]';
        else
            out = s;
        end
    elseif ischar(s)
        out = s;
    elseif isstring(s)
        if numel(s) > MAX_NUMEL
            out = sprintf('[string %s pruned]', local_size_str(size(s)));
        else
            out = char(strjoin(string(s(:))', ' | '));
        end
    elseif iscell(s)
        if numel(s) > MAX_NUMEL
            out = sprintf('{%s cell pruned}', local_size_str(size(s)));
        else
            out = cell(size(s));
            for i = 1:numel(s)
                out{i} = local_prune(s{i}, depth + 1);
            end
        end
    elseif isa(s, 'function_handle')
        out = ['@' char(func2str(s))];
    else
        out = sprintf('[%s object]', class(s));
    end
end

% ------------------------------------------------------------------------
function str = local_size_str(sz)
    str = strjoin(arrayfun(@(x) sprintf('%d', x), sz, 'UniformOutput', false), 'x');
end

% ------------------------------------------------------------------------
function txt = local_struct2text(s, indent)
% Minimal YAML-ish fallback serialiser (used when yaml.dump is unavailable).
    pad = repmat('  ', 1, indent);
    txt = '';
    if isstruct(s) && isscalar(s)
        f = fieldnames(s);
        for i = 1:numel(f)
            v = s.(f{i});
            if (isstruct(v) && isscalar(v)) || (iscell(v) && ~isempty(v))
                txt = [txt sprintf('%s%s:\n', pad, f{i}) local_struct2text(v, indent + 1)];
            else
                txt = [txt sprintf('%s%s: %s\n', pad, f{i}, local_val2str(v))];
            end
        end
    elseif iscell(s)
        for i = 1:numel(s)
            txt = [txt sprintf('%s- %s\n', pad, local_val2str(s{i}))];
        end
    else
        txt = sprintf('%s%s\n', pad, local_val2str(s));
    end
end

% ------------------------------------------------------------------------
function str = local_val2str(v)
    if ischar(v)
        str = v;
    elseif isstring(v)
        str = char(v);
    elseif islogical(v)
        str = mat2str(v);
    elseif isnumeric(v)
        if isscalar(v), str = num2str(v); else, str = mat2str(v); end
    elseif isstruct(v)
        str = sprintf('[%s struct]', local_size_str(size(v)));
    elseif iscell(v)
        str = sprintf('{%s cell}', local_size_str(size(v)));
    else
        str = class(v);
    end
end

% ------------------------------------------------------------------------
function html = local_dump_script()
% Client-side filter that keeps only matching lines of the config dump.
    js = ['(function(){var cd=document.getElementById("config-dump"),' ...
        'cs=document.getElementById("config-dump-search");if(!cd||!cs)return;var raw=cd.textContent;' ...
        'cs.addEventListener("input",function(e){var q=e.target.value.trim();if(!q){cd.textContent=raw;return;}' ...
        'var rx=new RegExp("("+q.replace(/[.*+?^${}()|[\]\\]/g,"\\$&")+")","gi");' ...
        'var lines=raw.split("\n").filter(function(l){return rx.test(l);});' ...
        'cd.innerHTML=lines.length?lines.join("\n").replace(rx,"<mark>$1</mark>"):' ...
        '"<span style=\"color:var(--ink-4)\">no matching keys</span>";});})();'];
    html = ['<script>' js '</script>'];
end
