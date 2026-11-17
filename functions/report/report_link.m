function href = report_link(target, report_dir)
% REPORT_LINK  href from a report to a file or folder
%
% A path relative to the report folder when both share the first path
% component (so the output folder can be copied together with what it links
% to), otherwise a file:// URI. Separators become '/'.
%
% Use as:
%   href = report_link(target, report_dir)
%
% Input:
%   target     - absolute path of the linked file or folder
%   report_dir - folder of the report
%
% Output:
%   href - char, for an href or src attribute (escape before inlining)
%
% See also: GENERATE_DEBUG_REPORT, GENERATE_SEGMENTATION_REPORT

arguments
    target     {mustBeTextScalar}
    report_dir {mustBeTextScalar}
end

    t = parts(char(target));
    r = parts(char(report_dir));
    n = 0;
    while n < min(numel(t), numel(r)) && strcmp(t{n+1}, r{n+1})
        n = n + 1;
    end
    if n == 0 || (n == 1 && isempty(t{1}))
        % no common root below the file-system root: absolute URI
        p = strrep(char(target), '\', '/');
        if ~startsWith(p, '/'); p = ['/' p]; end
        href = ['file://' p];
        return
    end
    up = repmat({'..'}, 1, numel(r) - n);
    href = strjoin([up, t(n+1:end)], '/');
    if isempty(href); href = '.'; end
end

function c = parts(p)
    p = strrep(p, '\', '/');
    while numel(p) > 1 && endsWith(p, '/')
        p = p(1:end-1);
    end
    c = strsplit(p, '/');
end
