classdef html_utils
% HTML_UTILS  Static HTML utility methods for PRESTUS report generation
%
% All methods are static; call as html_utils.method_name(...).
% Images are embedded as base64 data URIs so reports are self-contained.
%
% Methods:
%   html_utils.escape(str)                       — escape <, >, &, " for safe HTML
%   html_utils.base64(filepath)                  — encode a file as a base64 string
%   html_utils.embed_image(filepath, alt, cap)   — <figure> with inline base64 image
%   html_utils.collapsible(title, html, open, id)— <details>/<summary> collapsible block
%   html_utils.section_error(name, ME)           — error fallback <section> with message
%   html_utils.lightbox()                        — lightbox overlay markup + JS
%   html_utils.format_cell(val)                  — format a value for a table cell
%
% See also: GENERATE_SIMULATION_REPORT, GENERATE_UNCERTAINTY_REPORT, CSS_STYLES_BASE

    methods (Static)

        function str = escape(str)
        % escape  Escape <, >, &, and " for safe HTML embedding.
            if ~ischar(str), str = char(string(str)); end
            str = strrep(str, '&', '&amp;');
            str = strrep(str, '<', '&lt;');
            str = strrep(str, '>', '&gt;');
            str = strrep(str, '"', '&quot;');
        end

        function b64 = base64(filepath)
        % base64  Read image file and return base64 string.
            b64 = '';
            if ~isfile(filepath), return; end
            fid = fopen(filepath, 'r');
            if fid == -1, return; end
            raw = fread(fid, '*uint8');
            fclose(fid);
            b64 = matlab.net.base64encode(raw);
        end

        function html = embed_image(filepath, alt_text, caption)
        % embed_image  Returns <figure> with base64-encoded <img>, or empty string if file missing.
            html = '';
            if ~isfile(filepath), return; end

            b64 = html_utils.base64(filepath);
            if isempty(b64), return; end

            % Detect MIME type
            [~, ~, ext] = fileparts(filepath);
            switch lower(ext)
                case '.png',           mime = 'image/png';
                case {'.jpg','.jpeg'}, mime = 'image/jpeg';
                otherwise,             mime = 'image/png';
            end

            html = '<figure>';
            html = [html sprintf('<img src="data:%s;base64,%s" alt="%s">', ...
                mime, b64, html_utils.escape(alt_text))];
            html = [html sprintf('<figcaption>%s</figcaption>', html_utils.escape(caption))];
            html = [html '</figure>'];
        end

        function html = embed_image_dims(img_dir, base, suffix, alt_text, caption)
        % embed_image_dims  Embed a plot that may be saved either as a single
        % no-dimension file (2D grids: base+suffix.png) or as per-slice files
        % with an _x/_y/_z infix between the metric and the suffix (3D grids:
        % base+_x/_y/_z+suffix.png). Probes all four and embeds every file
        % found, annotating the caption with the slice dimension. Returns ''
        % when none exist. This makes the report robust to the per-dimension
        % naming used by thermal_analysis (maxT) and acoustic_analysis (intensity).
            html = '';
            variants = {'', '_x', '_y', '_z'};
            for k = 1:numel(variants)
                dim = variants{k};
                fpath = fullfile(img_dir, [base dim suffix '.png']);
                if isfile(fpath)
                    if isempty(dim)
                        cap = caption;
                    else
                        cap = sprintf('%s (%s-slice)', caption, dim(2:end));
                    end
                    html = [html html_utils.embed_image(fpath, alt_text, cap)];
                end
            end
        end

        function html = report_logo(logo_height)
        % report_logo  Inline base64 <img> of the vendored PRESTUS logo, or ''
        % if the asset is missing. Resolves the path relative to this file so
        % it works regardless of the working directory (e.g. on the HPC).
            if nargin < 1 || isempty(logo_height), logo_height = 30; end
            html = '';
            here = fileparts(mfilename('fullpath'));
            fpath = fullfile(here, 'assets', 'logo_PRESTUS.png');
            if ~isfile(fpath), return; end
            b64 = html_utils.base64(fpath);
            if isempty(b64), return; end
            html = sprintf(['<img class="brand-logo" alt="PRESTUS" ' ...
                'src="data:image/png;base64,%s" style="height:%dpx;width:auto;display:block">'], ...
                b64, round(logo_height));
        end

        function html = collapsible(title, content_html, is_open, section_id)
        % collapsible  Wraps content in a <details>/<summary> collapsible section.
        % is_open: true to default open, false to default collapsed.
            if is_open
                html = sprintf('<details class="report-section" id="%s" open>', section_id);
            else
                html = sprintf('<details class="report-section" id="%s">', section_id);
            end
            html = [html sprintf('<summary><h2>%s</h2></summary>', html_utils.escape(title))];
            html = [html '<div class="section-content">'];
            html = [html content_html];
            html = [html '</div></details>'];
        end

        function html = section_error(section_name, ME)
        % section_error  Returns an error-styled HTML section block.
            html = sprintf(['<section class="report-section error-section">' ...
                '<h2>%s</h2><p class="error-notice">Section failed: %s</p></section>'], ...
                html_utils.escape(section_name), html_utils.escape(ME.message));
        end

        function html = lightbox()
        % lightbox  Returns the lightbox overlay markup. The behaviour (open,
        % zoom, side-by-side compare, Esc-to-close) lives in
        % report_scripts.common() so there is a single source of JS.
            html = ['<div id="lightbox" class="lb-overlay" style="display:none">' ...
                    '<button class="lb-close" type="button" aria-label="Close">&times;</button>' ...
                    '<div class="lb-stage" id="lb-stage"></div>' ...
                    '<div class="lb-ctrls">' ...
                    '<button type="button" id="lb-zoom">Zoom +</button>' ...
                    '<button type="button" id="lb-reset">Reset</button>' ...
                    '<span style="color:#9aa7b4;font-size:.8rem">Tip: mark two figures &ldquo;compare&rdquo;, then click either</span>' ...
                    '</div></div>'];
        end

        function str = format_cell(val)
        % format_cell  Format a table cell value for HTML display.
            if isnumeric(val)
                if isscalar(val)
                    if isnan(val)
                        str = 'N/A';
                    elseif val == round(val) && abs(val) < 1e6
                        str = sprintf('%d', val);
                    else
                        str = sprintf('%.4g', val);
                    end
                else
                    str = ['[' strtrim(sprintf('%.4g ', val)) ']'];
                end
            elseif ischar(val) || isstring(val)
                str = html_utils.escape(char(val));
            elseif iscell(val)
                str = html_utils.escape(char(val{1}));
            else
                str = '—';
            end
        end

    end
end
