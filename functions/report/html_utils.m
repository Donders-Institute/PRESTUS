classdef html_utils
% html_utils  HTML utility methods for PRESTUS report generation.
%
%   html_utils.escape(str)                        — escape HTML entities
%   html_utils.base64(filepath)                   — encode file as base64
%   html_utils.embed_image(filepath,alt,caption)  — <figure> with embedded image
%   html_utils.collapsible(title,html,open,id)    — <details> collapsible section
%   html_utils.section_error(name,ME)             — error fallback section
%   html_utils.lightbox()                         — lightbox overlay markup + JS
%   html_utils.format_cell(val)                   — format table cell value

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
        % lightbox  Returns the lightbox overlay markup and JavaScript.
            html = ['<div id="lightbox" class="lb-overlay" style="display:none;" onclick="closeLightbox()">' ...
                    '<span class="lb-close" onclick="closeLightbox()">&times;</span>' ...
                    '<img id="lb-img" src="" alt="Enlarged image" onclick="event.stopPropagation()">' ...
                    '</div>' ...
                    '<script>' ...
                    'document.querySelectorAll(".image-grid figure img,.image-grid--3col figure img").forEach(function(img){' ...
                    'img.addEventListener("click",function(){' ...
                    'var lb=document.getElementById("lightbox");' ...
                    'document.getElementById("lb-img").src=this.src;lb.style.display="flex";});});' ...
                    'function closeLightbox(){document.getElementById("lightbox").style.display="none";}' ...
                    'document.addEventListener("keydown",function(e){if(e.key==="Escape")closeLightbox();});' ...
                    '</script>'];
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
