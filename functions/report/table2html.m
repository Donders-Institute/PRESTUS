function html = table2html(tbl, limits, color_columns)
% TABLE2HTML  Convert MATLAB table to HTML <table> with optional color-coded cells.
%
% Usage:
%   html = table2html(tbl)
%   html = table2html(tbl, limits)
%   html = table2html(tbl, limits, color_columns)
%
% Inputs:
%   tbl           - MATLAB table
%   limits        - (optional) struct of risk limits (from get_risk_limits)
%   color_columns - (optional) cell array of column names to color-code
arguments
    tbl
    limits        struct = struct()
    color_columns cell   = {}
end

    cols = tbl.Properties.VariableNames;

    html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    for c = 1:length(cols)
        html = [html sprintf('<th>%s</th>', html_utils.escape(cols{c}))];
    end
    html = [html '</tr></thead><tbody>'];

    for r = 1:height(tbl)
        html = [html '<tr>'];
        for c = 1:length(cols)
            col_name = cols{c};
            val = tbl{r, c};
            val_str = html_utils.format_cell(val);

            % Color coding for safety-relevant columns
            cell_class = '';
            if ismember(col_name, color_columns) && isfield(limits, col_name) && isnumeric(val) && isscalar(val)
                color = risk_color(val, limits.(col_name).limit);
                if ~strcmp(color, 'info')
                    cell_class = sprintf(' class="cell-%s"', color);
                end
            end

            html = [html sprintf('<td%s>%s</td>', cell_class, val_str)];
        end
        html = [html '</tr>'];
    end
    html = [html '</tbody></table></div>'];
end
