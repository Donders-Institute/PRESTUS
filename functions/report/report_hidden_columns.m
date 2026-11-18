function cols = report_hidden_columns()
% REPORT_HIDDEN_COLUMNS  Metrics-CSV columns that the HTML reports never show
%
% ACOUSTIC_ANALYSIS still writes these to the CSV, but READ_METRICS_TABLE
% drops them, so no report shows them (tables, dashboards, verdict): the
% ITRUSST mechanical criterion is the transcranial MI (MI_tc), and skull or
% skin MI only add noise to the safety view.
%
% Use as:
%   cols = report_hidden_columns()
%
% Output:
%   cols - cell array of CSV column names
%
% See also: READ_METRICS_TABLE, GET_RISK_LIMITS

    cols = {'MI_skull', 'MI_skin'};
end
