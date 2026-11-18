function t = read_metrics_table(csv_path)
% READ_METRICS_TABLE  Read a per-run metrics CSV with current column names
%
% Reads sub-XXX_<medium><affix>.csv and renames columns of earlier
% versions to their current names, so every reader sees one convention:
% *_bowl_mm columns count from the bowl apex (trans_pos), *_ep_mm columns
% from the exit plane.
%
%   real_focal_distance_mm     -> real_focal_distance_bowl_mm   (bowl apex)
%   focal_centroid_3dB_dist_mm -> focal_centroid_3dB_bowl_mm    (bowl apex)
%
% It also drops the columns the HTML reports never show (MI_skull, MI_skin;
% REPORT_HIDDEN_COLUMNS). They stay in the CSV.
%
% Use as:
%   t = read_metrics_table(csv_path)
%
% Input:
%   csv_path - path to the metrics CSV written by acoustic_analysis
%
% Output:
%   t        - table, column names preserved as written, without the
%              report-hidden columns
%
% See also: ACOUSTIC_ANALYSIS, GENERATE_GROUP_REPORT

    t = readtable(csv_path, 'VariableNamingRule', 'preserve');

    renames = {'real_focal_distance_mm',     'real_focal_distance_bowl_mm'; ...
               'focal_centroid_3dB_dist_mm', 'focal_centroid_3dB_bowl_mm'};
    names = t.Properties.VariableNames;
    for r = 1:size(renames, 1)
        if ismember(renames{r, 1}, names) && ~ismember(renames{r, 2}, names)
            t = renamevars(t, renames{r, 1}, renames{r, 2});
        end
    end
    t = removevars(t, intersect(report_hidden_columns(), t.Properties.VariableNames));
end
