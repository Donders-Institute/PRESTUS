%% EXAMPLE_REPORT_ERNIE  Head simulation reported against the ITRUSST consensus
% The tutorial's simulation (tutorial_ernie.yaml) with heating, the free-water run, targeting
% analysis, and CEM43 in the ISO form (example_report_ernie.yaml). After the run, the metrics
% table is printed against the ITRUSST non-significant-risk limits used by the HTML report
% (get_risk_limits.m). PRESTUS reports the comparison; it makes no safety recommendation.
%
% Before running: set the paths in tutorial_ernie.yaml and run its segmentation step
% (section 1 of tutorial_ernie.m). For your own image, change path.* there; nothing here.
% Not yet run end to end.

prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'config')));
addpath(genpath(fullfile(prestus_root, 'external')));

cfg_dir = fullfile(prestus_root, 'examples');
parameters = load_parameters('tutorial_ernie.yaml', cfg_dir, ...
                             'example_report_ernie.yaml', cfg_dir, {'subject_id', 1});

%% 1. Run (placement, acoustic, thermal, free-water run, report)
prestus_pipeline_start(parameters);

%% 2. Metrics against the ITRUSST limits
% Same colour rule as the HTML report: within (< 50 % of a limit), near (50-100 %), above.
csv_file = fullfile(get_output_dir(parameters), sprintf('sub-%03d_%s%s.csv', ...
    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
metrics = readtable(csv_file);
limits  = get_risk_limits(true);
iso     = struct('CEM43_brain', 'CEM43iso_brain', 'CEM43_skull', 'CEM43iso_skull', ...
                 'CEM43_skin', 'CEM43iso_skin');   % thermal.cem43_iso: 1

fprintf('\n%-28s %10s %10s  %s\n', 'Metric', 'Value', 'Limit', 'Status');
names = fieldnames(limits);
for i = 1:numel(names)
    col = names{i};
    if isfield(iso, col), col = iso.(col); end
    if ~ismember(col, metrics.Properties.VariableNames), continue, end
    value = metrics.(col)(1);
    lim   = limits.(names{i}).limit;
    if isinf(lim)
        status = 'no limit';
    elseif value < 0.5 * lim
        status = 'within';
    elseif value <= lim
        status = 'near';
    else
        status = 'above';
    end
    fprintf('%-28s %10.3g %10.3g  %s\n', limits.(names{i}).label, value, lim, status);
end
fprintf('\nThe thermal criteria are alternatives (temperature rise, maximum temperature, or dose);\n');
fprintf('see doc/doc_results.md. Full report: %s\n', get_output_dir(parameters));
