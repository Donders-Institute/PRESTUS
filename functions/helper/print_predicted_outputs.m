function plan = print_predicted_outputs(parameters)
% PRINT_PREDICTED_OUTPUTS  Print PREDICT_OUTPUTS' plan to the console.
%
% Call this on a resolved config (after load_parameters) before handing it to
% PRESTUS_PIPELINE_START, to see from the command line which pipeline mode it
% will dispatch to and which files it will write.
%
% Use as:
%   print_predicted_outputs(parameters)
%   plan = print_predicted_outputs(parameters)   % also returns the plan struct
%
% See also: PREDICT_OUTPUTS, PRESTUS_PIPELINE_START

arguments
    parameters (1,1) struct
end

plan = predict_outputs(parameters);

mode_label = plan.mode;
if ~strcmp(mode_label, 'single')
    mode_label = strrep(mode_label, '_', '-');
end
fprintf('PRESTUS output prediction — sub-%03d (%s mode)\n', parameters.subject_id, mode_label);
for w = 1:numel(plan.warnings)
    fprintf('  ! %s\n', plan.warnings{w});
end
if ~isempty(plan.jobs)
    fprintf('  Deployment:\n');
    for j = 1:numel(plan.jobs)
        fprintf('    - %s\n', plan.jobs{j});
    end
end
fprintf('  %s/\n', plan.output_dir);
width = 10;
for a = 1:numel(plan.artifacts)
    width = max(width, length(plan.artifacts(a).name) + 2);
end
for a = 1:numel(plan.artifacts)
    fprintf('      %-*s%s\n', width, plan.artifacts(a).name, plan.artifacts(a).note);
end
for d = 1:numel(plan.common_dirs)
    fprintf('      %-*s(regenerable intermediates)\n', width, plan.common_dirs{d});
end
end
