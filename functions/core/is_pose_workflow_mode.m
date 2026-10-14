function tf = is_pose_workflow_mode(parameters)
% IS_POSE_WORKFLOW_MODE  True when time_reversal.pose_workflow.enabled is set
%
% Use as:
%   tf = is_pose_workflow_mode(parameters)
%
% See also: POSE_WORKFLOW_PIPELINE

    tf = isfield(parameters, 'time_reversal') && isfield(parameters.time_reversal, 'pose_workflow') && ...
        isfield(parameters.time_reversal.pose_workflow, 'enabled') && ...
        logical(parameters.time_reversal.pose_workflow.enabled);
end
