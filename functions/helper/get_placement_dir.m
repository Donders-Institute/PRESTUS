function placement_dir = get_placement_dir(parameters)
% GET_PLACEMENT_DIR  Per-subject root of the placement outputs
%
% Placement results depend on the subject's images, the target, the
% transducer geometry and the placement method, but not on the acoustic
% setup, so they have their own root and can be shared across simulation
% setups (see doc_outputs.md, "Where outputs live"):
%
%   <path.placement>/sub-XXX        when path.placement is set
%   <path.sim>/placement/sub-XXX    otherwise (path.sim taken as the parent
%                                   of io.dir_output when only that is set)
%
% One placement run lives in <placement_dir>/<mode>/<target_name>/
% (PLACEMENT_RUN_DIR).
%
% Use as:
%   placement_dir = get_placement_dir(parameters)
%
% Input:
%   parameters - PRESTUS parameter struct (subject_id, path.placement or path.sim)
%
% Output:
%   placement_dir - char path (not created)
%
% See also: PLACEMENT_RUN_DIR, GET_OUTPUT_DIR

arguments
    parameters (1,1) struct
end

    root = '';
    if isfield(parameters, 'path') && isfield(parameters.path, 'placement') && ...
            (ischar(parameters.path.placement) || isstring(parameters.path.placement)) && ...
            strlength(string(parameters.path.placement)) > 0
        root = char(parameters.path.placement);
    elseif isfield(parameters, 'path') && isfield(parameters.path, 'sim') && ...
            (ischar(parameters.path.sim) || isstring(parameters.path.sim)) && ...
            strlength(string(parameters.path.sim)) > 0
        root = fullfile(char(parameters.path.sim), 'placement');
    elseif isfield(parameters, 'io') && isfield(parameters.io, 'dir_output') && ...
            ~isempty(parameters.io.dir_output)
        % io.dir_output is <path.sim>/sub-XXX (GET_OUTPUT_DIR)
        root = fullfile(fileparts(char(parameters.io.dir_output)), 'placement');
    end
    if isempty(root)
        error('get_placement_dir:missingPath', ...
            'Set path.placement or path.sim to locate the placement outputs.');
    end

    placement_dir = fullfile(root, sprintf('sub-%03d', parameters.subject_id));
end
