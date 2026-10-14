function run_dir = placement_run_dir(parameters, mode)
% PLACEMENT_RUN_DIR  Folder of one placement run (one mode, one target)
%
% <get_placement_dir>/<mode>/<target_name>, or placement.<mode>.work_dir when
% set. Everything the run produces lives here: the tool's inputs (masks,
% configs, runner scripts), its raw output, metrics, neuronavigation exports
% and the run record (placement_run.mat). The SimNIBS folders are only read.
%
% Modes with a run folder: heuristic, plantus, tus_entry. manual and mni are
% resolved from the config alone, localite from the recorded/ingested
% markers, so they have none (the simulation keeps its own record of the
% resolved positions, WRITE_PLACEMENT_POSITIONS).
%
% Use as:
%   run_dir = placement_run_dir(parameters, mode)
%
% Input:
%   parameters - PRESTUS parameter struct
%   mode       - 'heuristic' | 'plantus' | 'tus_entry'
%
% Output:
%   run_dir - char path (not created)
%
% See also: GET_PLACEMENT_DIR, PREPROC_TRANSDUCER_PLACEMENT

arguments
    parameters (1,1) struct
    mode       {mustBeTextScalar}
end

    mode = lower(char(mode));
    if ~any(strcmp(mode, {'heuristic', 'plantus', 'tus_entry'}))
        error('PRESTUS:placement:noRunDir', ...
            'placement.mode ''%s'' has no placement run folder.', mode);
    end

    s = struct();
    if isfield(parameters, 'placement') && isfield(parameters.placement, mode) && ...
            isstruct(parameters.placement.(mode))
        s = parameters.placement.(mode);
    end
    if isfield(s, 'work_dir') && ~isempty(s.work_dir)
        run_dir = char(s.work_dir);
        return
    end

    target_name = 'target';
    if isfield(s, 'target_name') && ~isempty(s.target_name)
        target_name = char(s.target_name);
    end
    run_dir = fullfile(get_placement_dir(parameters), mode, target_name);
end
