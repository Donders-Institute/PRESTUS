function [selected, known] = debug_stages(parameters)
% DEBUG_STAGES  Stages selected by simulation.debug_stages
%
% Each stage writes its debug outputs to its own subfolder of io.dir_debug
% (io.dir_debug_<stage>).
%
% Use as:
%   [selected, known] = debug_stages(parameters)
%
% Input:
%   parameters - PRESTUS config
%
% Output:
%   selected - cellstr of selected stages ('all' or unset = all known stages)
%   known    - cellstr of all stages: preproc, medium, source, babelbrain
%
% An unknown stage name is an error, so a typo does not silently switch
% debug outputs off.
%
% See also: DEBUG_ON, PATH_LOG_SETUP

arguments
    parameters (1,1) struct
end

    known = {'preproc', 'medium', 'source', 'babelbrain'};
    selected = known;
    if ~isfield(parameters, 'simulation') || ~isfield(parameters.simulation, 'debug_stages')
        return
    end
    v = parameters.simulation.debug_stages;
    if isempty(v) || (ischar(v) || isstring(v)) && any(strcmpi(cellstr(v), 'all'))
        return
    end
    if ischar(v) || isstring(v)
        v = cellstr(v);
    end
    if iscell(v) && any(strcmpi(v, 'all'))
        return
    end
    if ~iscellstr(v)
        error('PRESTUS:debug:stages', ...
            'simulation.debug_stages must be ''all'' or a list of stages (%s).', strjoin(known, ', '));
    end
    v = lower(v(:)');
    bad = setdiff(v, known);
    if ~isempty(bad)
        error('PRESTUS:debug:unknownStage', ...
            'simulation.debug_stages: unknown stage(s) %s; known stages: %s.', ...
            strjoin(bad, ', '), strjoin(known, ', '));
    end
    selected = v;
end
