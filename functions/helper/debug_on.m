function tf = debug_on(parameters, stage, kind)
% DEBUG_ON  Whether a stage writes a kind of debug output
%
% simulation.debug is the master switch. With it on, the finer settings
% select what is written:
%   simulation.debug_stages    - 'all' (default) or a list of stages
%                                (DEBUG_STAGES lists them)
%   simulation.debug_figures   - diagnostic figures (PNG; default 1)
%   simulation.debug_volumes   - intermediate volumes (NIfTI; default 1)
%
% Use as:
%   tf = debug_on(parameters)                  % master switch only
%   tf = debug_on(parameters, stage)           % stage selected
%   tf = debug_on(parameters, stage, kind)     % stage and kind selected
%
% Input:
%   parameters - PRESTUS config
%   stage      - 'preproc' | 'medium' | 'source' | 'babelbrain' ('' = any)
%   kind       - 'figures' | 'volumes' | 'other' ('' = any)
%
% Output:
%   tf - logical scalar
%
% See also: DEBUG_SAVE_FIGURE, DEBUG_REGISTER, DEBUG_STAGES, GENERATE_DEBUG_REPORT

arguments
    parameters (1,1) struct
    stage      {mustBeTextScalar} = ''
    kind       {mustBeTextScalar} = ''
end

    tf = isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') && ...
        isequal(double(parameters.simulation.debug), 1);
    if ~tf, return; end
    sim = parameters.simulation;

    if strlength(stage) > 0
        selected = debug_stages(parameters);
        tf = any(strcmp(char(stage), selected));
        if ~tf, return; end
    end

    switch char(kind)
        case {'', 'other'}
            tf = true;
        case 'figures'
            tf = flag(sim, 'debug_figures', 1);
        case 'volumes'
            tf = flag(sim, 'debug_volumes', 1);
        otherwise
            error('PRESTUS:debug:unknownKind', ...
                'Unknown debug output kind ''%s'' (figures, volumes, other).', kind);
    end
end

function tf = flag(sim, name, default)
    tf = logical(default);
    if isfield(sim, name) && ~isempty(sim.(name))
        tf = isequal(double(sim.(name)), 1) || isequal(sim.(name), true);
    end
end
