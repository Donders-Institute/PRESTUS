function varargout = live_view(action, parameters, varargin)
% LIVE_VIEW  Save the evolving fields of a simulation as frames for the Nexus
%
% With simulation.live_view = 1, frames of the evolving fields are saved at
% most every 5 s, also without the MATLAB desktop (batch runs, the Nexus,
% cluster jobs), each replacing the last:
%   - data frames ('data'): a 2D plane of the field as
%       <io.dir_output>/live/<stage>_sub-<id><affix>.f32  (float32, little-endian)
%     and its description, written last,
%       <io.dir_output>/live/<stage>_sub-<id><affix>.json
%     (quantity, units, spacing_mm, shape, order, min, max, markers as 0-based
%     plane indices). The Nexus draws them in its own theme. The thermal
%     simulation, jwave and BabelBrain write these.
%   - figure frames ('frame', 'watch'): the k-Wave MATLAB solvers plot the
%     pressure field while they run (PlotSim); its figure is saved as
%       <io.dir_output>/live/<stage>_sub-<id><affix>.png
% Each frame file is announced once with
%   [live] <stage> <path>
% from which the Nexus learns where to look.
%
% Use as:
%   on = live_view('enabled', parameters)
%   live_view('frame', parameters, fig, stage)          % save fig now (rate-limited)
%   live_view('frame', parameters, fig, stage, true)    % save fig now (always)
%   t = live_view('watch', parameters, stage)           % timer saving the newest figure
%   live_view('stop', parameters, t)                    % stop it and save a last frame
%   live_view('data', parameters, stage, plane, meta)    % save a data frame (rate-limited)
%   live_view('data', parameters, stage, plane, meta, true)   % ... always
%   file = live_view('file', parameters, stage)         % data frame (.json) path ('' when off),
%                                                       % for solvers that save frames themselves
%
% Input:
%   action     - 'enabled' | 'frame' | 'watch' | 'stop' | 'data' | 'file'
%   parameters - PRESTUS config (simulation.live_view, io.dir_output, subject_id,
%                io.output_affix)
%   plane      - 2D field, or a function handle returning it (evaluated only
%                when a frame is due, e.g. a gather from the GPU)
%   meta       - struct: quantity ('pressure' | 'temperature'), units, label,
%                spacing_mm ([2] mm), markers (struct of [2] 0-based indices)
%
% Output:
%   'enabled': logical; 'watch': timer object ([] when live view is off);
%   'file': char path of the stage's frame ('' when live view is off)
%
% Note:
%   'watch' captures k-Wave's own figure: MATLAB runs the timer when k-Wave
%   calls drawnow while plotting, so frames follow the solver's plot updates.
%   Plotting slows the solver down; live view is meant for checking a setup,
%   not for production runs. BabelBrain saves its frames itself (bb_progress.py,
%   environment variable PRESTUS_LIVE_FILE from 'file').
%
% See also: PRESTUS_PROGRESS, ACOUSTIC_SIMULATION, THERMAL_SIMULATION

    persistent LAST_WRITE ANNOUNCED
    if isempty(LAST_WRITE); LAST_WRITE = containers.Map(); end
    if isempty(ANNOUNCED); ANNOUNCED = containers.Map(); end
    min_interval_s = 5;

    switch action
        case 'enabled'
            varargout{1} = is_on(parameters);

        case 'file'
            varargout{1} = '';
            if is_on(parameters)
                varargout{1} = data_file(parameters, varargin{1});
            end

        case 'data'
            [stage, plane, meta] = varargin{1:3};
            force = numel(varargin) >= 4 && varargin{4};
            if ~is_on(parameters)
                return
            end
            file = data_file(parameters, stage);
            if ~force && LAST_WRITE.isKey(file) && toc(LAST_WRITE(file)) < min_interval_s
                return
            end
            LAST_WRITE(file) = tic;
            if write_data(plane, meta, file) && ~ANNOUNCED.isKey(file)
                ANNOUNCED(file) = true;
                fprintf('[live] %s %s\n', stage, file);
            end

        case 'frame'
            [fig, stage] = varargin{1:2};
            force = numel(varargin) >= 3 && varargin{3};
            if ~is_on(parameters) || isempty(fig) || ~isgraphics(fig, 'figure')
                return
            end
            file = frame_file(parameters, stage);
            if ~force && LAST_WRITE.isKey(file) && toc(LAST_WRITE(file)) < min_interval_s
                return
            end
            LAST_WRITE(file) = tic;
            if write_frame(fig, file) && ~ANNOUNCED.isKey(file)
                ANNOUNCED(file) = true;
                fprintf('[live] %s %s\n', stage, file);
            end

        case 'watch'
            stage = varargin{1};
            varargout{1} = [];
            if ~is_on(parameters)
                return
            end
            t = timer('Period', min_interval_s, 'ExecutionMode', 'fixedSpacing', ...
                'BusyMode', 'drop', 'StartDelay', min_interval_s, 'UserData', stage, ...
                'TimerFcn', @(~, ~) live_view('frame', parameters, newest_figure(), stage));
            start(t);
            varargout{1} = t;

        case 'stop'
            t = varargin{1};
            if isempty(t) || ~isvalid(t)
                return
            end
            stage = t.UserData;
            stop(t); delete(t);
            if ~isempty(stage)
                live_view('frame', parameters, newest_figure(), stage, true);
            end

        otherwise
            error('live_view:action', 'Unknown action ''%s''.', action);
    end
end

function on = is_on(parameters)
    on = isfield(parameters, 'simulation') && isfield(parameters.simulation, 'live_view') && ...
        ~isempty(parameters.simulation.live_view) && logical(parameters.simulation.live_view(1));
end

function file = frame_file(parameters, stage)
    file = fullfile(char(parameters.io.dir_output), 'live', ...
        sprintf('%s_sub-%03d%s.png', stage, parameters.subject_id, char(parameters.io.output_affix)));
end

function file = data_file(parameters, stage)
    file = frame_file(parameters, stage);
    file = [file(1:end-4) '.json'];
end

function ok = write_data(plane, meta, file)
% The plane as float32 (column-major) and its description, each via a
% temporary file and a rename; the description last, so a reader that sees it
% finds the matching data.
    ok = false;
    try
        if isa(plane, 'function_handle')
            plane = plane();
        end
        plane = single(gather(squeeze(plane)));
        folder = fileparts(file);
        if ~isfolder(folder); mkdir(folder); end
        raw = [file(1:end-5) '.f32'];
        fid = fopen([raw '.part'], 'w', 'ieee-le');
        fwrite(fid, plane, 'single');
        fclose(fid);
        movefile([raw '.part'], raw, 'f');
        [~, raw_name, raw_ext] = fileparts(raw);
        meta.data = [raw_name raw_ext];
        meta.shape = size(plane);
        meta.order = 'F';
        meta.min = double(min(plane(:)));
        meta.max = double(max(plane(:)));
        fid = fopen([file '.part'], 'w', 'n', 'UTF-8');
        fprintf(fid, '%s', jsonencode(meta));      % in the file's encoding ('°C')
        fclose(fid);
        movefile([file '.part'], file, 'f');
        ok = true;
    catch me
        warning('live_view:data', 'Live view frame not written (%s).', me.message);
    end
end

function fig = newest_figure()
% The most recently created figure that is not a waitbar (k-Wave opens one
% next to its plot); [] when there is none.
    fig = [];
    figs = findall(groot, 'Type', 'figure');
    figs = figs(~strcmp(get(figs, 'Tag'), 'TMWWaitbar'));
    if isempty(figs)
        return
    end
    nums = arrayfun(@figure_number, figs);   % NaN for figures without an integer handle
    if all(isnan(nums))
        fig = figs(1);           % findall lists the newest first
    else
        [~, i] = max(nums);
        fig = figs(i);
    end
end

function n = figure_number(f)
% The figure's number; NaN when it has none (IntegerHandle off: Number is
% empty, which arrayfun's uniform output cannot take).
    n = NaN;
    if isgraphics(f, 'figure') && ~isempty(f.Number)
        n = double(f.Number);
    end
end

function ok = write_frame(fig, file)
% Write via a temporary file and rename, so a reader never sees half a PNG.
    ok = false;
    try
        folder = fileparts(file);
        if ~isfolder(folder); mkdir(folder); end
        tmp = [file(1:end-4) '.part.png'];
        if exist('exportgraphics', 'file')
            exportgraphics(fig, tmp, 'Resolution', 96);
        else
            print(fig, tmp, '-dpng', '-r96');
        end
        movefile(tmp, file, 'f');
        ok = true;
    catch me
        warning('live_view:write', 'Live view frame not written (%s).', me.message);
    end
end
