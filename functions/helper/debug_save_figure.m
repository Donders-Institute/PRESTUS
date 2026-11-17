function png_files = debug_save_figure(parameters, stage, h, name, caption, options)
% DEBUG_SAVE_FIGURE  Save a debug figure, register it, and close it
%
% Writes <io.dir_debug_<stage>>/<name>.png, records it in the debug
% manifest (DEBUG_REGISTER) and closes the figure. Callers check
% DEBUG_ON(parameters, stage, 'figures') before drawing, so no figure is
% drawn that is not saved.
%
% A 3D figure is saved from several camera angles instead ('Views'): one
% PNG per view, <name>_view-<view>.png, registered as one figure group so
% the debug report shows the views side by side. Views:
%   beam    - along z, the beam axis in rotated grids (view [0 90])
%   side    - across the beam, x-z plane (view [0 0])
%   oblique - overview (view [30 25])
% Lights made with camlight follow the camera to each view.
%
% Use as:
%   png_files = debug_save_figure(parameters, stage, h, name, caption)
%   png_files = debug_save_figure(..., 'Views', {'beam', 'side', 'oblique'})
%   png_files = debug_save_figure(..., 'Resolution', 150, 'Close', false)
%
% Input:
%   parameters - PRESTUS config (io.dir_debug_<stage>)
%   stage      - debug stage (DEBUG_STAGES)
%   h          - figure handle
%   name       - file name without extension
%   caption    - one-line description for the debug report
%   options:
%     Views      - cellstr of views for a 3D figure ({} = save as drawn)
%     Resolution - dpi for exportgraphics; [] (default) saves with saveas
%     Close      - close the figure afterwards (default true)
%     Dir        - folder to write to instead of io.dir_debug_<stage>
%
% Output:
%   png_files - cellstr of the PNGs written ({} when none could be written)
%
% See also: DEBUG_ON, DEBUG_REGISTER

arguments
    parameters (1,1) struct
    stage      {mustBeTextScalar}
    h          (1,1)
    name       {mustBeTextScalar}
    caption    {mustBeTextScalar} = ''
    options.Views      cell = {}
    options.Resolution = []
    options.Close      (1,1) logical = true
    options.Dir        {mustBeTextScalar} = ''
end

    known = struct('beam', [0 90], 'side', [0 0], 'oblique', [30 25]);
    labels = struct('beam', 'view along z, the beam axis in rotated grids', ...
        'side', 'side view, x-z plane', 'oblique', 'oblique view');
    bad = setdiff(options.Views, fieldnames(known));
    if ~isempty(bad)
        error('PRESTUS:debug:unknownView', 'Unknown debug figure view(s) %s; known: %s.', ...
            strjoin(bad, ', '), strjoin(fieldnames(known)', ', '));
    end

    png_files = {};
    out_dir = char(options.Dir);
    if isempty(out_dir)
        field = ['dir_debug_' char(stage)];
        if isfield(parameters, 'io') && isfield(parameters.io, field)
            out_dir = parameters.io.(field);
        end
    end
    try
        if isempty(out_dir)
            error('PRESTUS:debug:noDir', 'io.dir_debug_%s is not set.', stage);
        end
        if ~isfolder(out_dir); mkdir(out_dir); end
        if isempty(options.Views)
            png_files{end+1} = save_png(h, fullfile(out_dir, [char(name) '.png']), options.Resolution);
            debug_register(parameters, stage, png_files{end}, caption, 'figure');
        else
            axs = findall(h, 'Type', 'axes');
            for v = options.Views(:)'
                for ax = axs(:)'
                    view(ax, known.(v{1}));
                    for lt = findall(ax, 'Type', 'light')'
                        camlight(lt, 'headlight');
                    end
                end
                drawnow;
                png_files{end+1} = save_png(h, fullfile(out_dir, sprintf('%s_view-%s.png', name, v{1})), ...
                    options.Resolution); %#ok<AGROW>
                debug_register(parameters, stage, png_files{end}, ...
                    sprintf('%s (%s)', caption, labels.(v{1})), 'figure', char(name));
            end
        end
    catch ME
        warning('PRESTUS:debug:figure', 'Could not save debug figure %s: %s', char(name), ME.message);
    end
    if options.Close && isgraphics(h)
        close(h);
    end
end

function f = save_png(h, f, resolution)
    if isempty(resolution)
        saveas(h, f, 'png');
    else
        exportgraphics(h, f, 'Resolution', resolution);
    end
end
