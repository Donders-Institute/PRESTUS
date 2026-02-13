function skull_rubber_wrap_visualize(parameters, SKULL, ZADDED, BALLOON, downsample_factor)
% skull_rubber_wrap_visualize
% Faster + more robust multi-view 3D visualization for HPC:
% - Default downsample_factor = 1
% - Uses tiledlayout (faster/cleaner axes creation) [web:126][web:139]
% - Uses explicit patch(Faces,Vertices,...) to avoid "Not enough input arguments" issues
% - Reduces rendering load: no edges, no gouraud, optional point cap
%
% Inputs:
%   parameters.debug_dir
%   parameters.results_filename_affix
%   SKULL, ZADDED, BALLOON : 3-D logical or numeric masks (same size)
%   downsample_factor      : (optional) integer >=1

    disp('Plotting skull rubber expansion');

    if nargin < 5 || isempty(downsample_factor)
        downsample_factor = 1;
    end

    %% Combine masks

    SKULL = logical(SKULL);
    ZADDED = logical(ZADDED);
    BALLOON = logical(BALLOON);

    %% Downsample (nearest neighbor for label/binary masks)

    if downsample_factor > 1
        if exist('imresize3','file') == 2
            s = 1 / downsample_factor;
            SKULL_DS   = imresize3(SKULL,   s, 'nearest');
            ZADDED_DS  = imresize3(ZADDED,  s, 'nearest');
            BALLOON_DS = imresize3(BALLOON, s, 'nearest');
        else
            SKULL_DS   = SKULL(  1:downsample_factor:end, 1:downsample_factor:end, 1:downsample_factor:end);
            ZADDED_DS  = ZADDED( 1:downsample_factor:end, 1:downsample_factor:end, 1:downsample_factor:end);
            BALLOON_DS = BALLOON(1:downsample_factor:end, 1:downsample_factor:end, 1:downsample_factor:end);
        end
    else
        SKULL_DS   = SKULL;
        ZADDED_DS  = ZADDED;
        BALLOON_DS = BALLOON;
    end

    SKULL_DISPLAY = SKULL_DS | ZADDED_DS;

    %% Isosurface once

    fv = isosurface(single(SKULL_DISPLAY), 0.5);  % fv.faces, fv.vertices

    % Guard: empty surface (e.g., mask empty)
    if ~isfield(fv,'faces') || isempty(fv.faces) || ~isfield(fv,'vertices') || isempty(fv.vertices)
        warning('skull_rubber_wrap_visualize:EmptySurface', 'Isosurface empty; skipping plot.');
        return;
    end

    F = fv.faces;
    V = fv.vertices;

    %% Balloon voxel coords (cap for speed)

    [iyB, ixB, izB] = ind2sub(size(BALLOON_DS), find(BALLOON_DS));
    xB = ixB; yB = iyB; zB = izB;

    maxPts = 8000;            % adjust to taste
    if numel(xB) > maxPts
        idx = randperm(numel(xB), maxPts);
        xB = xB(idx); yB = yB(idx); zB = zB(idx);
    end

    %% Views

    views_left = { ...
        [90    90],   'Front view';  ...
        [-90   0],    'Top view';    ...
        [90    0],    'Bottom view'  ...
    };

    views_right = { ...
        [0     45],   'Right side view'; ...
        [-180  45],   'Left side view'   ...
    };

    %% Figure + layout

    set(0, 'DefaultFigureRenderer', 'painters'); % HPC-safe, avoids OpenGL issues
    fig = figure('Color','k','Position',[100 100 1800 800], 'PaperPositionMode','auto', 'Visible','off');

    t = tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact'); % [web:126]
    % Left column
    axL(1) = nexttile(t, 1);  % [web:139]
    axL(2) = nexttile(t, 4);
    axL(3) = nexttile(t, 7);
    % Right: two big tiles (2x2 + 1x2)
    axR(1) = nexttile(t, 2, [2 2]); % spans rows 1-2, cols 2-3 [web:139]
    axR(2) = nexttile(t, 8, [1 2]); % spans row 3, cols 2-3 [web:139]

    set([axL axR], 'Color','k');

    %% Styling

    skullFaceColor = [0.92 0.92 0.92];

    blueColor      = [0 0.45 0.9];
    markerSmall    = 10;
    markerBig      = 16;

    %% Helper: draw one view into an axes (explicit Faces/Vertices)

    function draw_view(ax, viewvec, ttl, msize)
        cla(ax);

        pSkull = patch(ax, 'Faces', F, 'Vertices', V, ...
            'FaceColor', skullFaceColor, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.85, ...
            'FaceLighting', 'none', ...
            'BackFaceLighting', 'reverselit');

        % material(pSkull,'dull');

        hold(ax,'on');

        if ~isempty(xB)
            scatter3(ax, xB, yB, zB, msize, 'filled', 'MarkerFaceColor', blueColor);
        end

        axis(ax,'equal'); axis(ax,'tight'); axis(ax,'off');
        view(ax, viewvec);
        title(ax, ttl, 'Color','w');

        % If you want a tiny bit of shading without heavy cost:
        lighting(ax,'flat'); camlight(ax,'headlight');  % optional
    end

    %% Render
    
    for k = 1:3
        draw_view(axL(k), views_left{k,1}, views_left{k,2}, markerSmall);
    end
    for k = 1:2
        draw_view(axR(k), views_right{k,1}, views_right{k,2}, markerBig);
    end

    %% Legend (attach to bottom-right axes)

    lgd = legend(axR(2), {'Skull surface','Balloon voxels'}, ...
        'Location','southoutside', 'Orientation','horizontal');
    lgd.TextColor = 'w';
    lgd.Color = 'k';
    lgd.Box = 'off';

    %% Fast raster export (avoid slow vector pipeline)

    output_png = fullfile(parameters.debug_dir, ...
        sprintf('skull_visualization%s.png', parameters.results_filename_affix));

    drawnow limitrate;
    fr = getframe(fig);
    imwrite(fr.cdata, output_png, 'png', 'Compression','none');

    close(fig);
end
