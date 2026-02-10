function skull_rubber_wrap_visualize(parameters, SKULL, ZADDED, BALLOON)

    %% --- Inputs (already in your workspace) ---
    % --- Skull for display: original skull + Z-added voxels ---
    SKULL_DISPLAY = SKULL | ZADDED;

    fprintf("Skull (displayed) voxels : %d\n", nnz(SKULL_DISPLAY));
    fprintf("Balloon voxels (blue)    : %d\n", nnz(BALLOON));

    %% Skull surface
    fv = isosurface(single(SKULL_DISPLAY), 0.5);

    %% Balloon voxel coordinates (blue)
    [iyB, ixB, izB] = ind2sub(size(BALLOON), find(BALLOON));
    xB = ixB; yB = iyB; zB = izB;

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

    %% Figure + layout (black background; side views large, next to each other)
    fig = figure('Color','k','Position',[100 100 1800 800]);

    pLeft  = uipanel(fig, 'Units','normalized', 'Position',[0.02 0.08 0.28 0.90], 'BorderType','none', 'BackgroundColor','k');
    pRight = uipanel(fig, 'Units','normalized', 'Position',[0.32 0.08 0.66 0.90], 'BorderType','none', 'BackgroundColor','k');

    leftTL  = tiledlayout(pLeft,  3, 1, 'Padding','compact', 'TileSpacing','compact');
    rightTL = tiledlayout(pRight, 1, 2, 'Padding','compact', 'TileSpacing','compact');

    %% Visual styling
    skullFaceColor = [0.92 0.92 0.92];   % slightly off-white (better shading)
    skullEdgeColor = [0.15 0.15 0.15];   % subtle contour lines
    skullEdgeAlpha = 0.25;              % keep edges faint
    skullLineWidth = 0.2;               % thin outline

    blueColor      = [0 0.45 0.9];
    markerSmall    = 10;
    markerBig      = 16;

    %% --- Left column (3 smaller views) ---
    for k = 1:3
        ax = nexttile(leftTL);
        set(ax,'Color','k');

        % Skull surface
        pSkull = patch(ax, fv);
        pSkull.FaceColor = skullFaceColor;
        pSkull.EdgeColor = skullEdgeColor;   % <-- subtle edges for contours
        pSkull.LineWidth = skullLineWidth;
        pSkull.EdgeAlpha = skullEdgeAlpha;
        pSkull.FaceAlpha = 1.0;
        pSkull.BackFaceLighting = 'reverselit';
        material(pSkull,'dull');

        hold(ax,'on');

        % Balloon voxels (blue)
        scatter3(ax, xB, yB, zB, markerSmall, 'filled', 'MarkerFaceColor', blueColor);

        axis(ax,'equal'); axis(ax,'tight'); axis(ax,'off');
        view(ax, views_left{k,1});
        title(ax, views_left{k,2}, 'Color','w');

        % Lighting: two lights gives depth on dark bg
        camlight(ax,'headlight');
        camlight(ax,'right');
        lighting(ax,'gouraud');

        % Better normals for crisper shading
        isonormals(single(SKULL_DISPLAY), pSkull);
    end

    %% --- Right panel (2 BIG side views next to each other) ---
    for k = 1:2
        ax = nexttile(rightTL);
        set(ax,'Color','k');

        % Skull surface
        pSkull = patch(ax, fv);
        pSkull.FaceColor = skullFaceColor;
        pSkull.EdgeColor = skullEdgeColor;
        pSkull.LineWidth = skullLineWidth;
        pSkull.EdgeAlpha = skullEdgeAlpha;
        pSkull.FaceAlpha = 1.0;
        pSkull.BackFaceLighting = 'reverselit';
        material(pSkull,'dull');

        hold(ax,'on');

        % Balloon voxels (blue) - larger markers
        scatter3(ax, xB, yB, zB, markerBig, 'filled', 'MarkerFaceColor', blueColor);

        axis(ax,'equal'); axis(ax,'tight'); axis(ax,'off');
        view(ax, views_right{k,1});
        title(ax, views_right{k,2}, 'Color','w');

        camlight(ax,'headlight');
        camlight(ax,'right');
        lighting(ax,'gouraud');
        isonormals(single(SKULL_DISPLAY), pSkull);
    end

    %% Legend (white text on black), attach to one right-panel axis
    % Create invisible axes spanning the whole figure
    axRight = findall(pRight,'Type','axes'); 
    axLegend = axRight(end); 
    lgd = legend(axLegend, {'Input skull mask', 'Morphological shrink-wrap skull surface completion'}, ... 
        'Location','southoutside', 'Orientation','horizontal'); lgd.TextColor = 'w'; lgd.Color = 'k'; lgd.Box = 'off'; 

    lgd.FontSize = 16;          % try 16–18
    lgd.ItemTokenSize = [36 20];  % makes the legend symbols wider/taller

    %% Save figure as PNG
    output_png = fullfile(parameters.debug_dir, ...
        sprintf('skull_visualization%s.png', parameters.results_filename_affix));
    print(fig, output_png, '-dpng', '-r300');  % 300 DPI resolution
    fprintf("Figure saved: %s\n", output_png);

    close(fig);