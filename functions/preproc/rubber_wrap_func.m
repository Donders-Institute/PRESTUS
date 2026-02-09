function rubber_wrap_func()

    %fnGZ = "/Users/sjoerd.meijer/Documents/pCT/skull_mask_PVC.nii.gz";
    %subj = 'sub-098'
    
    subj = 'XXX';

    % Main directories
    main_dir = fullfile('/project/3017000.01/TUS/mareng/TUSdata_pilot', subj, 'masks');
    fnGZ = fullfile(main_dir, 'skull_mask_PVC.nii.gz');
   

    % ==============================================================================
    % Load NIfTI files
    % ==============================================================================

    %% Load NIfTI (.nii.gz supported)
    fn = fnGZ; tmpDir = "";
    if endsWith(lower(fnGZ), ".gz")
        tmpDir = tempname; mkdir(tmpDir);
        gunzip(fnGZ, tmpDir);
        [~, base, ~] = fileparts(fnGZ);         % base ends with .nii
        fnNii = fullfile(tmpDir, base);
        if ~isfile(fnNii)
            d = dir(fullfile(tmpDir, "*.nii"));
            fnNii = fullfile(tmpDir, d(1).name);
        end
        fn = string(fnNii);
    end

    info = niftiinfo(fn);
    BW = niftiread(info) ~= 0;     % 1 = skull
    sz = size(BW);
    fprintf("Volume size: %s\n", mat2str(sz));

    % ==============================================================================
    % CLEAN SKULL MASK
    % ==============================================================================
    
    %% Keep largest connected skull component (recommended)
    CC = bwconncomp(BW, 26);
    if CC.NumObjects == 0, error("Skull mask is empty."); end
    sizes = cellfun(@numel, CC.PixelIdxList);
    [~, iBig] = max(sizes);
    BW1 = false(sz);
    BW1(CC.PixelIdxList{iBig}) = true;

    % ==============================================================================
    % RUBBER WRAP AROUND SKULL = BALLOON MASK
    % ==============================================================================

    %% --- Rubber-wrap parameters ---
    wrapRadius = 10;   % voxels. Bigger = tighter rubber that ignores bigger dents (try 2–10)

    se = strel("sphere", wrapRadius);

    % "Shrink-wrap" / morphological hull:
    % Closing bridges concavities smaller than wrapRadius (rubber won't go in).
    BWwrap = imclose(BW1, se);

    % Make sure the wrap is a solid region (optional but usually good)
    BWwrap = imfill(BWwrap, "holes");

    % Optional: remove any stray bits and keep largest component again
    CCw = bwconncomp(BWwrap, 26);
    sizesw = cellfun(@numel, CCw.PixelIdxList);
    [~, iw] = max(sizesw);
    BWwrap2 = false(sz);
    BWwrap2(CCw.PixelIdxList{iw}) = true;

    BWwrap = BWwrap2;

    % Choose balloon mask to save
    BWsave = BWwrap;   % or BWwrap_half

    % Ensure logical / uint8
    BWsave = uint8(BWsave);   % NIfTI likes numeric types

    %outNii = "/Users/sjoerd.meijer/Documents/pCT/balloon_mask.nii";
    outNii = fullfile(main_dir, 'balloon_mask.nii')


    infoOut = info;                 % reuse original header
    infoOut.Datatype = 'uint8';
    infoOut.BitsPerPixel = 8;
    infoOut.Filename = outNii;
    infoOut.Filemoddate = char(datetime("now"));

    niftiwrite(BWsave, outNii, infoOut, 'Compressed', false);

    gzip(outNii);
    delete(outNii);   % optional: remove uncompressed .nii


    %% balloon_touch_GM_and_SKIN_voxelwise.m

    %balloonPath = "/Users/sjoerd.meijer/Documents/pCT/balloon_mask.nii.gz";
    balloonPath = fullfile(main_dir, 'balloon_mask.nii.gz')

    %tissuePath  = "/Users/sjoerd.meijer/Documents/pCT/final_tissues_resamp1mm_crop.nii.gz";
    tissuePath = fullfile(main_dir, 'final_tissues_resamp1mm_crop.nii.gz')

    GM_LABEL   = 57;
    SKIN_LABEL = 142;

    touchRadius = 1;   % voxels: 1 = immediate neighbors (26-neighborhood)

    %outNii = "/Users/sjoerd.meijer/Documents/pCT/balloon_mask_touch_GM_AND_SKIN_noGM.nii";
    outNii = fullfile(main_dir, 'balloon_mask_touch_GM_AND_SKIN_noGM.nii')

    %% Load balloon
    infoB = niftiinfo(balloonPath);
    B = niftiread(infoB) ~= 0;

    %% Load tissues
    T = niftiread(niftiinfo(tissuePath));

    if ~isequal(size(B), size(T))
        error("Size mismatch: balloon %s vs tissues %s", mat2str(size(B)), mat2str(size(T)));
    end

    GM   = (T == GM_LABEL);
    SKIN = (T == SKIN_LABEL);

    %% "Touch" masks: balloon voxel touches GM if it lies within 1-voxel dilation of GM
    se = strel("sphere", touchRadius);

    GM_touch_region   = imdilate(GM,   se);
    SKIN_touch_region = imdilate(SKIN, se);

    % Keep balloon voxels that touch BOTH
    B_touch_both = B & GM_touch_region & SKIN_touch_region;

    % Remove overlap with GM from the remaining mask
    B_out = B_touch_both & ~GM;

    %% Include MEN(85) voxels at the MEN–SKIN interface (strict 6-neigh), near B_out

    MENINGES_LABEL = 85;
    SKIN_LABEL     = 142;

    MEN  = (T == MENINGES_LABEL);
    SKIN = (T == SKIN_LABEL);

    % --- Strict 6-neighborhood kernel (faces only) ---
    K6 = zeros(3,3,3,'logical');
    K6(2,2,2) = true;
    K6(1,2,2) = true; K6(3,2,2) = true;
    K6(2,1,2) = true; K6(2,3,2) = true;
    K6(2,2,1) = true; K6(2,2,3) = true;

    % SKIN dilated by 6-neighborhood => voxels face-adjacent to SKIN
    SKIN_touch_region_6 = convn(single(SKIN), single(K6), 'same') > 0;

    % MEN voxels that touch SKIN (faces only)
    MEN_touch_SKIN_6 = MEN & SKIN_touch_region_6;

    fprintf("MEN voxels touching SKIN (6-neigh): %d\n", nnz(MEN_touch_SKIN_6));

    % --- Restrict to MEN voxels near existing balloon mask ---
    nearRadius = 2;                 % 1–2 suggested
    seNear = strel("sphere", nearRadius);
    nearBalloon = imdilate(B_out, seNear);

    MEN_touch_SKIN_6_near = MEN_touch_SKIN_6 & nearBalloon;

    fprintf("...of those, near existing B_out (r=%d): %d\n", nearRadius, nnz(MEN_touch_SKIN_6_near));

    % --- Add MEN interface voxels to balloon mask ---
    B_out2 = B_out | MEN_touch_SKIN_6_near;

    % Keep your original constraint: no GM overlap
    B_out2 = B_out2 & ~GM;

    B_out = BW1 | B_out2;

    %% 1) Fill along Z where there is skull in front AND behind (same Y,X)
    B0 = (B_out ~= 0);   % canonical logical

    BWm1 = circshift(B0, [0 0  1]);
    BWp1 = circshift(B0, [0 0 -1]);

    BWm1(:,:,1)   = false;
    BWp1(:,:,end) = false;

    fill1 = ~B0 & BWm1 & BWp1;

    B1 = B0 | fill1;     % final filled mask (logical)

    %% Save as .nii.gz

    BWsave = uint8(B1);

    infoOut = infoB;
    infoOut.Datatype = 'uint8';
    infoOut.BitsPerPixel = 8;
    infoOut.Filename = outNii;
    infoOut.Filemoddate = char(datetime("now"));

    niftiwrite(BWsave, outNii, infoOut, 'Compressed', false);
    gzip(outNii);
    delete(outNii);

    fprintf("Saved: %s.gz\n", outNii);

    %% --- Inputs (already in your workspace) ---
    SKULL   = BW1;
    BALLOON = B0 & ~SKULL;   % balloon-approach voxels
    ZADDED  = fill1;         % Z-axis added voxels

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
    lgd = legend(axLegend, {'UTE-derived pseudo-CT', 'Morphological shrink-wrap skull surface completion'}, ... 
        'Location','southoutside', 'Orientation','horizontal'); lgd.TextColor = 'w'; lgd.Color = 'k'; lgd.Box = 'off'; 

    lgd.FontSize = 16;          % try 16–18
    lgd.ItemTokenSize = [36 20];  % makes the legend symbols wider/taller

    %% Save figure as PNG
    output_png = fullfile(main_dir, 'skull_visualization.png');
    print(fig, output_png, '-dpng', '-r300');  % 300 DPI resolution
    fprintf("Figure saved: %s\n", output_png);

    close all
end