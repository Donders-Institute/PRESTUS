function tp_plot_candidate_positions(img, target, trans_candidate, ...
    pixel_size, parameters, subject_id, target_name)
    % TP_PLOT_CANDIDATE_POSITIONS Plot transducer candidate positions on skull surface slice
    %
    % Creates bounds visualization showing target (red), random transducer candidate (blue),
    % and search sphere intersection with skull surface.
    %
    % INPUT
    %   img             - 3D tissue segmentation image
    %   target          - 1x3 target coordinates [x,y,z]
    %   trans_candidate - Structure with info on transducer candidate
    %   pixel_size      - Voxel size (mm)
    %   parameters      - Transducer parameters for get_transducer_box
    %   subject_id      - Scalar ID for filename
    %   target_name     - String for filename
    %
    % OUTPUT
    %   Saves: sub-XXX_bounds_TARGET.png

    disp("[TP] Plotting initial candidate position ...")

    trans_xz = trans_candidate.trans_xz;
    target_xyz = trans_candidate.trans_pos;

    % Get image slice
    img_slice = ind2rgb(squeeze(img(:,target_xyz(2),:)), viridis(max(img(:))+1));

    h = figure;
    colormap([0.3 0.3 0.3; lines(12)])
    imagesc(img_slice); axis image; hold on;

    % Target (red rectangle)
    target_xz = target([1,3]);
    rectangle('Position',[flip(target_xz)-2, 4, 4], 'Curvature',[0,0], ...
        'EdgeColor','r', 'LineWidth',2, 'LineStyle','-');

    % Candidate position (blue rectangle)
    rectangle('Position',[flip(trans_xz)-2, 4, 4], 'Curvature',[0,0], ...
        'EdgeColor','b', 'LineWidth',2, 'LineStyle','-');

    % Transducer visualization box
    get_transducer_box(trans_xz, target_xz, [], pixel_size, parameters);

    % Save
    output_plot = fullfile(parameters.io.output_dir, ...
        sprintf('sub-%03d_bounds_%s.png', subject_id, target_name));
    saveas(h, output_plot, 'png');
    close(h);
end
