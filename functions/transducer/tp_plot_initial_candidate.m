function tp_plot_initial_candidate_positions(segm_img_slice, target, outer_sphere_3d, t1_x, t1_y, t1_z, ...
    pixel_size, parameters, subject_id, target_name)
    % TP_PLOT_INITIAL_CANDIDATE_POSITIONS Plot transducer candidate positions on skull surface slice
    %
    % Creates bounds visualization showing target (red), random transducer candidate (blue),
    % and search sphere intersection with skull surface.
    %
    % INPUT
    %   segm_img_slice  - RGB slice through target(y) with sphere highlighted
    %   target          - 1x3 target coordinates [x,y,z]
    %   outer_sphere_3d - 3D logical mask of skull surface candidates
    %   t1_x/y/z        - ndgrid coordinate arrays
    %   pixel_size      - Voxel size (mm)
    %   parameters      - Transducer parameters for get_transducer_box
    %   subject_id      - Scalar ID for filename
    %   target_name     - String for filename
    %
    % OUTPUT
    %   Saves: sub-XXX_bounds_TARGET.png

    h = figure;
    colormap([0.3 0.3 0.3; lines(12)])
    imagesc(segm_img_slice); axis image; hold on;

    % Target (red rectangle)
    target_xz = target([1,3]);
    rectangle('Position',[flip(target_xz)-2, 4, 4], 'Curvature',[0,0], ...
        'EdgeColor','r', 'LineWidth',2, 'LineStyle','-');

    % Random transducer candidate from skull surface intersection
    outer_idx = find(outer_sphere_3d&t1_y==target(2));
    trans_idx = randsample(outer_idx, 1);
    trans_pos = [t1_x(trans_idx), t1_y(trans_idx), t1_z(trans_idx)];
    trans_xz = trans_pos([1,3]);

    % Candidate position (blue rectangle)
    rectangle('Position',[flip(trans_xz)-2, 4, 4], 'Curvature',[0,0], ...
        'EdgeColor','b', 'LineWidth',2, 'LineStyle','-');

    % Transducer visualization box
    get_transducer_box(trans_xz, target_xz, pixel_size, parameters);

    % Save
    output_plot = fullfile(parameters.output_dir, sprintf('sub-%03d_bounds_%s.png', subject_id, target_name));
    saveas(h, output_plot, 'png');
    close(h);
end
