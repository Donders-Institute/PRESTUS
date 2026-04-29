function tp_plot_candidate_positions(img, target, trans_candidate, ...
    pixel_size, parameters, subject_id, target_name)
    % TP_PLOT_CANDIDATE_POSITIONS  Plot transducer candidate positions on a skull surface slice
    %
    % Creates a bounds visualisation showing the target (red), a random
    % transducer candidate (blue), and the search sphere intersection with
    % the skull surface. Saves the figure to disk.
    %
    % Use as:
    %   tp_plot_candidate_positions(img, target, trans_candidate, ...
    %       pixel_size, parameters, subject_id, target_name)
    %
    % Input:
    %   img             - [Nx x Ny x Nz] tissue segmentation volume
    %   target          - [1x3] target coordinates in voxel space
    %   trans_candidate - struct with transducer candidate position fields
    %   pixel_size      - voxel size [mm]
    %   parameters      - (1,1) simulation parameters struct
    %   subject_id      - subject identifier used in output filename
    %   target_name     - target label string used in output filename
    %
    % See also: TP_FIND_INITIAL_CANDIDATE, TP_PLOT_GEOMETRY_OVERLAY

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
    output_plot = fullfile(parameters.io.figures_preproc_dir, ...
        sprintf('sub-%03d_bounds_%s.png', subject_id, target_name));
    saveas(h, output_plot, 'png');
    close(h);
end
