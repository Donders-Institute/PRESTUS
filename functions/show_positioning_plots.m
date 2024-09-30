function show_positioning_plots(segmented_img_orig, t1_pixel_size, trans_pos_orig, focus_pos_orig, segmented_img_final, trans_pos_final, focus_pos_final, parameters, output_plot)
    arguments
        segmented_img_orig (:,:,:) 
        t1_pixel_size (1,1)
        trans_pos_orig (1,3)
        focus_pos_orig (1,3)
        segmented_img_final (:,:,:) 
        trans_pos_final (1,3)
        focus_pos_final (1,3)
        parameters struct
        output_plot string
    end

    view_pos = [-180, 0];
    slice_cap = [1,0,0];
    if trans_pos_final(3) > size(segmented_img_final,3)/2
        view_pos = [0,0];
        slice_cap = [-1,0,0];
    end
    coord_mesh_xyz = get_xyz_mesh(segmented_img_orig);
    if gpuDeviceCount>0
        coord_mesh_xyz = gpuArray(coord_mesh_xyz);
    end

    h = figure('Position', [200 200 1000 300]);
    subplot(1,3,1)
    show_3d_head(segmented_img_orig, focus_pos_orig, trans_pos_orig, parameters, ...
        t1_pixel_size, coord_mesh_xyz, [0 0 0], view_pos, 0)

    subplot(1,3,2)
    show_3d_head(segmented_img_orig, focus_pos_orig, trans_pos_orig, parameters, ...
        t1_pixel_size, coord_mesh_xyz, slice_cap, view_pos, 0)

    ax3 = subplot(1,3,3);
    imagesc(squeeze(segmented_img_final(:,round(trans_pos_final(2)),:)))
    set(ax3,'dataAspectRatio',[1 1 1])
    rectangle('Position',[focus_pos_final([3,1]) - 2, 4 4],...
              'Curvature',[0,0], 'EdgeColor','r',...
             'LineWidth',2,'LineStyle','-');

    rectangle('Position',[trans_pos_final([3,1]) - 2, 4 4],...
              'Curvature',[0,0], 'EdgeColor','b',...
             'LineWidth',2,'LineStyle','-');

    line([trans_pos_final(3) focus_pos_final(3)], [trans_pos_final(1) focus_pos_final(1)], 'Color', 'white')
    get_transducer_box(trans_pos_final([1,3])', focus_pos_final([1,3])', parameters.grid_step_mm, parameters)
    colormap(ax3, [0.3 0.3 0.3; lines(12)])

    if output_plot ~= ""
        saveas(h, output_plot, 'png')
        close(h);
    end
end