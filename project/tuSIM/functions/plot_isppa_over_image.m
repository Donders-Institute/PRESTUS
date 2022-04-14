function [res_image, transducer_pars] = plot_isppa_over_image(Isppa_map, bg_image, transducer_bowl, after_exit_plane_mask, slice, trans_pos, focus_pos, max_isppa_pos)
    slice_x = 1:size(Isppa_map,1);
    slice_y = 1:size(Isppa_map,2);
    slice_z = 1:size(Isppa_map,3);

   
    if slice{1} == 'x'
        slice_x = slice{2};
        focus_pos = focus_pos(2:3);
        trans_pos = trans_pos(2:3);
        max_isppa_pos = max_isppa_pos(2:3);
    elseif slice{1} == 'y'
        slice_y = slice{2};
        focus_pos = focus_pos([1,3]);
        trans_pos = trans_pos([1,3]);
        max_isppa_pos = max_isppa_pos([1,3]);
    elseif slice{1} == 'z'
        slice_z = slice{2};
        focus_pos = focus_pos(1:2);
        trans_pos = trans_pos([1,2]);
        max_isppa_pos = max_isppa_pos([1,2]);
    else
        error("slice must be a cell array with the first element of 'x','y', or 'z' and a second element a slice number")
    end
    
    bg_slice = mat2gray(squeeze(bg_image(slice_x, slice_y, slice_z)));
    transducer_bowl = mat2gray(squeeze(transducer_bowl(slice_x, slice_y, slice_z)));
    before_exit_plane_mask = ~squeeze(after_exit_plane_mask(slice_x, slice_y, slice_z));
    Isppa_map = squeeze(Isppa_map(slice_x, slice_y, slice_z));

    figure;
    ax1 = axes;
    imagesc(bg_slice+transducer_bowl+0.2*before_exit_plane_mask);
    colormap(ax1,'gray');
    axis image
    axis off;
    ax2 = axes;
    imagesc(ax2,Isppa_map,'alphadata', (Isppa_map-min(Isppa_map))/max(Isppa_map(:)));
    colormap(ax2,'viridis');
    %caxis(ax2,[min(Isppa_map(:)) max(Isppa_map(:))]);
    ax2.Visible = 'off';
    %linkprop([ax1 ax2],'Position');
    axis image
    axis off;

    rect_size = 2;
    rectangle('Position', [trans_pos(2)-rect_size/2  trans_pos(1)-rect_size/2  rect_size*2+1 rect_size*2+1], 'EdgeColor', 'g', 'LineWidth',1,'LineStyle','-')
    rectangle('Position', [focus_pos(2)-rect_size/2 focus_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'r', 'LineWidth',1,'LineStyle','-')
    rectangle('Position', [max_isppa_pos(2)-rect_size/2 max_isppa_pos(1)-rect_size/2 rect_size*2+1 rect_size*2+1], 'EdgeColor', 'b', 'LineWidth',1,'LineStyle','-')
    %colorbar;

end
