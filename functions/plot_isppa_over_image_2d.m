function [res_image, transducer_pars, h] = plot_isppa_over_image_2d(...
    Isppa_map, ...
    bg_image, ...
    transducer_bowl, ...
    after_exit_plane_mask, ...
    trans_pos, focus_pos, ...
    max_isppa_pos)

    bg_slice = mat2gray(bg_image);
    transducer_bowl = mat2gray(transducer_bowl);
    before_exit_plane_mask = ~after_exit_plane_mask;

    h = figure;
    ax1 = axes;
    imagesc(bg_slice+transducer_bowl+0.2*before_exit_plane_mask);
    colormap(ax1,'gray');
    axis image
    axis off;
    ax2 = axes;
    imagesc(ax2,Isppa_map,'alphadata', (Isppa_map-min(Isppa_map))/max(Isppa_map(:)));
    colormap(ax2,'viridis');
%     if exist("clim")==2 % renamed in R2022a
%         clim(ax2,[min(Isppa_map(:)) max(Isppa_map(:))]);
%     elseif exist("caxis")==2
%         caxis(ax2,[min(Isppa_map(:)) max(Isppa_map(:))]);
%     end
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
