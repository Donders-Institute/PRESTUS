function plot_coronal_slices(img, coord)
arguments
    img (:,:,:)
    coord (3,1) = round(size(img)/2)
end
img_rescaled = round(rescale(img, 1, 255));
figure
montage({squeeze(img_rescaled(coord(1),:,:))',squeeze(img_rescaled(:,coord(2),:))',squeeze(img_rescaled(:,:,coord(3)))'},gray, 'Size',[1 3])
end