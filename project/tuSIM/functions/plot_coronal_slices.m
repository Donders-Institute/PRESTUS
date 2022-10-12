function plot_coronal_slices(img)
im_center = round(size(img)/2);
%img = (img-min(img(:)))/(max(img(:))-min(img(:)))*255;
figure
montage({squeeze(img(im_center(1),:,:)),squeeze(img(:,im_center(2),:)),squeeze(img(:,:,im_center(3)))},gray, 'Size',[1 3])
end