function plot_coronal_slices(img, coord, options)
arguments
    img (:,:,:)
    coord (3,1) = round(size(img)/2)
    options.cmap = []
    options.rescale = 0
    options.labels = ["Background","White Matter", "Gray Matter", "CSF", "Bone", "Scalp", "Eye balls", "Compact bone", "Spongy bone", "Blood", "Muscle", "Electrode", "Saline or gel"] % defaults based on SimNIBS 4
end
default_colors = [230,230,230;129,129,129;104,163,255;255,239,179;255,166,133;255,240,0;255,239,179;255,138,57;0,65,142;0,118,14;37,79,255;103,255,226;];
if isinteger(img) && isempty(options.cmap)
    N = length(unique(img(:)));
    options.cmap = [[0,0,0];default_colors(1:N,:)/255];
end
if options.rescale
img= round(rescale(img, 1, 255));
end
figure;
montage({imrotate(squeeze(img(coord(1),:,:)),90),imrotate(squeeze(img(:,coord(2),:)),90),squeeze(img(:,:,coord(3)))'}, options.cmap, 'Size',[1 3])

if isinteger(img)
hold on
L = line(ones(N),ones(N), 'LineWidth',2);               % generate line 
set(L,{'color'},mat2cell(options.cmap(1:N,:),ones(1,N),3));            % set the colors according to cmap
legend(options.labels(unique(img(:))+1))                                 % add as many legend entries as data
end
end