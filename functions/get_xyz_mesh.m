function coord_mesh_xyz = get_xyz_mesh(img)
%% coord_mesh_xyz 
% create an Nx3 array of point coordinates for a 3D image
[x, y, z] = ndgrid(1:size(img, 1),1:size(img, 2),1:size(img, 3));

coord_mesh_xyz = [reshape(x,[],1) reshape(y,[],1) reshape(z,[],1)];

end