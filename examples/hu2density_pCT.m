
% Convert the HU of the pseudoCT into density values, using the
% hounsfield2density k-Wave function

cd /home/action/elecar/PRESTUS
addpath(genpath('toolboxes')) 

% Load pseudoCT
filename_pseudoCT = fullfile('/project/3023001.06/eleonora/new_segmentations/m2m_sub-022/pseudoCT14062023.nii.gz');
ct_data = niftiread(filename_pseudoCT);
pseudoCT_header = niftiinfo(filename_pseudoCT); 

filename_tissues_mask = fullfile('/project/3023001.06/eleonora/new_segmentations/m2m_sub-022/tissues_mask14062023.nii.gz');
tissues_mask = niftiread(filename_tissues_mask); 
tissues_mask_header = niftiinfo(filename_tissues_mask); 

 
% Add 1000 to the pCT data because the calibration in the hounsfield2density
% function requires an offset of 1000.
ct_data(tissues_mask == 4) = ct_data(tissues_mask == 4) + 1000;
% Mask of the positive values
ct_data(tissues_mask == 4) = max(ct_data(tissues_mask == 4),0);

% Create an empty initial density matrix
density = zeros(size(ct_data));

% Apply the function
density(tissues_mask == 4) = hounsfield2density(ct_data(tissues_mask == 4), true);

% Make histograms of the results
%histogram(ct_data_offset(mask))
%figure
%histogram(density(mask))
% It seems that it works because the histograms are similar but different
% and they seem to resemble the paper equation

density = max(density,0);

filename_ct_offset = fullfile('/project/3023001.06/eleonora/new_segmentations/m2m_sub-022/ct_offset');
niftiwrite(ct_data, filename_ct_offset, pseudoCT_header, 'Compressed',true);

density_header = pseudoCT_header;
density_header.Datatype = 'double';
density_header.BitsPerPixel = 64;
filename_density = fullfile('/project/3023001.06/eleonora/new_segmentations/m2m_sub-022/density');
niftiwrite(density, filename_density, density_header, 'Compressed',true);

niftiwrite(density, filename_density, 'Compressed',true);

