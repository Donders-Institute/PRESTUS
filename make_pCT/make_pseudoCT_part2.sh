# developed by Eleonora Carpino
# eleonora.carpino@icloud.com
# 14 June 2023

# Code to create a pseudoCT, second part

# Normalisation
# Normalise dividing for the soft tissue peak intensity value
# Use the Matlab code 'soft_tissue_peak.m' to plot the histogram distribution of the intensities and identify the soft tissue peak. The value is automatically stored in the 'soft_tissue_value.txt' file

# Divide for the soft tissue peak intensity
read -p "Enter the subject number with three digits (e.g., 001 or 022): " subject_number
cd "/home/action/elecar/MyProject/new_segmentations/m2m_sub-${subject_number}/"

# Automatically read the soft tissue peak value
read peak_value < soft_tissue_value.txt
echo "The soft tissue peak value obtained from the Matlab code is: $peak_value"
fslmaths UTE_reg_thr0_corr -div "$peak_value" UTE_reg_thr0_corr_norm

# Linear mapping
# Mapping from intensities to HU
# Map in the skull with the equation: pCT = -2194 UTE + 2236
fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull pCT_skull

# Partial volume effects correction
# For the skull/air and the skin/air interface 
# Modify the final tissues mask with a dilated skull and a contour of the outer head
# Dilating the skull mask of 3x3x3 voxels volume
fslmaths final_tissues_skull -dilF final_tissues_skull_dil
# Adding it to the final tissues mask so that it contains the dilated skull
fslmaths final_tissues -bin -sub final_tissues_skull_dil -mul final_tissues final_tissues_dil
fslmaths final_tissues_skull_dil -mul 7 -add final_tissues_dil final_tissues_dil
#fslmaths final_tissues_skull_dil -mul 10 -add final_tissues final_tissues_dil
# Create a contour of the head for the skin/air interface
fslmaths final_tissues_dil -bin -dilF final_tissues_dil_d
fslmaths final_tissues_dil -bin -ero final_tissues_dil_e
fslmaths final_tissues_dil_d -sub final_tissues_dil_e outer_head
# Adding it to the final tissues mask so that it also contains the head contour
fslmaths final_tissues_dil -bin -sub outer_head -mul final_tissues_dil final_tissues_dil
fslmaths outer_head -mul 11 -add final_tissues_dil final_tissues_dil
#fslmaths outer_head -mul 50 -add final_tissues_dil final_tissues_dil
# Create a skull/air interface? To isolate the sinuses with respect to the rest of the outer head

# Create the new tissues masks in the dilated mask
fslmaths final_tissues_dil -thr 9 -uthr 10 -bin final_tissues_blood_muscle_dil
fslmaths final_tissues_dil -thr 4 -uthr 5 -bin final_tissues_skin_dil #it would also contain a few voxels of other types of bone but they are not usually present
fslmaths final_tissues_dil -thr 1 -uthr 3 -add final_tissues_blood_muscle_dil -add final_tissues_skin_dil -bin final_tissues_soft_tissue_dil
fslmaths final_tissues_dil -add 1 -uthr 1 -bin final_tissues_air_dil

# Separate the mask for the dilated skull and the mask for the outer head attributing the overlapping parts to the skull
fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 3 -uthr 4 -bin final_tissues_skull_dil
fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 1 -uthr 1 outer_head

# Apply the mapping equations
# Mapping from intensities to HU
# Give the value -1000 HU to air 
fslmaths final_tissues_air_dil -mul -1000 pCT_air_PV 
# Give the value 42 HU to soft tissue
fslmaths final_tissues_soft_tissue_dil -mul 42 pCT_soft_tissues_PV 
# Equation for the skin/air interface (outer head):  
fslmaths outer_head -mul UTE_reg_thr0_corr_norm -mul 1389.333 -sub 1000 -mul outer_head outer_head_PV
# Thresholding on the outer head to correct values below 42 (soft tissue threshold) and make them become 42
fslmaths outer_head_PV -thr 42 -bin -mul 42 outer_head1_PV
fslmaths outer_head_PV -uthr 42 outer_head2_PV
# Equation for the skull: pCT = -2194 UTE + 2236
fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull_dil pCT_skull_PV 
# Add up the different parts to create the initial PVE corrected pCT
fslmaths outer_head1_PV -add outer_head2_PV -add pCT_skull_PV -add pCT_air_PV -add pCT_soft_tissues_PV pseudoCT_PV

# General smoothing 
/opt/ANTs/20180607/build/bin/SmoothImage 3 pseudoCT_PV.nii.gz 0.8 pseudoCT_PV_smooth.nii.gz 

# Re-apply the initial values in the skull part not to loose resolution
fslmaths final_tissues_skull -add 1 -thr 1 -uthr 1 -mul pseudoCT_PV_smooth -add pCT_skull pseudoCT  # this is the final pseudoCT

# Final tissues mask
# To run the simulations with PRESTUS it is necessary to adapt the pseudoCT and create a final 'tissues_mask.nii'  file
fslmaths final_tissues_dil -thr 1 -uthr 2 -bin -mul 2 brain_mask
fslmaths final_tissues_dil -thr 10 -uthr 10 -bin -add final_tissues_skin_dil -add outer_head -bin -mul 3 skin_mask #it also contains some small muscle parts and a few voxels of other bone structures (very rarely)
# we consider the whole outer head as skin (skin is also scalp anyway)
#fslmaths final_tissues_dil -thr 7 -uthr 7 -bin -mul 4 skull_cortical_mask
#fslmaths final_tissues_dil -thr 8 -uthr 8 -bin -mul 5 skull_trabecular_mask
fslmaths final_tissues_skull_dil -bin -mul 4 skull_mask
# the skull doesn't have a distinction between cortical and trabecular anymore, but we consider a unique mask with value 4
fslmaths final_tissues_dil -thr 3 -uthr 3 -bin -mul 6 csf_mask #without air csf needs to be 6

fslmaths brain_mask -add skin_mask -add skull_mask -add csf_mask tissues_mask # this should have value 0 for the places where there is air

# Move the files obtained during the process in a sub-folder. Eliminate other files that are not useful. Everything in the sub-folder can be eliminated.
mkdir process_backup

rm final_tissues_dil_d.nii.gz
rm final_tissues_dil_e.nii.gz
rm outer_head.nii.gz
rm final_tissues_blood_muscle_dil.nii.gz
rm final_tissues_soft_tissue_dil.nii.gz
rm final_tissues_air_dil.nii.gz
rm final_tissues_skin_dil.nii.gz
rm pCT_air_PV.nii.gz
rm pCT_soft_tissues_PV.nii.gz
rm outer_head_PV.nii.gz
rm outer_head1_PV.nii.gz
rm outer_head2_PV.nii.gz
rm brain_mask.nii.gz
rm skin_mask.nii.gz
rm csf_mask.nii.gz

mv final_tissues_skull.nii.gz process_backup
mv UTE_reg_thr0.nii.gz process_backup
mv UTE_reg_thr0_corr.nii.gz process_backup
mv UTE_thr0_BiasField.nii.gz process_backup
mv final_tissues_dil.nii.gz process_backup
mv final_tissues_skull_dil.nii.gz process_backup
mv pCT_skull.nii.gz process_backup
mv pCT_skull_PV.nii.gz process_backup
mv pseudoCT_PV.nii.gz process_backup
mv pseudoCT_PV_smooth.nii.gz process_backup
mv skull_mask.nii.gz process_backup


# You can now run the simulation with the pseudoCT. It is important that both 'pseudoCT.nii.gz' and 'tissues_mask.nii.gz' are in the m2m segmentation subject folder
