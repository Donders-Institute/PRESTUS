# developed by Eleonora Carpino
# eleonora.carpino@icloud.com
# 14 June 2023

# Code to create a pseudoCT, starting from a T1w image and a PETRA UTE image

# The first step is to run a segmentation using SimNIBS (charm) where the T1 is a normal weighted T1 and instead of the T2 we use a PETRA UTE scan, so that the UTE will be registered to the T1

read -p "Enter the subject number with three digits (e.g., 001 or 022): " subject_number

# Modify the directory so that it is the directory where the segmentation output is stored (SimNIBS m2m folder)
cd "/home/action/elecar/MyProject/new_segmentations/m2m_sub-${subject_number}/"
mv T2_reg.nii.gz UTE_reg.nii.gz # change name to the T2_reg

# Create a skull mask that will be used later
fslmaths final_tissues -thr 7 -uthr 8 -bin final_tissues_skull

# Thresholding at 0
fslmaths UTE_reg -thr 0 UTE_reg_thr0

# Bias Field correction
# For this step ANTs needs to be installed. Modify the directory so that it contains the ANTs folder where the N4BiasFieldCorrection function is stored
/opt/ANTs/20180607/build/bin/N4BiasFieldCorrection \
--image-dimensionality 3 \
--input-image UTE_reg_thr0.nii.gz \
--convergence [50x50x50x50,0.0000001] \
--bspline-fitting [180] \
--output [UTE_reg_thr0_corr.nii.gz,UTE_thr0_BiasField.nii.gz]

# Normalisation
# Normalise dividing for the soft tissue peak intensity value
# Create a txt file to save the soft tissue peak intensity
touch soft_tissue_value.txt

# Use the Matlab code 'soft_tissue_peak.m' to plot the histogram distribution of the intensities and identify the soft tissue peak and then run the second part of the script

