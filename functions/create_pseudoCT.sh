#!/bin/bash

# Code to create a pseudoCT, starting from a T1w image and a PETRA UTE image

# Example usage:
# create_pseudoCT "001" "${rootpath}/data/simnibs"

# developed by Eleonora Carpino
# eleonora.carpino@icloud.com
# 14 June 2023

function create_pseudoCT() 
{
    echo "starting pseudoCT creation"

    local sub_id="$1"       # Subject ID (without 'sub-')
    local m2m_path="$2"     # Path to the SimNIBS folder
    local scriptpath="$3"   # Path to the MATLAB function find_soft_tissue_peak.m

    # Display subject ID and SimNIBS path
    echo "Subject ID: $sub_id"
    echo "M2M Path: $m2m_path"

    # The first step is to run a segmentation using SimNIBS (charm) where the T1 is a normal weighted T1 
    # and instead of the T2 we use a PETRA UTE scan, so that the UTE will be registered to the T1
    # Change directory to the SimNIBS segmentation folder for the subject
    cd "${m2m_path}/m2m_sub-${sub_id}/" || { echo "Directory not found"; exit 1; }

    # Check if UTE_reg.nii.gz already exists. If it exists, skip renaming the file.
    if [ -f "UTE_reg.nii.gz" ]; then
        echo "UTE_reg.nii.gz already exists. Skipping renaming."
    else
        # Rename T2_reg.nii.gz to UTE_reg.nii.gz to use the PETRA UTE image instead of T2
        mv T2_reg.nii.gz UTE_reg.nii.gz
        echo "Renamed T2_reg.nii.gz to UTE_reg.nii.gz"
    fi

    # Create a binary skull mask from final_tissues, thresholding to include skull (value 7 or 8)
    fslmaths final_tissues -thr 7 -uthr 8 -bin final_tissues_skull

    # Apply a threshold to the UTE image at 0, eliminating negative values
    fslmaths UTE_reg -thr 0 UTE_reg_thr0

    # Perform bias field correction using ANTs' N4BiasFieldCorrection to remove intensity inhomogeneity
    # This step requires ANTs to be installed and available in the environment (e.g., module load ants)
    if [ -f "UTE_reg_thr0_corr.nii.gz" ]; then
        echo "N4BiasFieldCorrection already run; reusing estimates."
    else
        echo "Running N4BiasFieldCorrection."
        N4BiasFieldCorrection \
        --image-dimensionality 3 \
        --input-image UTE_reg_thr0.nii.gz \
        --convergence [50x50x50x50,0.0000001] \
        --bspline-fitting [180] \
        --output [UTE_reg_thr0_corr.nii.gz,UTE_thr0_BiasField.nii.gz]
    fi

    # Create a file to store the soft tissue peak intensity value (for later normalization)
    touch soft_tissue_value.txt

    # Call MATLAB to run the soft_tissue_peak.m script to find the peak intensity of soft tissue
    echo "Running MATLAB script soft_tissue_peak.m..."
    cd "${scriptpath}" || { echo "Directory not found"; exit 1; }
    matlab_command="find_soft_tissue_peak(\"$sub_id\",\"$m2m_path\")"
    matlab -nodisplay -batch "$matlab_command" || { echo "MATLAB function failed"; exit 1; }

    # Check exit status of MATLAB command
    # If MATLAB script completes successfully, proceed to normalization
    if [ $? -eq 0 ]; then
        echo "MATLAB function completed successfully."
        # Return to SimNIBS m2m folder (if not there already)
        cd "${m2m_path}/m2m_sub-${sub_id}/" || { echo "Directory not found"; exit 1; }

        # Read the soft tissue peak value from the file created by the MATLAB script
        read peak_value < soft_tissue_value.txt
        echo "The soft tissue peak value obtained from the Matlab code is: $peak_value"

        # Normalize UTE image using the soft tissue peak value -> soft tissue peak = 1 
        fslmaths UTE_reg_thr0_corr -div "$peak_value" UTE_reg_thr0_corr_norm

        # Linear mapping: map image intensities to Hounsfield units (HU): pCT = -2194 UTE + 2236
        fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull pCT_skull

        # Correct partial volume effects in different tissue regions
        fslmaths final_tissues_skull -dilF final_tissues_skull_dil
        fslmaths final_tissues -bin -sub final_tissues_skull_dil -mul final_tissues final_tissues_dil
        fslmaths final_tissues_skull_dil -mul 7 -add final_tissues_dil final_tissues_dil
        fslmaths final_tissues_dil -bin -dilF final_tissues_dil_d
        fslmaths final_tissues_dil -bin -ero final_tissues_dil_e
        fslmaths final_tissues_dil_d -sub final_tissues_dil_e outer_head
        fslmaths final_tissues_dil -bin -sub outer_head -mul final_tissues_dil final_tissues_dil
        fslmaths outer_head -mul 11 -add final_tissues_dil final_tissues_dil
        fslmaths final_tissues_dil -thr 9 -uthr 10 -bin final_tissues_blood_muscle_dil
        fslmaths final_tissues_dil -thr 4 -uthr 5 -bin final_tissues_skin_dil
        fslmaths final_tissues_dil -thr 1 -uthr 3 -add final_tissues_blood_muscle_dil -add final_tissues_skin_dil -bin final_tissues_soft_tissue_dil
        fslmaths final_tissues_dil -add 1 -uthr 1 -bin final_tissues_air_dil
        fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 3 -uthr 4 -bin final_tissues_skull_dil
        fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 1 -uthr 1 outer_head
        
        # Mapping: set air fraction to 1000 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
        fslmaths final_tissues_air_dil -mul -1000 pCT_air_PV 
        # Mapping: set soft-tissue fraction to 42 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
        fslmaths final_tissues_soft_tissue_dil -mul 42 pCT_soft_tissues_PV 

        # These steps may not yet be documented?
        fslmaths outer_head -mul UTE_reg_thr0_corr_norm -mul 1389.333 -sub 1000 -mul outer_head outer_head_PV
        fslmaths outer_head_PV -thr 42 -bin -mul 42 outer_head1_PV
        fslmaths outer_head_PV -uthr 42 outer_head2_PV
        # Mapping: set dilated skull to pCT = -2194 UTE + 2236
        fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull_dil pCT_skull_PV 
        fslmaths outer_head1_PV -add outer_head2_PV -add pCT_skull_PV -add pCT_air_PV -add pCT_soft_tissues_PV pseudoCT_PV
        
        # Smooting: (kernel size = 3x3x3 voxels) at the skin-air and skull-air interfaces.
        SmoothImage 3 pseudoCT_PV.nii.gz 0.8 pseudoCT_PV_smooth.nii.gz 
        
        # Final processing to combine all masks and images into the final pseudoCT
        fslmaths final_tissues_skull -add 1 -thr 1 -uthr 1 -mul pseudoCT_PV_smooth -add pCT_skull pseudoCT

        # Create and store a final tissue mask for simulations
        fslmaths final_tissues_dil -thr 1 -uthr 2 -bin -mul 2 brain_mask
        fslmaths final_tissues_dil -thr 10 -uthr 10 -bin -add final_tissues_skin_dil -add outer_head -bin -mul 5 skin_mask
        fslmaths final_tissues_skull_dil -bin -mul 7 skull_mask
        fslmaths final_tissues_dil -thr 3 -uthr 3 -bin -mul 3 csf_mask
        fslmaths brain_mask -add skin_mask -add skull_mask -add csf_mask tissues_mask

        # Move intermediate and temporary files to a backup directory and clean up
        mkdir pseudoCT_backup
        mv final_tissues_skull.nii.gz pseudoCT_backup
        mv UTE_reg_thr0.nii.gz pseudoCT_backup
        mv UTE_reg_thr0_corr.nii.gz pseudoCT_backup
        mv UTE_thr0_BiasField.nii.gz pseudoCT_backup
        mv final_tissues_dil.nii.gz pseudoCT_backup
        mv final_tissues_skull_dil.nii.gz pseudoCT_backup
        mv pCT_skull.nii.gz pseudoCT_backup
        mv pCT_skull_PV.nii.gz pseudoCT_backup
        mv pseudoCT_PV.nii.gz pseudoCT_backup
        mv pseudoCT_PV_smooth.nii.gz pseudoCT_backup
        mv skull_mask.nii.gz pseudoCT_backup
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
        echo "Process complete. You can now run the simulation with the pseudoCT. Ensure both 'pseudoCT.nii.gz' and 'tissues_mask.nii.gz' are in the m2m segmentation subject folder."
    else
        echo "MATLAB function failed"
        exit 1
    fi
    }