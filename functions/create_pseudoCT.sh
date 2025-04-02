#!/bin/bash

# Create a pseudoCT from a PETRA UTE image and a SimNIBS 4 segmentation

# Example usage:
# create_pseudoCT "001" "${rootpath}/data/simnibs" "/opt/matlab/R2022b/bin/matlab" "kosciessa" "1"

function create_pseudoCT() 
{
    echo "starting pseudoCT creation"

    local sub_id="$1"        # Subject ID (without 'sub-')
    local path_simnibs="$2"  # Path to SimNIBS folder
    local path_matlab="$3"   # Path to MATLAB executable (e.g., "/opt/matlab/R2022b/bin/matlab")
    local skullmapping=${4:-"kosciessa"}  # Algorithm for linear skull mapping ["miscouridou"; "carpino"; "wiesinger"; "treeby"; "kosciessa" (default)]
    local debug=${5:-"1"}    # Save intermediate images [yes: "1" (default); no: "0"]
    local path_script="${6:-$(dirname "$(realpath "$0")")}"   # Path to PRESTUS/functions; specify when used in HPC job [default: Use path to this function]

    # Display subject ID and SimNIBS path
    echo "ID: $sub_id"
    echo "SimNIBS path: $path_simnibs"
    echo "MATLAB binary: $path_matlab"
    echo "Function path: $path_script"

    # The first step is to run a segmentation using SimNIBS (charm) where the T1 is a normal weighted T1 
    # and instead of the T2 we use a PETRA UTE scan, so that the UTE will be registered to the T1
    # Change directory to the SimNIBS segmentation folder for the subject
    cd "${path_simnibs}/m2m_sub-${sub_id}/" || { echo "Directory not found"; exit 1; }

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

    # create dedicated pCT directory for intermediate backups
    if [ ! -d "pseudoCT" ]; then
        mkdir "pseudoCT"
    fi
    path_pct="${path_simnibs}/m2m_sub-${sub_id}/pseudoCT"

    # Create a file to store the soft tissue peak intensity value (for later normalization)
    touch ${path_pct}/pCT_soft_tissue_value.txt

    # Call MATLAB to run the soft_tissue_peak.m script to find the peak intensity of soft tissue
    echo "Running MATLAB script pct_soft_tissue_peak.m..."
    cd "${path_script}" || { echo "Directory not found"; exit 1; }
    matlab_command="pct_soft_tissue_peak(\"$sub_id\",\"$path_simnibs\")"
    $path_matlab -nodisplay -r "$matlab_command" || { echo "MATLAB function failed"; exit 1; }

    # Check exit status of MATLAB command
    # If MATLAB script completes successfully, proceed to normalization
    if [ $? -eq 0 ]; then
        echo "MATLAB function completed successfully."
        # Move to simnibs folder
        cd "${path_simnibs}/m2m_sub-${sub_id}" || { echo "Directory not found"; exit 1; }

        # Read the soft tissue peak value from the file created by the MATLAB script
        read peak_value < ${path_pct}/pCT_soft_tissue_value.txt
        echo "The soft tissue peak value obtained from the Matlab code is: $peak_value"

        # Normalize UTE image using the soft tissue peak value -> soft tissue peak = 1 
        fslmaths UTE_reg_thr0_corr -div "$peak_value" UTE_reg_thr0_corr_norm

        # The skull fraction should now be between the noise/air peak (around 0) and the soft tissue peak (1).

        # Create a file to store the skull mapping
        touch ${path_pct}/pCT_skull_mapping.txt
        # Linear mapping: map image intensities in skull to Hounsfield units (HU)
        case "$skullmapping" in
            "miscouridou")
                echo "Mapping skull ~ miscouridou... pHU = −2085 UTE + 2329 "
                fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul -2085 -add 2329 -mul final_tissues_skull pCT_skull
                ;;
            "carpino")
                echo "Mapping skull ~ carpino... pHU = -2194 UTE + 2236"
                fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull pCT_skull
                ;;
            "wiesinger")
                echo "Mapping skull ~ wiesinger... pHU = -2000 (UTE-1) + 42"
                fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -add -1 -mul -2000 -add 42 -mul final_tissues_skull pCT_skull
                ;;
            "treeby")
                echo "Mapping skull ~ treeby... pHU = -2929.6 UTE + 3247.9"
                fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul -2929.6 -add 3247.9 -mul final_tissues_skull pCT_skull
                ;;
            "kosciessa")
                echo "Mapping skull ~ kosciessa (default)... establishing linear skull mapping"
                # Call MATLAB to run the pct_skullmapping.m script
                echo "Running MATLAB script pct_skullmapping.m..."
                cd "${path_script}" || { echo "Directory not found"; exit 1; }
                matlab_command="pct_skullmapping(\"$sub_id\",\"$path_simnibs\")"
                $path_matlab -nodisplay -r "$matlab_command" || { echo "MATLAB function failed"; exit 1; }
                # Move to SimNIBS m2m folder (if not there already)
                cd "${path_simnibs}/m2m_sub-${sub_id}/" || { echo "Directory not found"; exit 1; }
                # Read the mapping coefficients
                if [ -f ${path_pct}/pCT_skull_mapping.txt ]; then
                    read m1 <<< $(awk 'NR==1 {print $1}' ${path_pct}/pCT_skull_mapping.txt)
                    read m2 <<< $(awk 'NR==2 {print $1}' ${path_pct}/pCT_skull_mapping.txt)
                else
                    echo "Error: File pCT_skull_mapping.txt does not exist."
                    exit 1
                fi
                echo ""
                echo "Processing kosciessa... pCT = $m1 UTE + $m2"
                fslmaths final_tissues_skull -mul UTE_reg_thr0_corr_norm -mul $m1 -add $m2 -mul final_tissues_skull pCT_skull
                ;;
        esac
        # Record which mapping variant has been used
        echo "$skullmapping" >> ${path_pct}/pCT_skull_mapping.txt

        # Correct partial volume effects in tissues
        # Warning: some hard-coded labels (SimNibs 4 tissue labels)

        # Dilate Skull Mask
        fslmaths final_tissues_skull -dilF final_tissues_skull_dil
        # Refine Soft Tissue Mask (subtract dilated skull)
        fslmaths final_tissues -bin -sub final_tissues_skull_dil -mul final_tissues final_tissues_dil
        # Label (7) [Dilated Skull]
        fslmaths final_tissues_skull_dil -mul 7 -add final_tissues_dil final_tissues_dil
        # Dilate All Tissues
        fslmaths final_tissues_dil -bin -dilF final_tissues_dil_d
        # Erode All Tissues
        fslmaths final_tissues_dil -bin -ero final_tissues_dil_e
        # Create Outer Head Mask (dilated minus eroded)
        fslmaths final_tissues_dil_d -sub final_tissues_dil_e outer_head
        # Refine Tissue Mask by Removing Outer Head
        fslmaths final_tissues_dil -bin -sub outer_head -mul final_tissues_dil final_tissues_dil
        # Label (11) [Outer Head]
        fslmaths outer_head -mul 11 -add final_tissues_dil final_tissues_dil
        # Create Blood/Muscle Mask
        fslmaths final_tissues_dil -thr 9 -uthr 10 -bin final_tissues_blood_muscle_dil
        # Create Skin Mask
        fslmaths final_tissues_dil -thr 4 -uthr 5 -bin final_tissues_skin_dil
        # Create Soft Tissue Mask
        fslmaths final_tissues_dil -thr 1 -uthr 3 -add final_tissues_blood_muscle_dil -add final_tissues_skin_dil -bin final_tissues_soft_tissue_dil
        # Create Air Mask
        fslmaths final_tissues_dil -add 1 -uthr 1 -bin final_tissues_air_dil
        # Refine Skull Mask [assign 3 to dilated skull, + outer head (dilated minus eroded)); threshold between 3 (skull) and 4(skull+dilated)]
        fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 3 -uthr 4 -bin final_tissues_skull_dil
        # Refine Outer Head Mask [assign 3 to dilated skull, + outer head (dilated minus eroded)); select only label 1]
        fslmaths final_tissues_skull_dil -mul 3 -add outer_head -thr 1 -uthr 1 outer_head
        
        # Mapping: Air Fraction = -1000 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
        fslmaths final_tissues_air_dil -mul -1000 pCT_air_PV 
        # Mapping: Soft-Tissue Fraction = 42 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
        fslmaths final_tissues_soft_tissue_dil -mul 42 pCT_soft_tissues_PV 
        # Mapping: Dilated Skull Fraction
        case "$skullmapping" in 
            "miscouridou")
                fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul -2085 -add 2329 -mul final_tissues_skull_dil pCT_skull_PV 
                ;;
            "carpino")
                fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul -2194 -add 2236 -mul final_tissues_skull_dil pCT_skull_PV 
                ;;
            "wiesinger")
                fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -add -1 -mul -2000 -add 42 -mul final_tissues_skull_dil pCT_skull_PV 
                ;;
            "treeby")
                fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul -2929.6 -add 3247.9 -mul final_tissues_skull_dil pCT_skull_PV 
                ;;
            "kosciessa")
                fslmaths final_tissues_skull_dil -mul UTE_reg_thr0_corr_norm -mul $m1 -add $m2 -mul final_tissues_skull_dil pCT_skull_PV 
                ;;
        esac
        
        # Mapping: Outer Head Fraction: scale UTE_norm by 1389.333 and subtract 1000 (not documented)
        fslmaths outer_head -mul UTE_reg_thr0_corr_norm -mul 1389.333 -sub 1000 -mul outer_head outer_head_PV
        # Thresholds Outer Head Fraction to pHU ≥42
        fslmaths outer_head_PV -thr 42 -bin -mul 42 outer_head1_PV
        # Thresholds Outer Head Fraction to pHU ≤42
        fslmaths outer_head_PV -uthr 42 outer_head2_PV

        # Combine Tissue Maps into Final PseudoCT
        fslmaths outer_head1_PV -add outer_head2_PV -add pCT_skull_PV -add pCT_air_PV -add pCT_soft_tissues_PV pCT_PV
    
        # Smooting: kernel size = 3x3x3 voxels; smoothing factor = 0.8
        SmoothImage 3 pCT_PV.nii.gz 0.8 pCT_PV_smooth.nii.gz 
        
        # Final processing to combine all masks and images into the final pseudoCT
        fslmaths final_tissues_skull -add 1 -thr 1 -uthr 1 -mul pCT_PV_smooth -add pCT_skull pseudoCT

        # Create and store a final tissue mask for simulations
        fslmaths final_tissues_dil -thr 1 -uthr 2 -bin -mul 2 brain_mask
        fslmaths final_tissues_dil -thr 10 -uthr 10 -bin -add final_tissues_skin_dil -add outer_head -bin -mul 5 skin_mask
        fslmaths final_tissues_skull_dil -bin -mul 7 skull_mask
        fslmaths final_tissues_dil -thr 3 -uthr 3 -bin -mul 3 csf_mask
        fslmaths brain_mask -add skin_mask -add skull_mask -add csf_mask tissues_mask

        # Copy pCT to backup directory
        cp pseudoCT.nii.gz ${path_pct}

        # Move intermediate files to backup directory
        mv final_tissues_skull.nii.gz ${path_pct}
        mv UTE_reg_thr0.nii.gz ${path_pct}
        mv UTE_reg_thr0_corr.nii.gz ${path_pct}
        mv UTE_thr0_BiasField.nii.gz ${path_pct}
        mv final_tissues_dil.nii.gz ${path_pct}
        mv final_tissues_skull_dil.nii.gz ${path_pct}
        mv pCT_skull.nii.gz ${path_pct}
        mv pCT_skull_PV.nii.gz ${path_pct}
        mv pCT_PV.nii.gz ${path_pct}
        mv pCT_PV_smooth.nii.gz ${path_pct}
        mv skull_mask.nii.gz ${path_pct}
        mv final_tissues_dil_d.nii.gz ${path_pct}
        mv final_tissues_dil_e.nii.gz ${path_pct}
        mv outer_head.nii.gz ${path_pct}
        mv final_tissues_blood_muscle_dil.nii.gz ${path_pct}
        mv final_tissues_soft_tissue_dil.nii.gz ${path_pct}
        mv final_tissues_air_dil.nii.gz ${path_pct}
        mv final_tissues_skin_dil.nii.gz ${path_pct}
        mv pCT_air_PV.nii.gz ${path_pct}
        mv pCT_soft_tissues_PV.nii.gz ${path_pct}
        mv outer_head_PV.nii.gz ${path_pct}
        mv outer_head1_PV.nii.gz ${path_pct}
        mv outer_head2_PV.nii.gz ${path_pct}
        mv brain_mask.nii.gz ${path_pct}
        mv skin_mask.nii.gz ${path_pct}
        mv csf_mask.nii.gz ${path_pct}

        # If DEBUG=0, remove nifti images in pseudoCT folder
        if [ "$debug" -eq 0 ]; then
            echo "Cleanup: Removing nifti images in pseudoCT folder"
            rm ${path_pct}/*.nii.gz
        fi

        echo "Process complete. You can now run the simulation with the pseudoCT. Ensure both 'pseudoCT.nii.gz' and 'tissues_mask.nii.gz' are in the m2m segmentation subject folder."
    else
        echo "MATLAB function failed"
        exit 1
    fi
    }