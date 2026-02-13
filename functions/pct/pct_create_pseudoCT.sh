#!/bin/bash

# Create a pseudoCT from a PETRA UTE image and a SimNIBS 4 segmentation

# Example usage:
# pct_create_pseudoCT "001" "${rootpath}/data/simnibs" "/opt/matlab/R2022b/bin/matlab" "kosciessa" "1"

function pct_create_pseudoCT() 
{
    echo "starting pseudoCT creation"

    local sub_id="$1"        # Subject ID (without 'sub-')
    local path_simnibs="$2"  # Path to SimNIBS folder
    local path_matlab="$3"   # Path to MATLAB executable (e.g., "/opt/matlab/R2022b/bin/matlab")
    local skullmapping=${4:-"kosciessa"}  # Algorithm for linear skull mapping ["miscouridou"; "carpino"; "wiesinger"; "treeby"; "kosciessa" (default)]
    local debug=${5:-"1"}    # Save intermediate images [yes: "1" (default); no: "0"]
    local path_fun="${6:-$(dirname "$(realpath "$0")")}" # Path to PRESTUS/functions; specify when used in HPC job [default: Use path to this function]
    local path_fun_pct="${path_fun}/pct"
    local path_fun_head="${path_fun}/head"

    # shift simnibs path to sub-id
    path_simnibs="${path_simnibs}/m2m_sub-${sub_id}"

    # create dedicated pCT directory for intermediate backups under the SimNIBS path
    path_pct="${path_simnibs}/pseudoCT"
    if [ ! -d $path_pct ]; then
        mkdir "$path_pct"
    fi

    # Create a log
    exec > >(tee -a "${path_pct}/pct_create_${sub_id}_$(date +%Y%m%d_%H%M%S).log") 2>&1

    echo "=== PCT Creation Log ==="
    echo "Date: $(date)"
    echo "Subject: $sub_id"
    echo "Log dir: $log_date_dir"
    echo "Starting pseudoCT creation"
    echo "SimNIBS path: $path_simnibs"
    echo "pCT path: $path_pct"
    echo "MATLAB binary: $path_matlab"
    echo "Function path: $path_fun"
    echo "Function PCT path: $path_fun_pct"
    echo "Function Head path: $path_fun_head"

    # Key file names
    segmentation="${path_simnibs}/final_tissues.nii.gz"
    skull_mask="${path_pct}/skull_mask.nii.gz"
    pCT="${path_pct}/pseudoCT.nii.gz"

    # Temporary file names
    t2w_reg="${path_simnibs}/T2_reg.nii.gz"
    ute_reg="${path_simnibs}/UTE_reg.nii.gz"
    ute_reg_thr0="${path_pct}/UTE_reg_thr0.nii.gz"
    ute_reg_thr0_corr="${path_pct}/UTE_reg_thr0_corr.nii.gz"
    ute_reg_bias="${path_pct}/UTE_thr0_BiasField.nii.gz"
    ute_stnorm="${path_pct}/UTE_STnorm.nii.gz"

    # Externaly run a segmentation using SimNIBS (charm) where the T1 is a normal weighted T1 
    # and instead of the T2 we use a PETRA UTE scan, so that the UTE will be registered to the T1

    # Change directory to the pCT folder
    cd "${path_pct}" || { echo "Directory not found"; exit 1; }

    # Check if UTE_reg.nii.gz already exists. If it exists, skip renaming the file.
    if [ -f $ute_reg ]; then
        echo "UTE_reg.nii.gz already exists. Skipping renaming."
    else
        # Rename T2_reg.nii.gz to UTE_reg.nii.gz to use the PETRA UTE image instead of T2
        mv $t2w_reg $ute_reg
        echo "Renamed T2_reg.nii.gz to UTE_reg.nii.gz"
    fi

    # Create a binary skull mask from final_tissues, thresholding to include skull (value 7 or 8)
    fslmaths $segmentation -thr 7 -uthr 8 -bin $skull_mask

    # Refine skull mask
    echo "Running MATLAB script pct_skullexpand.m..."
    cd "${path_fun_pct}" || { echo "Directory not found: $path_fun_pct"; exit 1; }
    matlab_command="addpath(genpath('${path_fun}')); pct_skullexpand('${path_simnibs}', '${path_pct}'); exit"
    "$path_matlab" -nodisplay -nosplash -batch "$matlab_command" || { echo "MATLAB pct_skullexpand failed"; exit 1; }
    # Move back to pCT folder
    cd "${path_pct}" || { echo "Directory not found"; exit 1; }

    ################################################################################
    ##### pCT processing #####
    ################################################################################

    ################################################################################
    ##### STEP 1: UTE THRESHOLD #####
    ################################################################################

    # Apply a threshold to the UTE image at 0, eliminating negative values
    fslmaths $ute_reg -thr 0 $ute_reg_thr0

    ################################################################################
    ##### STEP 2: UTE BIAS FIELD CORRECTION #####
    ################################################################################

    # Perform bias field correction using ANTs' N4BiasFieldCorrection to remove intensity inhomogeneity
    # This step requires ANTs to be installed and available in the environment (e.g., module load ants)
    if [ -f $ute_reg_thr0_corr ]; then
        echo "N4BiasFieldCorrection already run; reusing estimates."
    else
        echo "Running N4BiasFieldCorrection."
        N4BiasFieldCorrection \
        --image-dimensionality 3 \
        --input-image $ute_reg_thr0 \
        --convergence [50x50x50x50,0.0000001] \
        --bspline-fitting [180] \
        --output [$ute_reg_thr0_corr, $ute_reg_bias]
    fi

    ################################################################################
    ##### STEP 3: UTE SOFT TISSUE NORMALISATION #####
    ################################################################################

    # Create a file to store the soft tissue peak intensity value (for later normalization)
    touch ${path_pct}/pCT_soft_tissue_value.txt

    # Call MATLAB to run the soft_tissue_peak.m script to find the peak intensity of soft tissue
    echo "Running MATLAB script pct_soft_tissue_peak.m..."
    cd "${path_fun_pct}" || { echo "Directory not found"; exit 1; }
    matlab_command="addpath(genpath('${path_fun}')); pct_soft_tissue_peak(\"$path_simnibs\",\"$path_pct\"); exit"
    $path_matlab -nodisplay -r "$matlab_command" || { echo "MATLAB function failed"; exit 1; }
    # Move to pCT folder
    cd "${path_pct}" || { echo "Directory not found"; exit 1; }

    # Read the soft tissue peak value from the file created by the MATLAB script
    read peak_value < ${path_pct}/pCT_soft_tissue_value.txt
    echo "The soft tissue peak value obtained from the Matlab code is: $peak_value"

    # Normalize UTE image using the soft tissue peak value -> soft tissue peak = 1 
    fslmaths $ute_reg_thr0_corr -div "$peak_value" $ute_stnorm

    # The skull fraction should now be between the noise/air peak (around 0) and the soft tissue peak (1).
    
    ################################################################################
    ##### STEP 4: SEGMENTATION PARTIAL VOLUME CORRECTION #####
    ################################################################################

    # Set filenames for intermediate pCT files
    # Masks
    skull_mask_inv="${path_pct}/skull_mask_inv.nii.gz"
    skull_mask_7="${path_pct}/skull_mask_7.nii.gz"
    skull_mask_PVC="${path_pct}/skull_mask_PVC.nii.gz"
    soft_tissues_mask="${path_pct}/soft_tissues_mask.nii.gz"
    air_mask="${path_pct}/air_mask.nii.gz"
    segmentation_no_skull="${path_pct}/segmentation_no_skull.nii.gz"
    segmentation_mask="${path_pct}/segmentation_mask.nii.gz"
    segmentation_mask_no_dil_skull="${path_pct}/segmentation_mask_no_dil_skull.nii.gz"
    segmentation_mask_no_dil_skull_final="${path_pct}/segmentation_mask_no_dil_skull_final.nii.gz"

    # Eroded/Dilated Masks
    skull_mask_dil="${path_pct}/skull_mask_dil.nii.gz"
    soft_tissues_mask_dil="${path_pct}/soft_tissues_mask_dil.nii.gz"
    air_mask_dil="${path_pct}/air_mask_dil.nii.gz"
    segmentation_mask_dil="${path_pct}/segmentation_mask_dil.nii.gz"
    segmentation_mask_ero="${path_pct}/segmentation_mask_ero.nii.gz"
    segmentation_mask_no_dil_skull="${path_pct}/segmentation_mask_no_dil_skull.nii.gz"

    # Skull/Air INTERFACE
    skull_air_interface="${path_pct}/skull_air_interface.nii.gz"
    segmentation_skull_air_PVC="${path_pct}/segmentation_skull_air_PVC.nii.gz"
    segmentation_skull_air_PVC_binv="${path_pct}/segmentation_skull_air_PVC_binv.nii.gz"
    segmentation_skull_air_PVC_skin_air_PVC="${path_pct}/segmentation_skull_air_PVC_skin_air_PVC.nii.gz"
    segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_1to6="${path_pct}/segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_1to6.nii.gz"
    segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_9to13="${path_pct}/segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_9to13.nii.gz"

    # Outer Head
    outer_head_contour="${path_pct}/outer_head_contour.nii.gz"
    outer_head_contour_PVC_inside="${path_pct}/outer_head_contour_PVC_inside.nii.gz"
    outer_head_contour_PVC_outside="${path_pct}/outer_head_contour_PVC_outside.nii.gz"

    # pCT variants
    pCT_air="${path_pct}/pCT_air.nii.gz"
    pCT_soft_tissues="${path_pct}/pCT_soft_tissues.nii.gz"
    pCT_skull="${path_pct}/pCT_skull.nii.gz"
    pCT_PVC="${path_pct}/pCT_PVC.nii.gz"
    pCT_PV_smooth="${path_pct}/pCT_PV_smooth.nii.gz"

    # Partial Volume Correction
    # For the [1] skull/air and the [2] skin/air interface

    # SimNIBS final_tissues.nii.gz labels (from m2m/charm LUT):
    # 0=Air  1-4=WM/GM/CSF  5=Skin  6=Connective/CSF  7-8=Bone(Skull)  9-10=Vessels/CSF 11=L-eye 12=R-eye  13=Optic  >13=Custom

    # [1] SKULL/AIR INTERFACE PVC
    fslmaths $skull_mask -dilF -bin $skull_mask_dil                                    # Dilate skull (catch bone PVE)
    fslmaths $segmentation -bin $segmentation_mask                                     # Binary head mask
    fslmaths $segmentation_mask -sub $skull_mask_dil $segmentation_mask_no_dil_skull   # Head minus dilated skull
    fslmaths $segmentation_mask_no_dil_skull -add 1 $segmentation_mask_no_dil_skull    # Non-skull head → 1+, skull→0
    fslmaths $segmentation_mask_no_dil_skull -binv -mul 7 $skull_air_interface         # Skull PVE voxels=7
    fslmaths $skull_mask -binv $skull_mask_inv                                         # Inverse skull
    fslmaths $segmentation -mul $skull_mask_inv $segmentation_no_skull                 # Segmentation without skull
    fslmaths $skull_mask -mul 7 $skull_mask_7                                          # Pure skull=7
    fslmaths $skull_mask_7 -add $skull_air_interface -add $segmentation_no_skull $segmentation_skull_air_PVC  # Rebuild w/ skull PVE

    # [2] SKIN/AIR INTERFACE PVC
    fslmaths $segmentation_mask -bin -dilF $segmentation_mask_dil                      # Dilate head contour
    fslmaths $segmentation_mask -bin -ero $segmentation_mask_ero                       # Erode head contour  
    fslmaths $segmentation_mask_dil -sub $segmentation_mask_ero $outer_head_contour    # Outer head PVE ring
    fslmaths $segmentation_skull_air_PVC -binv $segmentation_skull_air_PVC_binv        # Non-skullair mask
    fslmaths $outer_head_contour -mul $segmentation_skull_air_PVC_binv $outer_head_contour  # Clean outer ring
    fslmaths $outer_head_contour -mul 5 -add $segmentation_skull_air_PVC $segmentation_skull_air_PVC_skin_air_PVC  # Skin PVE=5

    # === MASKS FOR HU MAPPING ===
    # AIR (0 HU): Everything outside head
    fslmaths $segmentation_skull_air_PVC_skin_air_PVC -binv $air_mask

    # SOFT TISSUES (~42 HU): ALL non-skull/non-air (WM/GM/CSF/Skin/Eyes/Vessels)
    # NOTE: Includes 1-6 (brain/CSF/skin) + 9-13 (eyes/optic/vessels)
    fslmaths $segmentation_skull_air_PVC_skin_air_PVC -thr 1 -uthr 6 -bin $segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_1to6
    fslmaths $segmentation_skull_air_PVC_skin_air_PVC -thr 9 -uthr 13 -bin $segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_9to13
    fslmaths $segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_1to6 -add $segmentation_skull_air_PVC_skin_air_PVC_SoftTissuesMask_9to13 $soft_tissues_mask

    # SKULL: Pure cortical/trabecular bone
    fslmaths $segmentation_skull_air_PVC_skin_air_PVC -thr 7 -uthr 7 -bin $skull_mask_PVC

    ################################################################################
    ##### STEP 5: UTE-HU MAPPING #####
    ################################################################################
    # Mapping from intensities to pHU

    pCT_soft_tissues_value="42"
    pCT_air_value="-1000"

    # Mapping: Air Fraction = -1000 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
    fslmaths $air_mask -mul $pCT_air_value -uthr -1 $pCT_air -odt float  # -odt prevents int overflow

    # Mapping: Soft-Tissue Fraction = 42 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
    fslmaths $soft_tissues_mask -mul $pCT_soft_tissues_value $pCT_soft_tissues -odt float

    # Mapping: Dilated Skull Fraction

    # Create a file to store the skull mapping
    touch ${path_pct}/pCT_skull_mapping.txt

    # Linear mapping: map image intensities in skull to Hounsfield units (HU)
    case "$skullmapping" in
        "miscouridou")
            echo "Mapping skull ~ miscouridou... pHU = -2085 UTE + 2329 "
            fslmaths $ute_stnorm -mul -2085 -add 2329 -mul $skull_mask_PVC $pCT_skull
            ;;
        "carpino")
            echo "Mapping skull ~ carpino... pHU = -2194 UTE + 2236"
            fslmaths $ute_stnorm -mul -2194 -add 2236 -mul $skull_mask_PVC $pCT_skull
            ;;
        "wiesinger")
            echo "Mapping skull ~ wiesinger... pHU = -2000 (UTE-1) + 42"
            fslmaths $ute_stnorm -add -1 -mul -2000 -add 42 -mul $skull_mask_PVC $pCT_skull
            ;;
        "treeby")
            echo "Mapping skull ~ treeby... pHU = -2929.6 UTE + 3247.9"
            fslmaths $ute_stnorm -mul -2929.6 -add 3247.9 -mul $skull_mask_PVC $pCT_skull
            ;;
        "kosciessa")
            echo "Mapping skull ~ kosciessa (default)... establishing linear skull mapping"
            # Call MATLAB to run the pct_skullmapping.m script
            echo "Running MATLAB script pct_skullmapping.m..."
            cd "${path_fun_pct}" || { echo "Directory not found"; exit 1; }
            matlab_command="pct_skullmapping(\"$sub_id\",\"$path_simnibs\")"
            $path_matlab -nodisplay -r "$matlab_command" || { echo "MATLAB function failed"; exit 1; }
            # Move to pCT folder (if not there already)
            cd "${path_pct}" || { echo "Directory not found"; exit 1; }
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
            fslmaths $ute_stnorm -mul $m1 -add $m2 -mul $skull_mask_PVC $pCT_skull
            ;;
    esac
    # Record mapping variant
    echo "$skullmapping" >> ${path_pct}/pCT_skull_mapping.txt

    # Add up the different parts to create the initial PVE corrected pCT
    fslmaths $pCT_skull -add $pCT_air -add $pCT_soft_tissues $pCT_PVC

    # Smooting: kernel size = 3x3x3 voxels; smoothing factor = 0.8
    SmoothImage 3 $pCT_PVC 0.8 $pCT_PV_smooth
    
    # Use smoothed mask at skull boundary

    # === BOUNDARY SMOOTHING VARS (add with others in STEP 4) ===
    skull_mask_PVC_ero="${path_pct}/skull_mask_PVC_ero.nii.gz"
    skull_mask_PVC_dil="${path_pct}/skull_mask_PVC_dil.nii.gz"
    skull_boundary_mask="${path_pct}/skull_boundary_mask.nii.gz"
    boundary_smooth_region="${path_pct}/boundary_smooth_region.nii.gz"
    pCT_boundary_smooth_temp="${path_pct}/pCT_boundary_smooth_temp.nii.gz"

    # === BOUNDARY-ONLY SMOOTHING ===
    echo "Creating skull boundary smoothing..."

    # TRUE 1-VOXEL BOUNDARY: Dilated(outer) - Eroded(inner) = symmetric shell
    fslmaths $skull_mask_PVC -ero $skull_mask_PVC_ero
    fslmaths $skull_mask_PVC -dilF $skull_mask_PVC_dil  
    fslmaths $skull_mask_PVC_dil -sub $skull_mask_PVC_ero $skull_boundary_mask

    # Smooth buffer region (+1vox for kernel)
    fslmaths $skull_boundary_mask -dilF $boundary_smooth_region
    fslmaths $pCT_PVC -mas $boundary_smooth_region $pCT_boundary_smooth_temp

    # Smooth only boundary
    SmoothImage 3 $pCT_boundary_smooth_temp 0.8 $pCT_PV_smooth 1 -mask $skull_boundary_mask

    # === FINAL COMBINE ===
    # Raw interiors + boundary smooth + skull override
    fslmaths $skull_mask_PVC -add 1 -thr 1 -uthr 1 -mul $pCT_PVC -add $pCT_PV_smooth -add $pCT_skull $pCT

    # Copy pCT to simnibs directory
    cp $pCT ${path_simnibs}

    # If DEBUG=0, remove nifti images in pseudoCT folder
    if [ "$debug" -eq 0 ]; then
        echo "Cleanup: Removing nifti images in pseudoCT folder"
        rm ${path_pct}/*.nii.gz
    fi

    echo "Process complete. You can now run the simulation with the pseudoCT. Ensure both 'pseudoCT.nii.gz' and 'tissues_mask.nii.gz' are in the m2m segmentation subject folder."
    }