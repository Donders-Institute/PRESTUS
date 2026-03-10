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
    tissues_mask="${path_pct}/tissues_mask.nii.gz"

    # Temporary file names
    t2w_reg="${path_simnibs}/T2_reg.nii.gz"
    ute_reg="${path_simnibs}/UTE_reg.nii.gz"
    ute_reg_thr0="${path_pct}/UTE_reg_thr0.nii.gz"
    ute_reg_thr0_corr="${path_pct}/UTE_reg_thr0_corr.nii.gz"
    ute_reg_bias="${path_pct}/UTE_thr0_BiasField.nii.gz"
    ute_norm="${path_pct}/UTE_STnorm.nii.gz"

    # echo "========================================================================="
    # echo "SEGMENTATION"
    # echo "========================================================================="

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

    echo "========================================================================="
    echo "SKULL MASK"
    echo "========================================================================="

    # Create a binary skull mask from final_tissues, thresholding to include skull (value 7 or 8)
    fslmaths th -thr 7 -uthr 8 -bin $skull_mask

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
    
    echo "========================================================================="
    echo "STEP 1: UTE THRESHOLD"
    echo "========================================================================="

    # Apply a threshold to the UTE image at 0, eliminating negative values
    fslmaths $ute_reg -thr 0 $ute_reg_thr0

    ################################################################################
    ##### STEP 2: UTE BIAS FIELD CORRECTION #####
    ################################################################################
    
    echo "========================================================================="
    echo "STEP 2: UTE BIAS FIELD CORRECTION"
    echo "========================================================================="

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

    echo "========================================================================="
    echo "STEP 3: UTE SOFT TISSUE NORMALISATION"
    echo "========================================================================="

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
    fslmaths $ute_reg_thr0_corr -div "$peak_value" $ute_norm

    # The skull fraction should now be between the noise/air peak (around 0) and the soft tissue peak (1).
    
    ###############################################################################
    ##### STEP 4: ADVANCED PVC + MAPPING + BOUNDARY SMOOTHING BLOCK
    # Purpose: Generate publication-quality pseudoCT for SimNIBS electromagnetic/acoustic simulations
    # 
    # WHY PARTIAL VOLUME CORRECTION (PVC) AT BOUNDARIES?
    # Problem: SimNIBS segmentation has ca. 1mm³ voxels. At tissue interfaces (skull-CSF, 
    #          skull-air, skin-air), MOST boundary voxels contain 30-70% mixed tissues.
    # Impact:  Without PVC, these "half skull, half CSF" voxels can get mislabeled as 
    #          PURE skull → wrong HU → 50-200% error in acoustic pressure/E-field.
    # Solution: Dilate/erode operations detect PVE voxels, label them separately (5=skin-air, 
    #           7=skull-air), then smooth ONLY these regions → smooth simulation gradients while 
    #           preserving bulk tissue accuracy.
    #
    # ORDERING (PVC → mapping → boundary-smooth → final):
    # 1. PVC: Correct skull-air + skin-air partial volume → PVC skull masks
    # 2. Mapping: UTE intensity → HU (air/soft tissue/skull)
    # 3. Smoothing of composite PCT (with a dilated skull mask)
    # 4. Replace PCV skull mask values with unsmoothed values
    ###############################################################################
        
    # Inputs
    segmentation="${path_simnibs}/final_tissues.nii.gz"                 # SimNIBS labels 0-15
    
    # PVC masks: skull-air + skin-air partial volume correction
    skull_mask="${path_pct}/skull_mask.nii.gz"                          # Raw skull (labels 7+8)
    skull_mask_dil="${path_pct}/skull_mask_dil.nii.gz"                  # Dilated skull (catches bone-air PVE voxels)
    skull_mask_inv="${path_pct}/skull_mask_inv.nii.gz"                  # Non-skull voxels
    skull_mask_7="${path_pct}/skull_mask_7.nii.gz"                      # Pure skull interior → label 7
    segmentation_mask="${path_pct}/segmentation_mask.nii.gz"            # Binary head (non-zero labels)
    segmentation_mask_no_dil_skull="${path_pct}/segmentation_mask_no_dil_skull.nii.gz"  # Head minus dilated skull
    skull_air_interface="${path_pct}/skull_air_interface.nii.gz"        # Skull-air PVE boundary (label=7)
    segmentation_no_skull="${path_pct}/segmentation_no_skull.nii.gz"    # Segmentation with no skull
    seg_skull_air_PVC="${path_pct}/seg_skull_air_PVC.nii.gz"            # Segmentation with skull-air PVC applied
    seg_skull_air_PVC_binv="${path_pct}/seg_skull_air_PVC_binv.nii.gz"  # Non-(skull+air) mask
    seg_skull_air_skin_air_PVC="${path_pct}/seg_skull_air_skin_air_PVC.nii.gz"  # Full PVC (skull-air + skin-air, label=5)
    segmentation_mask_dil="${path_pct}/segmentation_mask_dil.nii.gz"    # Dilated head contour
    segmentation_mask_ero="${path_pct}/segmentation_mask_ero.nii.gz"    # Eroded head contour
    outer_head_contour="${path_pct}/outer_head_contour.nii.gz"          # Skin-air PVE ring (1 voxel wide)
    soft_tissues_mask_1to6="${path_pct}/soft_tissues_mask_1to6.nii.gz"  # Brain/CSF/skin/muscle (SimNIBS 1-6)
    soft_tissues_mask_9to13="${path_pct}/soft_tissues_mask_9to13.nii.gz" # Eyes/optic/vessels (SimNIBS 9-13)
    soft_tissues_mask="${path_pct}/soft_tissues_mask.nii.gz"            # Combined soft tissues (1-6 + 9-13)
    air_mask="${path_pct}/air_mask.nii.gz"                              # Outside head (label=0)
    skull_mask_PVC="${path_pct}/skull_mask_PVC.nii.gz"                  # PURE cortical/trabecular bone (PVC-corrected label=7)
    skull_mask_PVC_dil="${path_pct}/skull_mask_PVC_dil.nii.gz"          # Dilated pure skull (outer surface)

    # pCT HU images
    pCT_air="${path_pct}/pCT_air.nii.gz"                                # Air fraction → -1000 HU
    pCT_soft_tissues="${path_pct}/pCT_soft_tissues.nii.gz"              # Soft tissues → 42 HU (literature)
    pCT_skull="${path_pct}/pCT_skull.nii.gz"                            # UTE-calibrated skull → 200-1500 HU
    pCT_PVC="${path_pct}/pCT_PVC.nii.gz"                                # Sum of pure tissue fractions (no smoothing)
    pCT_PV_smooth="${path_pct}/pCT_PV_smooth.nii.gz"                    # Smoothed pCT
    pCT="${path_pct}/pseudoCT.nii.gz"                                   # FINAL: bulk raw + boundaries smooth

    echo "========================================================================="
    echo "STEP 4.1: PARTIAL VOLUME CORRECTION - skull-air + skin-air interfaces"
    echo "========================================================================="
    echo "INTERFACES CORRECTED:"
    echo "  1. Skull-air:"
    echo "  2. Skull-CSF"
    echo "  3. Skin-air"
    echo ""

    cd "$path_pct" || { echo "ERROR: $path_pct not found"; exit 1; }

    # ---------------------------------------------------------------------------
    # STEP 4.1.1: Initial skull mask from SimNIBS segmentation
    # ---------------------------------------------------------------------------

    fslmaths $segmentation -thr 7 -uthr 8 -bin "$skull_mask"

    # ---------------------------------------------------------------------------
    # STEP 4.1.2: SKULL-AIR INTERFACE PVC
    # Problem: Voxels at skull outer/inner surface are mixed (bone+air, bone+CSF)
    # Solution: Dilate skull by 1 voxel → captures all PVE boundary voxels,
    #           then label them separately (value=7) for HU smoothing later
    # ---------------------------------------------------------------------------

    fslmaths "$skull_mask" -dilF -bin "$skull_mask_dil"
    fslmaths $segmentation -bin "$segmentation_mask"
    fslmaths "$segmentation_mask" -sub "$skull_mask_dil" "$segmentation_mask_no_dil_skull"
    fslmaths "$segmentation_mask_no_dil_skull" -add 1 "$segmentation_mask_no_dil_skull"
    fslmaths "$segmentation_mask_no_dil_skull" -binv -mul 7 "$skull_air_interface"

    fslmaths "$skull_mask" -binv "$skull_mask_inv"
    fslmaths $segmentation -mul "$skull_mask_inv" "$segmentation_no_skull"
    fslmaths "$skull_mask" -mul 7 "$skull_mask_7"
    fslmaths "$skull_mask_7" -add "$skull_air_interface" -add "$segmentation_no_skull" "$seg_skull_air_PVC"

    # ---------------------------------------------------------------------------
    # STEP 4.1.3: SKIN-AIR INTERFACE PVC
    # Problem: Scalp surface voxels are mixed (skin+air)
    # Solution: Head contour (dilate-erode) = 1-voxel ring at scalp surface,
    #           label these as skin-air PVE (value=5) for smooth HU transition
    # Impact:  Prevents sharp skin(42 HU) → air(-1000 HU) jumps → realistic CT
    # ---------------------------------------------------------------------------

    fslmaths "$segmentation_mask" -bin -dilF "$segmentation_mask_dil"
    fslmaths "$segmentation_mask" -bin -ero "$segmentation_mask_ero"
    fslmaths "$segmentation_mask_dil" -sub "$segmentation_mask_ero" "$outer_head_contour"

    fslmaths "$seg_skull_air_PVC" -binv "$seg_skull_air_PVC_binv"
    fslmaths "$outer_head_contour" -mul "$seg_skull_air_PVC_binv" "$outer_head_contour"
    fslmaths "$outer_head_contour" -mul 5 -add "$seg_skull_air_PVC" "$seg_skull_air_skin_air_PVC"

    # ---------------------------------------------------------------------------
    # STEP 4.1.4: PURE TISSUE MASKS (for accurate UTE→HU mapping)
    # Now that PVE boundaries are separated, extract PURE tissue regions:
    #   - air_mask = outside head
    #   - soft_tissues_mask = brain/CSF/skin (non-PVE, non-skull)
    #   - skull_mask_PVC = only 100% bone voxels (no PVE contamination)
    # These pure masks ensure UTE intensity maps to correct physics (especially skull!)
    # ---------------------------------------------------------------------------

    fslmaths "$seg_skull_air_skin_air_PVC" -binv "$air_mask"
    fslmaths "$seg_skull_air_skin_air_PVC" -thr 1 -uthr 6 -bin "$soft_tissues_mask_1to6"
    fslmaths "$seg_skull_air_skin_air_PVC" -thr 9 -uthr 13 -bin "$soft_tissues_mask_9to13"
    fslmaths "$soft_tissues_mask_1to6" -add "$soft_tissues_mask_9to13" "$soft_tissues_mask"
    # Remove dilated skull from soft tissues to prevent dilated skull overlap in pCT_PVC
    fslmaths "$soft_tissues_mask" -sub "$skull_mask_dil" -max 0 "$soft_tissues_mask"

    # Create skull masks
    fslmaths "$seg_skull_air_skin_air_PVC" -thr 7 -uthr 7 -bin "$skull_mask_PVC"
    fslmaths "$skull_mask_PVC" -dilF "$skull_mask_PVC_dil"

    echo "========================================================================="
    echo "STEP 4.2: UTE→HU MAPPING on PVC-corrected pure tissue masks"
    echo "========================================================================="

    cd "$path_pct"
    pCT_soft_val=42     # Literature: water-equivalent soft tissues (Wiesinger 2016, Miscouridou 2022)
    pCT_air_val=-1000   # Standard air HU

    fslmaths "$air_mask" -mul $pCT_air_val -uthr -1 "$pCT_air" -odt float
    fslmaths "$soft_tissues_mask" -mul $pCT_soft_val "$pCT_soft_tissues" -odt float

    # Skull: UTE→HU linear mapping (use dilated skull mask to avoid rough transitions!)
    case "$skullmapping" in
    miscouridou)
        fslmaths "$ute_norm" -mul -2085 -add 2329 -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    carpino)
        fslmaths "$ute_norm" -mul -2194 -add 2236 -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    wiesinger)
        fslmaths "$ute_norm" -add -1 -mul -2000 -add 42 -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    treeby)
        fslmaths "$ute_norm" -mul -2929.6 -add 3247.9 -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    kosciessa|*)
        m1=$(awk 'NR==1{print $1}' "${path_pct}/pCTskullmapping.txt")
        m2=$(awk 'NR==2{print $1}' "${path_pct}/pCTskullmapping.txt")
        fslmaths "$ute_norm" -mul "$m1" -add "$m2" -mul "$skull_mask_PVC_dil" "$pCT_skull" ;;
    esac

    # Combine pure tissues (NO smoothing yet → sharp transitions at boundaries)
    fslmaths "$pCT_air" -add "$pCT_soft_tissues" -add "$pCT_skull" "$pCT_PVC"

    echo "========================================================================="
    echo "STEP 4.3: SMOOTHING using ANTs: sigma=0.8 voxel "
    echo "========================================================================="
    echo ""

    SmoothImage 3 "$pCT_PVC" 0.8 "$pCT_PV_smooth"

    echo "========================================================================="
    echo "STEP 4.4: FINAL COMPOSITION: smoothed everywhere except pure skull interior"
    echo "========================================================================="
    echo ""
    echo "FINAL ASSEMBLY LOGIC:"
    echo "  1. Smooth base"
    echo "  2. Override ONLY skull_mask_PVC voxels with raw pCT_skull"
    echo ""

    fslmaths "$skull_mask_PVC" -binv tmp_non_skull_mask.nii.gz                              # Non PVC-skull mask
    fslmaths "$pCT_PV_smooth" -mas tmp_non_skull_mask.nii.gz tmp_smooth_non_skull.nii.gz    # Use smoothed values in non-skull mask
    fslmaths "$skull_mask_PVC" -mul "$pCT_skull" tmp_skull_raw.nii.gz                       # Use unsmoothed values in PVC skull mask
    fslmaths tmp_smooth_non_skull.nii.gz -add tmp_skull_raw.nii.gz "$pCT" -odt float        # Final: smooth(non-skull) + raw(skull)
    rm tmp_non_skull_mask.nii.gz tmp_smooth_non_skull.nii.gz tmp_skull_raw.nii.gz           # Cleanup

    # Copy pCT to simnibs root
    cp "$pCT" "${path_simnibs}/pseudoCT.nii.gz"

    # Create and store a final tissue mask for simulations
    # Step 1: Inverse mask (original labels outside skull, 0 inside)
    fslmaths "$skull_mask_PVC_dil" -binv tmp_invmask.nii.gz
    fslmaths "$segmentation" -mas tmp_invmask.nii.gz "$tissues_mask"
    # Step 2: Add bone=4 only in skull regions
    fslmaths "$skull_mask_PVC_dil" -bin -mul 4 tmp_skull4.nii.gz
    fslmaths "$tissues_mask" -add tmp_skull4.nii.gz "$tissues_mask"
    rm tmp_invmask.nii.gz tmp_skull4.nii.gz

    # Copy tissues_mask to simnibs root
    cp "$tissues_mask" "${path_simnibs}/tissues_mask.nii.gz"

    # If DEBUG=0, remove nifti images in pseudoCT folder
    if [ "$debug" -eq 0 ]; then
        echo "Cleanup: Removing nifti images in pseudoCT folder"
        rm ${path_pct}/*.nii.gz
    fi

    echo "Process complete. You can now run the simulation with the pseudoCT. Ensure both 'pseudoCT.nii.gz' and 'tissues_mask.nii.gz' are in the m2m segmentation subject folder."
}