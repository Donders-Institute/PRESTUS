# PseudoCT

PRESTUS supports the use of UTE-based images as a source of pseudo-Hounsfield units. 
This is currently only implemented for layered setups.

### Creating a pseudoCT from UTE scans

1) Perform a SimNIBS segmentation using a T1w and a PETRA UTE scan (instead of T2) as inputs. 
    - You can either use the SimNIBS GUI or PRESTUS (**default**).
2) Run in bash: `create_pseudoCT.sh`.
    - functions/create_pseudoCT.sh calls the MATLAB function find_soft_tissue_peak

The pseudoCT and an associated mask file will be deposited in the `m2m` folder alongside the SimNIBS segmentation. 

To run the code you will need to install (or load the following modules):
- SimNIBS
- FSL
- ANTs

**Note: Example PETRA UTE parameters** (Carpino et al., 2023)

- TR: 3.32 ms
- TE: 0.07 ms
- voxel size: 0.8 mm3
- 352 sagittal slices
- flip angle: 2°
- FoV: 294 mm

#### Conceptual overview

The following steps are used to create the pseudoCT:

- Threshold UTE image at 0
- Bias field correction to remove inhomogeneity
- (Identify soft tissue peak intensity value)
- Normalization: Divide by soft tissue peak value, i.e., normalised UTE intensity = 1
- Linear mapping: image intensities -> Hounsfield units (HU)
    - Skull *(mask defined as SimNibs tissue layers 7+8)*: -2194 UTE + 2236
        - This is motivated as follows (Carpino et al., 2023): *"The interface between cortical and trabecular bone, identified at 0.7 normalised UTE intensity, is mapped to 700 HU (Lim Fat et al., 2012)."*
    - Soft-tissue (normalised UTE intensity = 1): 42 HU (Wiesinger et al., 2016)
    - Air: 1000 HU (Wiesinger et al., 2016; Miscouridou et al., 2022)
- Thresholding/smoothing across various tissue masks
    - Gaussian smoothing kernel (3x3x3 voxels) is applied at the skin-air and skull-air interfaces
    - *currently not well documented*

The following is an example script that you can use for a `create_pseudoCT.sh` call.

```
#!/bin/bash

rootpath="$(pwd)/.."
rootpath=$(builtin cd $rootpath; pwd)

# Load necessary modules
module load ants
module load matlab
module load fsl

# Change directory to the scripts directory
scriptpath=${rootpath}/tools/PRESTUS/functions
cd "${scriptpath}" || { echo "Directory not found"; exit 1; }

# Source the script containing the function
source ${scriptpath}/create_pseudoCT.sh

subject_id="001" # 'sub-' will automatically be added
m2m_path=${rootpath}/data/simnibs

# Call the create_pseudoCT function
create_pseudoCT "$subject_id" "$m2m_path" "${scriptpath}"
```

### Starting acoustic + thermal simulations using pseudoCT

see [this issue](https://github.com/Donders-Institute/PRESTUS/issues/43)

To inform skull properties by pCTs in simulations, set `parameters.usepseudoCT = 1`, define `parameters.t2_path_template` as the `pseudoCT.nii.gi` in the simnibs output directory, and choose `parameters.pseudoCT_variant`. The current code supports the following variants to use pCTs to inform skull tissue parameters. Note that this affects only the pCT-to-tissueproperty conversion, only the above described procedure to derive pCTs is currently supported.

- `carpino` | (**default**) Algorithm described in Carpino et al. (2024). <br>

    All subsequent steps are only applied inside the skull mask.

    ```
    ρ_skull = hounsfield2density(HU+1000)
    c_skull = 1.33 * ρ_skull + 167
    α_skull = α_bone_min + (α_bone_max − α_bone_min) * (1 − (HU − HU_min) / (HU_max − HU_min))^0.5 [Mueller et al., 2017]
    ```

    - The k-Wave function hounsfield2density converts pseudo-HUs to density. Original pseudoCT values are initially shifted by 1000 and thresholded at 300 to align with [hounsfield2density](http://www.k-wave.org/documentation/hounsfield2density.php). 
    - Resulting density values are regularized to a minimum of the specified water density, and a maximum density of 2100 kg/m3. 
    - The sound speed c is calculated from density 𝜌 using the linear relationship: c = 1.33𝜌 + 167. This is identical to the [k-Plan estimation](https://dispatch.k-plan.io/static/docs/simulation-pipeline.html#evaluating-plans). Note that due to the density regularization, sound speed is implicitly regularized.
    - The absorption coefficient 𝛼 is derived from HU values according to formula (3) in Yaakub et al., with `α_bone_min` = 4 and `α_bone_max` = 8.7. For both `carpino` and `yakuub` variants, these 𝛼 bounds are based on estimates made at 500 kHz (i.e., 𝛼(f); see Aubry, J.-F., 2022 for prior benchmark simulations). However, we require ```alpha_0``` in ```𝛼(f) = alpha_0 x f[MHz] ^ y```. We therefore estimate ```alpha_0 = 𝛼(f)/0.5^y``` with  ```y``` being the specified ```alpha_power_true```for the skull tissue. 

    *Reference:* Adapted from Carpino et al. (2024). Transcranial ultrasonic stimulation of the human amygdala to modulate threat learning. MSc thesis.

- `yakuub` | Algorithm specified in Yaakub et al. (2023). <br>

    ```
    ρ_skull = ρ_water + (ρ_bone − ρ_water) * (HU − HU_min) / (HU_max − HU_min) [Marsac et al., 2017]
    c_skull = c_water + (c_bone − c_water) * (ρ_skull − ρ_water) / (ρ_bone − ρ_water) [Marsac et al., 2017]
    α_skull = α_bone_min + (α_bone_max − α_bone_min) * (1 − (HU − HU_min) / (HU_max − HU_min))^0.5 [Mueller et al., 2017]
    ```

    *References:* 
    - Yaakub, S. N. et al. Pseudo-CTs from T1-Weighted MRI for Planning of Low-Intensity Transcranial Focused Ultrasound Neuromodulation: An Open-Source Tool. Brain Stimulation. 16. 75–78 (2023). 
    - Marsac, L. et al. Ex Vivo Optimisation of a Heterogeneous Speed of Sound Model of the Human Skull for Non-Invasive Transcranial Focused Ultrasound at 1 MHz. International Journal of Hyperthermia. 33. 635–645 (2017).
    - Mueller, J. K., Ai, L., Bansal, P. & Legon, W. Numerical Evaluation of the Skull for Human Neuromodulation with Transcranial Focused Ultrasound. Journal of Neural Engineering. 14. 066012 (2017).

- `k-plan` | Algorithm described in Carpino et al. (2024) with more fixed skull tissue properties akin to k-Plan <br>

    ```
    ρ_skull = hounsfield2density(HU+1000)
    c_skull = 1.33 * ρ_skull + 167
    α_skull = alpha_0
    ```

    k-Plan fixes the absorption coeff. to `13.3` and power law to `1`. To allow more flexibility, this variant reads in the alpha power values specified for the respective bone segmentation (`trabecular` or `cortical`) from the current config. To replicate k-Plan's setup, specify `alpha_0 = 13.3` and `alpha_power_true = 1` in the configuration of all bone segmentations.
    
    Heating simulations will initially run the acoustic simulation as specified above, and then overwrite bone density (= `1850`) and sound speed (= `1.33 x 1850 + 166.7`) prior to starting the heating simulation.

    ```
    ρ_skull = 1850
    c_skull = 1.33 * ρ_skull + 167
    α_skull = alpha_0
    ```


