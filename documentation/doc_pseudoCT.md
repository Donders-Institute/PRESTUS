# (pseudo-)CT

![pct_ute_density](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/ute_pct_density.png)

PRESTUS supports the use of UTE-based images as a source of pseudo-Hounsfield units. These can replace the skull layer of a layered simulation for continuous tissue property mapping. 

> [!WARNING]
> (pseudo-)HU mapping is an experimental feature in active development.

### Creating pseudoCTs from UTE scans

0. **MR acquisition** 

    Acquire an MR sequence with detailed bone contrast 

    >**Example PETRA UTE parameters** (Carpino et al., 2023)
    TR: 3.32 ms
    TE: 0.07 ms
    voxel size: 0.8 mm3
    352 sagittal slices
    flip angle: 2°
    FoV: 294 mm

1. **SimNIBS segmentation**. 

    Perform a SimNIBS segmentation using a T1w and a PETRA UTE scan (instead of T2w) as inputs. You can either use the SimNIBS GUI or PRESTUS (**default**).

2. **pseudoCT generation**. 

Run `pct_create_pseudoCT.sh` in bash. An example script to create a pseudoCT is provided at `examples/createPseudoCT.sh`. To run the code you will need to install (or load the following modules):

- SimNIBS
- FSL
- ANTs

The pseudoCT and an associated mask file will be deposited in the `m2m` folder alongside the SimNIBS segmentation. 

#### pseudoCT generation

The following steps are used to create the pseudoCT:

- UTE: Threshold at 0, apply bias field correction to remove inhomogeneity, normalize (divide by soft tissue peak value > i.e., normalised UTE intensity = 1)
- Soft tissue value is identified with the MATLAB function `pct_soft_tissue_peak`
- pHU Soft-tissue (normalised UTE intensity = 1): 42 HU (Wiesinger et al., 2018)
- pHU Air: -1000 HU (Wiesinger et al., 2018; Miscouridou et al., 2022)
- pHU Skull: Linear mapping
    - `miscouridou` | pHU = −2085 UTE + 2329
    - `carpino`     | pHU = -2194 UTE + 2236
    - `wiesinger`   | pHU = -2000 (UTE-1) + 42
    - `treeby`      | pHU = -2929.6 UTE + 3247.9
    - `kosciessa`   | individualized: pHU_trabecular=300 & pHU_trabecular_cortical = 700
- 3D Gaussian smoothing (σ=0.8 voxels via ANTs’ SmoothImage) prevents abrupt tissue transitions in the simulation grid. To retain information in the skull layer, unsmoothed values are retained within an eroded skull mask that attempts to correct for partial volume effects.

<span style="font-size:11px;">
**References:**    
Wiesinger, F. et al. Zero TE-based pseudo-CT image conversion in the head and its application in PET/MR attenuation correction and MR-guided radiation therapy planning. Magn. Reson. Med. 80, 1440–1451 (2018).    
Miscouridou, M., Pineda-Pardo, J. A., Stagg, C. J., Treeby, B. E. & Stanziola, A. Classical and Learned MR to Pseudo-CT Mappings for Accurate Transcranial Ultrasound Simulation. IEEE Trans. Ultrason., Ferroelectr., Freq. Control 69, 2896–2905 (2022).  
Fat, D. L. et al. The Hounsfield value for cortical bone geometry in the proximal humerus—an in vitro study. Skelet. Radiol. 41, 557–568 (2012).    
Carpino et al. (2024). Transcranial ultrasonic stimulation of the human amygdala to modulate threat learning. MSc thesis.   
</span>

### Using (pseudo-)CTs to inform acoustic properties

Hounsfield Units can be used to inform acoustic properties in the bone layer of a SimNIBS segmentation. This is a subtype of a `layered` medium.

To inform skull properties by pCTs in simulations, set `parameters.use_pseudoCT = 1`, define `parameters.t2_path_template` as the `pseudoCT.nii.gz` in the simnibs output directory, and choose the desired acoustic parameter mapping algorithms. The current code supports the following mappings optionally model density, speed of sounds, and attenuation in the skull bone. When set to `"none"` or when the parameter field is deleted/missing, each property will revert to the uniform `skull` value. The defaults highlighted below are the ones specified in the `default_config.yaml`. 

>[!warning] The most suitable mapping remains an active area of research. All mappings should therefore be treated as experimental. See [this issue](https://github.com/Donders-Institute/PRESTUS/issues/43). Mappings will exclusively be applied within the uniform skull mask.

#### Mapping skull density

Mapping algorithm is goverened by `pct_mapping_density`.

- `k-plan` | 4-part piecewise linear fit based on k-Plan defaults (**default**) <br>

    **Acoustic simulations**    
    Piece-wise linear mapping between HU and mass density in kg/m^3 using k-Plan default calibration [see k-Plan documentation](https://dispatch.k-plan.io/static/docs/planning-images.html#ct-calibration):
    
    $$ \rho_\text{skull} = \text{fit_pairwiselinear}(\text{HU}, \text{density}) $$
    ```
    HU = [-990, 60, 1000, 1950]
    density = [1.2, 1060, 1530, 2150]
    ```

    **Thermal simulations**     
    Overwrite bone density prior to the thermal simulation:
   
    $$ \rho_\text{skull} = 1850 \, \text{kg/m}^3 $$

    <br>

- `k-wave` |  4-part piecewise linear fit based on Schneider et al. (1996) <br>

    PseudoCT values are initially shifted by +1000 and thresholded at 300 to align with [hounsfield2density](http://www.k-wave.org/documentation/hounsfield2density.php).   

    $$ \rho_\text{skull} = \text{hounsfield2density}(\text{HU} + 1000) $$

    Density estimates are regularized at the minimum to the water density specified in the configuration, and a max. skull density `rho_bone = 2100` kg/m3.

    <span style="font-size: 12px;">
    **References:**    
    Schneider, U., Pedroni, E., and Lomax A., "The calibration of CT Hounsfield units for radiotherapy treatment planning," Phys. Med. Biol., 41, pp. 111-124 (1996).
    </span>
    <br>

- `marsac` | Marsac et al. (2017) <br>

    This algorithm initially regularizes skull pHU to a range of `HU_min = 300` and `HU_max = 2000`. Note that in contrast to Marsac et al., 2017, PRESTUS regularizes pHU values instead of excluding pHU < HU_min from the skull mask as the latter is prone to create skull holes.

    $$ \rho_\text{skull} = \rho_\text{water} + (\rho_\text{bone} - \rho_\text{water}) \cdot \frac{\text{HU} - \text{HU}_\text{min}}{\text{HU}_\text{max} - \text{HU}_\text{min}} $$

    with max. skull density `rho_bone = 2100` kg/m3.

    <span style="font-size: 12px;">
    **References:**    
    Marsac, L. et al. Ex Vivo Optimisation of a Heterogeneous Speed of Sound Model of the Human Skull for Non-Invasive Transcranial Focused Ultrasound at 1 MHz. International Journal of Hyperthermia. 33. 635–645 (2017). <br>
    </span>
    <br>

- `aubry` | "Porosity-based" mapping based on Aubry et al. (2003) <br>

    This algorithm first estimates the porosity as:

    $$ \phi_\text{skull} = 1 - \frac{\text{HU}}{\text{HU}_\text{max}} $$

    This uses the max. skull HU denominator from Guo et al., 2019.  
    Skull density is then estimated as the mixture of bone and water density composites:
    
    $$ \rho_\text{skull} = \rho_\text{water} \cdot \phi_\text{skull} + \rho_\text{bone} \cdot (1 - \phi_\text{skull}) $$

    Respective tissue sound speed properties are extracted from the configuration.

    <span style="font-size: 12px;">
    **References:**    
    Aubry, J.-F., et al. Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans. *J Acoust Soc Am* **113**, 84–93 (2003).  
    Guo, S., et al. Feasibility of ultrashort echo time images using full-wave acoustic and thermal modeling for transcranial MRI-guided focused ultrasound (tcMRgFUS) planning. *Phys. Med. Biol.* **64**, 095008 (2019).  
    </span>
    <br>
    
#### Mapping skull sound speed

Mapping algorithm is goverened by `pct_mapping_soundspeed`.

- `k-plan` | Density-based mapping (**default**) <br>

    Implements [the k-Plan mapping](https://dispatch.k-plan.io/static/docs/simulation-pipeline.html#converting-the-primary-planning-image-to-material-properties):

    $$ c_\text{skull} = 1.33 \cdot \rho_\text{skull} + 167 $$

- `marsac` | Density-based mapping according to Marsac et al. (2017) <br>
    
    $$ c_\text{skull} = c_\text{water} + (c_\text{bone} - c_\text{water}) \cdot \frac{\rho_\text{skull} - \rho_\text{water}}{\rho_\text{bone} - \rho_\text{water}} $$

    with max. speed of sound in skull `c_skull = 3360` m/s and max. skull density `rho_bone = 2100` kg/m3. Respective water properties are extracted from the configuration.

    <span style="font-size: 12px;">
    **References:**    
    Marsac, L. et al. Ex Vivo Optimisation of a Heterogeneous Speed of Sound Model of the Human Skull for Non-Invasive Transcranial Focused Ultrasound at 1 MHz. International Journal of Hyperthermia. 33. 635–645 (2017). <br>
    </span>
    <br>

- `aubry` | "Porosity"-based mapping based on Aubry et al. (2003) <br>

    This algorithm first estimates the porosity as:

    $$ \phi_\text{skull} = 1 - \frac{\text{HU}}{\text{HU}_\text{max}} $$

    This uses the max. skull HU denominator from Guo et al., 2019.  
    Skull sound speed is then estimated as the mixture of bone and water sound speed composites:
    
    $$ c_\text{skull} = c_\text{water} \cdot \phi_\text{skull} + c_\text{bone} \cdot (1 - \phi_\text{skull}) $$

    Respective tissue sound speed properties are extracted from the configuration.

    <span style="font-size: 12px;">
    **References:**    
    Aubry, J.-F., et al. Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans. *J Acoust Soc Am* **113**, 84–93 (2003).  
    Guo, S., et al. Feasibility of ultrashort echo time images using full-wave acoustic and thermal modeling for transcranial MRI-guided focused ultrasound (tcMRgFUS) planning. *Phys. Med. Biol.* **64**, 095008 (2019).
    </span>
    <br>

#### Mapping skull attenuation

Mapping algorithm is goverened by `pct_mapping_attenuation`.

- `k-plan` | Uniform fixed (**default**)

    $$ \alpha_\text{skull} = \alpha_\text{coeff} $$

    k-Plan fixes the attenuation coeff. to `13.3` and frequency power law to `1`. To allow more flexibility, this variant reads in the alpha power values specified for the respective bone segmentation (`trabecular` or `cortical`) from the current config. To replicate k-Plan's setup, specify `alpha_coeff = 13.3` and `alpha_power = 1` in the skull bone medium configuration.

    <br>

- `mueller` | Algorithm specified in Mueller et al. (2017). <br>

    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \left(1 - \frac{\text{HU} - \text{HU}_\text{min}}{\text{HU}_\text{max} - \text{HU}_\text{min}}\right)^{0.5} $$

    with `α_min` = 4 and `α_max` = 8.7 (see Yakuub et al., 2023).     
    
    These 𝛼 estimates are based on 500 kHz (i.e., 𝛼(f); see Aubry, J.-F., 2022 for prior benchmark simulations). However, we require `alpha_0` in `𝛼(f) = alpha_0 x f[MHz] ^ y`. We therefore estimate `alpha_0 = 𝛼(f)/0.5^y` with  `y` being the specified `alpha_power`for the skull tissue.
    
    <span style="font-size: 12px;">
    **References:**    
    Aubry, J.-F. et al. Benchmark problems for transcranial ultrasound simulation: Intercomparison of compressional wave modelsa).  J. Acoust. Soc. Am. 152, 1003–1019 (2022).  
    Mueller, J. K., Ai, L., Bansal, P. & Legon, W. Numerical Evaluation of the Skull for Human Neuromodulation with Transcranial Focused Ultrasound. Journal of Neural Engineering. 14. 066012 (2017).  
    Yaakub, S. N. et al. Pseudo-CTs from T1-Weighted MRI for Planning of Low-Intensity Transcranial Focused Ultrasound Neuromodulation: An Open-Source Tool. Brain Stimulation. 16. 75–78 (2023).
    </span>     
    <br>

- `aubry` | "Porosity-based" mapping based on Aubry et al. (2003) <br>

    This algorithm first estimates the porosity as:

    $$ \phi_\text{skull} = 1 - \frac{\text{HU}}{\text{HU}_\text{max}} $$

    This uses the max. skull HU denominator from Guo et al., 2019.  
    Skull attenuation is then estimated as:
    
    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \cdot \phi_\text{skull}^{0.5} $$

    with `α_min = 0.2` and `α_max = 8` (Aubry et al., 2003).

    <span style="font-size: 12px;">
    **References:**    
    Aubry, J.-F., et al. Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans. *J Acoust Soc Am* **113**, 84–93 (2003).  
    Guo, S., et al. Feasibility of ultrashort echo time images using full-wave acoustic and thermal modeling for transcranial MRI-guided focused ultrasound (tcMRgFUS) planning. *Phys. Med. Biol.* **64**, 095008 (2019).  
    </span>
    <br>