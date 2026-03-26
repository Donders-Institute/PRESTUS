# Medium properties

### Medium specification

PRESTUS supports homo- and heterogeneous medium configurations. These can be specified with `parameter.simulation_medium`. The following options are supported out of the box.

- `water`

    Places a transducer into homogeneous water tissue (at the inner edge of the requested PML layer).
    <br>

- `layered`

    Heterogeneous tissue composition based on charm segmentation. By default, the following media are included: 
    
    - water
    - brain
    - skin
    - cortical skull bone 
    - trabecular skull bone
    
    Tissue can be removed (or added if more detailed segmentations and medium acoustic properties are included) via the ```layers``` configuration fields. The skull layer of a layered simulation can be informed by (pseudo-)CT images (see ```doc_pseudoCT```). The pCT implementation is a subtype of a `layered` simulation with a unified skull layer. If the segmentation for any requested layer is not available, requested layers will be removed from the specification. The output HTML provides an overview of which layers were modeled, along with their acoustic properties.
    <n>
    `Layered` media will be [preprocessed](doc_preproc.md) during grid setup.
    <br>

- `phantom`

    A 3D imaging-based simulation is not always necessary. For debugging and benchmarking, one may for instance externally design a 2D phantom, and run simulations in a well-defined space. As the 2D phantom is explicitly designed, no preprocessing (e.g., rotation, cropping, etc.) is necessary. Using the ```phantom``` flag allows the same tissue flexibility as the ```layered``` version, but maps segmentation files onto medium masks without any further image processing. The 2D segmentation phantom must be provided as a ```final_tissues.nii.gz``` file in a (dummy) SimNIBS output folder.
    <n>
    A `phantom` Nifti is expected to contain artificial layers of a SimNIBS segmentation, including the PML padding.
    <br>

### Medium acoustic properties

![medium_properties](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/medium_properties.png)

Accurate modeling of wave propagation through heterogeneous media requires valid assignment of acoustic parameters such as speed of sound, density, and attenuation to the separate tissue media. PRESTUS provides defaults based on the literature, while providing flexibility to alter assigned parameter values, e.g., to implement more conservative or liberal assumptions. 

The following default medium properties applt in a full layered simulation. Parameters are specified in the default_config.yaml file (see doc_parameters.md). While PRESTUS’s default parameters are chosen to align with commercially available tools and prior benchmarks, these should not be regarded as ground truths, and may need to be adjusted according to the goals at hand (e.g., a more conservative assessment of intensity and/or heating).

| Tissue           | Density [kg/m³] | Sound Speed [m/s] | Attenuation Coefficient [dB/(cm·MHzy)]      | Attenuation Power Law (y)| Thermal Conductivity [W/(m·K)]| Specific Heat [J/(kg·K)] | Perfusion [mL/min/kg] | Absorption Fraction |
|-----------------------|--------|----------|----------|-------|--------|------|-------|---|
| Water                 | 994    | 1500     | 0.00217  | 2     | 0.6    | 4178 | 0     | 1 |
| Brain                 | 1046   | 1546     | 0.59     | 1.2   | 0.51   | 3630 | 559   | 1 |
| Skin                  | 1090   | 1610     | 0.4      | 1     | 0.37   | 3391 | 106   | 1 |
| Skull (single layer)  | 1850   | 2800     | 13.3     | 1     | 0.32   | 1793 | 20    | 0.28 |
| Skull (cortical)      | 1850   | 2800     | 13.3     | 1     | 0.32   | 1313 | 10    | 0.28 |
| Skull (trabecular)    | 1700   | 2300     | 13.3     | 1     | 0.32   | 2274 | 30    | 0.28 |

> **Note: Absorption Fraction.** For thermal simulations, PRESTUS models the fraction of attenuation that is assumed to be absorbed in the bioheat equation. The current default estimate is informed as follows: The bulk attenuation coefficient includes contributions from both longitudinal and shear waves. However, longitudinal waves contribute most (~ 80%,  Wang et al., 2018) to absorption and heating in the skull (Pinton et al., 2012; White et al., 2006). For longitudinal and shear waves, absorption values of 2.7 dB/cm2/MHz, and 5.4 dB/cm2/MHz have been reported (Pinton et al., 2012), with a bulk attenuation of 13.3 dB/cm2/MHz. As such, we assume: absorption fraction = (longitudinal wave absorption + 0.2 * shear wave absorption) / bulk attenuation. This yields (2.7+0.2*5.4)/13.3 = 0.28. This estimate of the bone absorption fraction (28%) is substantially lower than the default of 100% implemented in k-Plan, but larger than the default implemented in BabelBrain (16%). More measurements are needed to decrease uncertainty about the relative acoustic energy absorption in the skull layer.

<br>
<span style="font-size: 12px;">
  <strong>References:</strong><br>
  Pinton, G., Aubry, J., Bossy, E., Muller, M., Pernot, M., & Tanter, M. Attenuation, scattering, and absorption of ultrasound in the skull bone. <em>Méd Phys</em> <strong>39</strong>, 299–307 (2012). https://doi.org/10.1118/1.3668316.<br>
  Wang, X.-D., Lin, W.-J., Su, C., & Wang, X.-M. Influence of mode conversions in the skull on transcranial focused ultrasound and temperature fields utilizing the wave field separation method: A numerical study. <em>Chin Phys B</em> <strong>27</strong>, 024302 (2018). https://doi.org/10.1088/1674-1056/27/2/024302.<br>
  White, P. J., Clement, G. T., & Hynynen, K. Longitudinal and shear mode ultrasound propagation in human skull bone. <em>Ultrasound Med Biol</em> <strong>32</strong>, 1085–1096 (2006). https://doi.org/10.1016/j.ultrasmedbio.2006.03.015.
</span>
<br>