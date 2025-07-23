# Specifying homo- or heterogeneous medium properties

### Choose a simulation medium

PRESTUS supports various medium configurations. These can be specified with ```parameter.simulation_medium```. The following options are supported out of the box.

- ```water```

    Places a transducer into homogeneous water tissue.

- ```water_and_skull``` or ```brain_and_skull``` (deprecated)

    Use a baseline medium of either water or brain, and insert skull based on imaging masks. Constrained variant of the ```layered``` setup that may be deprecated in future releases.

- ```layered```

    Heterogeneous tissue composition based on charm segmentation. By default, the following media are included: 
    
    - water
    - brain
    - skin
    - cortical skull bone 
    - trabecular skull bone
    
    Tissue can be removed (or added if more detailed segmentations and tissue properties are included) via the ```layer_labels``` and ```medium``` configuration fields. The skull layer of a layered simulation can be informed by (pseudo-)CT images (see ```doc_pseudoCT```).

- ```phantom```

    A 3D imaging-based simulation is not always necessary. For benchmarking, one may for instance externally design a 2D phantom, and run simulations in a well-defined space. As the 2D phantom is explicitly designed, no preprocessing (e.g., rotation, cropping, etc.) is necessary. Using the ```phantom``` flag allows the same tissue flexibility as the ```layered``` version, but maps segmentation files onto medium masks without any further image processing. The 2D segmentation phantom must be provided as a ```final_tissues.nii.gz``` file in a (dummy) SimNIBS output folder.

### Default medium properties

The following are the current default medium properties in a full layered simualtion. Parameters are specified in the default_config.yaml file (see doc_config.md).
When a homogeneous skull layer is requested, the default assumes the medium properties of cortical bone.

| Tissue           | Density [kg/m³] | Sound Speed [m/s] | Attenuation Coefficient [dB/(cm·MHzy)]      | Attenuation Power Law (y)| Thermal Conductivity [W/(m·K)]| Specific Heat [J/(kg·K)] | Perfusion [mL/min/kg] | Absorption Fraction [Note1] |
|-----------------------|--------|----------|----------|-------|--------|------|-------|---|
| Water                 | 994    | 1500     | 0.00217  | 2     | 0.6    | 4178 | 0     | 1 |
| Brain                 | 1046   | 1546     | 0.59     | 1.2   | 0.51   | 3630 | 559   | 1 |
| Skin                  | 1090   | 1610     | 0.4      | 1     | 0.37   | 3391 | 106   | 1 |
| Skull (single layer)  | 1850   | 2800     | 13.3     | 1     | 0.32   | 1793 | 20    | 0.28 |
| Skull (cortical)      | 1850   | 2800     | 13.3     | 1     | 0.32   | 1313 | 10    | 0.28 |
| Skull (trabecular)    | 1700   | 2300     | 13.3     | 1     | 0.32   | 2274 | 30    | 0.28 |

[Note1] The fraction of attenuation (alpha) that is assumed to be absorbed in the bioheat equation (heating simulations). This is a best guess that is informed as follows: The bulk attenuation coefficient includes contributions from both longitudinal and shear waves. However, longitudinal waves contribute most (~ 80%,  Wang et al., 2018) to absorption and heating in the skull (Pinton et al., 2012; White et al., 2006). For longitudinal and shear waves, absorption values of 2.7 dB/cm2/MHz, and 5.4 dB/cm2/MHz have been reported (Pinton et al., 2012), with a bulk attenuation of 13.3 dB/cm2/MHz. As such, we assume: absorption fraction = (longitudinal wave absorption + 0.2 * shear wave absorption) / bulk attenuation. This yields (2.7+0.2*5.4)/13.3 = 0.28.

*References*

- Pinton G, Aubry J, Bossy E, Muller M, Pernot M, Tanter M. Attenuation, scattering, and absorption of ultrasound in the skull bone. Méd Phys 2012;39:299–307. https://doi.org/10.1118/1.3668316.
- Wang X-D, Lin W-J, Su C, Wang X-M. Influence of mode conversions in the skull on transcranial focused ultrasound and temperature fields utilizing the wave field separation method: A numerical study. Chin Phys B 2018;27:024302. https://doi.org/10.1088/1674-1056/27/2/024302.
- White PJ, Clement GT, Hynynen K. Longitudinal and shear mode ultrasound propagation in human skull bone. Ultrasound Med Biol 2006;32:1085–96. https://doi.org/10.1016/j.ultrasmedbio.2006.03.015.
