# Specifying home- & heterogeneous medium properties

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

| Tissue           | Density [kg/m³] | Sound Speed [m/s] | Absorption Coefficient [dB/(cm·MHzy)]      | Absorption Power Law (y)| Thermal Conductivity [W/(m·K)]| Specific Heat [J/(kg·K)] | Perfusion (WIP) | Absorption Fraction* (WIP) |
|------------------|--------|----------|----------|-------|--------|------|-------|---|
| Water            | 994    | 1500     | 0.00217  | 2     | 0.6    | 4178 | 0     | 1 |
| Brain            | 1046   | 1546     | 0.59     | 1.2   | 0.51   | 3630 | 559   | 1 |
| Skin             | 1090   | 1610     | 0.4      | 1     | 0.37   | 3391 | 106   | 1 |
| Skull-trabecular | 1700   | 2300     | 13.3     | 1     | 0.32   | 2274 | 30    | 1 |
| Skull-cortical   | 1850   | 2800     | 13.3     | 1     | 0.3    | 1313 | 10    | 1 |

* The fraction of attenuation (alpha) that is assumed to be absorbed in the bioheat equation (heating simulations).

When only a binary skull layer is requested, the default assumes the medium properties of cortical bone.