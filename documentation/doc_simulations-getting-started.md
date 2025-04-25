# How to design your own processing pipeline

### Packages to install (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4
- k-wave 1.4 (place in PRESTUS' toolbox folder)

### Start with the tutorial (documentation/PRESTUS_intro_tutorial.md)
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- The structural data is segemented with SimNIBS 4's charm.
- The segmented data is translated to a matrix of acoustic properties.
- k-wave performs acoustic and thermal simulations.
- Select outputs are plotted.

### Run the single_subject_pipeline
Ideally, the full pipeline only needs a subject id (a number consisting of no more than 3 digits). For this, the configuration metadata must be fully specified in (a set of) config files placed in the 'configs' folder (see ```doc_config.md```). The ```load_parameters.m``` function imports the conifguration. The different parameters that can be used are found in the 'default_config.yaml'. The default config should generally not be altered. Instead, parameter changes should be specified via an application-specific config file or before performing the subject call in MATLAB).

### Create a config file
- Config files are created in the .yaml format.
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
- The default config shows which ones should and don't have to be changed.

### Importing a transducer location
- These can be specified in your config [T1_grid_space].
- Another option is to import transducer locations from Localite and translate them to the correct space using the `get_trans_pos_from_trigger_markers` and `ras_to_grid` functions.

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