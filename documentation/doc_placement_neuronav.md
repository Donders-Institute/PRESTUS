# Neuronavigation read-in

Extracting neuronavigation coordinates for post-hoc ultrasound simulations is critical to verify and predict the focus of the realized stimulation sites during neuromodulation experiments. Neuronavigation systems provide spatial tracking of the ultrasound transducer position relative to the individual subject’s anatomy in real-time. 

PRESTUS provides functions that aim to facilitate the read-in and preprocessing of positions recorded with the Localite system. 
To get started, see an example demo (without provided data) in `/examples/demo_localite.m`.

Multiple Localite filetypes are supported:

### `TriggerMarkers`

Multiple triggers can be sent during stimulation to acquire updates to the transducer location over time. Recorded locations will be recorded by localite to files named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>`. Recorded locations will be averaged within a position (associated with a continuous time stamp stream). Multiple positions per recording session are supported. 

### `GUMMarkers`

These files contain a static location that can be directly read in (either via `position_transducer_localite` or using the more comprehensive `neuronav_` function suite as in `demo_localite.m`).

Average transducer and target locations can be provided in multiple coordinate spaces (native image space being the most relevant for PRESTUS). 
The functions also facilitate standard-space definitons (experimental), e.g., to plot average transducer positions in standard space.

### Neuronavigation Coordinate Transforms

See the [coordinate system documentation](doc_coordinate_systems.md).

Recorded localite coordinates are relative to a neuronavigation planning image [1] (e.g., `T1_forneuronav`) that was loaded in the neuronavigation software. 

SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image however, so PRESTUS initially applies a transform from the planning to the segmentation image to the coordinates before they can be ported to MNI space.

1.	Trigger to Voxel (`neuronav_convert_trigger_to_voxels`): → subject voxel indices
2.	Voxel → mm (`neuronav_grid_to_mm_batch`): → subject mm
3.	Native mm → MNI mm (`neuronav_apply_deformation`): → MNI mm

[Optional]: Interpolate average positions across subjects

4.	Aggregate means in MNI mm across subjects (`neuronav_get_group_mean_mni`)
5.	Group mean MNI mm → subject mm (`neuronav_apply_inverse_deformation`)
6.	Subject mm → subject voxel (`neuronav_mm_to_voxel`)

#### Coordinate Spaces in neuronav functions

- `neuronav_convert_trigger_to_voxels`
	- Localite/trigger data (RAS mm) → subject (native) voxel space
- `neuronav_grid_to_mm_batch`
	- Input: Subject voxel indices (i,j,k)
	- Operation: Uses subject’s T1 NIfTI affine transform
	- Output: Subject native mm coordinates (in scanner/world space)
- `neuronav_apply_deformation` (forward warp)
	- Input: Native mm locations
	- Operation: Queries subject’s Conform2MNI_nonl.nii.gz (subject-native → MNI (mm))
	- Output: MNI mm coordinates (for reporting/plotting/averaging)
- `neuronav_get_group_mean_mni`
	- Input: CSVs containing columns like Mtrans_pos_MNI_x/y/z (all from previous function outputs)
	- Operation: Aggregates these MNI mm locations across subjects/sessions, computes mean (in mm)
	- Output: MNI mm group-mean locations (still in millimeters, not voxels)
- `neuronav_apply_inverse_deformation` (inverse warp)
	- Input: MNI mm group-mean locations (x, y, z in mm)
	- Operation: Uses subject’s MNI2Conform_nonl.nii.gz nonlinear field to map standard points back into that subject’s native space.
	- Output: Native mm location (in subject’s world/scanner space, for that subject)
- `neuronav_mm_to_voxel` (helper)
	- Input: Native mm (x, y, z) location
	- Operation: Uses subject T1 affine to map mm → voxel (i, j, k)
	- Output: Subject voxel indices for subsequent image work or storage
- `neuronav_export_session_csv`
	- Columns include:
	- Native (subject) voxel indices for target/transducer: (i, j, k from subject’s grid)
	- Native mm coordinates (if desired): (x, y, z in subject’s scanner space)
	- MNI mm coordinates for each (from deformation step): (x, y, z in standard space)
	- MNI voxel indices (may be inaccurate with limited FOVs)