Extracting neuronavigation coordinates for post-hoc ultrasound simulations is critical to verify and predict the focus of the realized stimulation sites during neuromodulation experiments. Neuronavigation systems provide spatial tracking of the ultrasound transducer position relative to the individual subject’s anatomy in real-time. 

PRESTUS provides functions that aim to facilitate the read-in and preprocessing of positions recorded with the Localite system. Multiple triggers can be sent during stimulation to acquire updates to the transducer location over time. Recorded locations will be recorded by localite to files named `TMSTrigger/TriggerMarkers_Coil0_<timestamp>'.

Recorded locations will be averaged within a position (associated with a continuous time stamp stream). Multiple positions per recording session are supported. Average transducer and target locations can be provided in multiple coordinate spaces (native image space being the most relevant for PRESTUS). The functions also facilitate standard-space definitons (experimental), e.g., to plot average transducer positions in standard space.

### Relevant spaces and coordinate systems

1. **Subject Space [planning image] (voxels)**
	- Type: Integer (i, j, k) indices. E.g., 93 116 78
	- Reference: Specific image grid (subject’s T1 MRI, subject’s own segmentation, or MNI template image)
	- Usage: Access data arrays; relates directly to matrix dimensions.
	- Orientation: Based on image storage; mapping real anatomy depends on image header “transform”.
    - Columns in output CSV:
	    - `Mtrans_pos_x`, `Mtrans_pos_y`, `Mtrans_pos_z`
	    - `Mtarget_pos_x`, `Mtarget_pos_y`, `Mtarget_pos_z`
2. **Subject Space [segmentation image] (voxels)**
    - same as above, but for the simnibs segmentation image
3. **World (Scanner) Space / RAS (Right-Anterior-Superior) (mm)**
	- Type: Real (x, y, z) coordinates, in millimeters.
	- Reference: Physical space of scanner (native) or MNI brain.
	- Usage: Anatomical localization; input for nonlinear warps and overlays across images.
	- Orientation: RAS (Right-Anterior-Superior, SPM/FreeSurfer/FSL convention); origin location varies by header (e.g., AC).
4. **MNI Space RAS (Right-Anterior-Superior) (mm)**
	- Type: Real (x, y, z) in millimeters.
	- Reference: Standard stereotactic space, e.g., MNI152 (often referred to as “world” space in standard template).
	- Usage: Cross-subject reporting, group analysis, publication, plotting in standard space.
	- Note: Not usually 0-centered; (0,0,0) refers to AC in MNI, which can be off-center in template volume.
    - Columns in output CSV (nonlinear deformation field (e.g., `Conform2MNI_nonl.nii.gz`) using the SimNIBS/CHARM pipeline):
	    - `Mtrans_pos_MNI_x`, `Mtrans_pos_MNI_y`, `Mtrans_pos_MNI_z`
	    - `Mtarget_pos_MNI_x`, `Mtarget_pos_MNI_y`, `Mtarget_pos_MNI_z`
5. **MNI Space (voxels)**
	- Type: Integer voxel (i, j, k) indices within a template MNI image.
	- Reference: Specific image (often 181×217×181 grid), for operations like drawing spheres or extracting image data.
	- Conversion: Use template’s affine matrix to convert MNI mm ↔ MNI voxel indices.
    ```
    template_affine = niftiinfo('MNI152_T1_1mm.nii.gz').Transform.T;
    inv_affine = inv(template_affine);
    mni_vox = (inv_affine * [mni_mm 1]')'; % yields [i j k 1], take [1:3]
    ```

### Coordinate Spaces used in PRESTUS neuronav functions

- neuronav_convert_trigger_to_voxels
	- Input: Localite/trigger data (usually in RAS mm)
	- Converts: RAS mm → subject (native) voxel space
	- Output: grid/voxel indices (i,j,k for subject’s T1) for targets/transducers
- neuronav_grid_to_mm_batch
	- Input: Subject voxel indices (i,j,k)
	- Operation: Uses subject’s T1 NIfTI affine transform
	- Output: Subject native mm coordinates (in scanner/world space)
- neuronav_apply_deformation (forward warp)
	- Input: Native mm locations
	- Operation: Queries subject’s Conform2MNI_nonl.nii.gz (subject-native → MNI (mm))
	- Output: MNI mm coordinates (for reporting/plotting/averaging)
- neuronav_get_group_mean_mni
	- Input: CSVs containing columns like Mtrans_pos_MNI_x/y/z (all from previous function outputs)
	- Operation: Aggregates these MNI mm locations across subjects/sessions, computes mean (in mm)
	- Output: MNI mm group-mean locations (still in millimeters, not voxels)
- neuronav_apply_inverse_deformation (inverse warp)
	- Input: MNI mm group-mean locations (x, y, z in mm)
	- Operation: Uses subject’s MNI2Conform_nonl.nii.gz nonlinear field to map standard points back into that subject’s native space.
	- Output: Native mm location (in subject’s world/scanner space, for that subject)
- neuronav_mm_to_voxel (helper)
	- Input: Native mm (x, y, z) location
	- Operation: Uses subject T1 affine to map mm → voxel (i, j, k)
	- Output: Subject voxel indices for subsequent image work or storage
- neuronav_export_session_csv
	- Columns include:
	- Native (subject) voxel indices for target/transducer: (i, j, k from subject’s grid)
	- Native mm coordinates (if desired): (x, y, z in subject’s scanner space)
	- MNI mm coordinates for each (from deformation step): (x, y, z in standard space)
	- MNI voxel indices (may be inaccurate with limited FOVs)

#### Example Transform order (per subject)
	1.	Trigger to Voxel (neuronav_convert_trigger_to_voxels): → subject voxel indices
	2.	Voxel → mm (neuronav_grid_to_mm_batch): → subject mm
	3.	Native mm → MNI mm (neuronav_apply_deformation): → MNI mm
	# [Optional]: Interpolate average positions across subjects
	4.	Aggregate means in MNI mm across subjects (neuronav_get_group_mean_mni)
	5.	Group mean MNI mm → subject mm (neuronav_apply_inverse_deformation)
	6.	Subject mm → subject voxel (neuronav_mm_to_voxel)

#### Available transforms

PRESTUS reuses available nonlinear transform matrices from SimNIBS.

- CHARM’s `Conform2MNI_nonl.nii.gz` is a nonlinear deformation field that maps points from an individual SUBJECT’S SCANNER/CONFORM (native) space to MNI (template) space, expressed in millimeters.
- `MNI2Conform_nonl.nii.gz` is the inverse of `Conform2MNI_nonl.nii.gz`. It is a 4D nonlinear deformation field created by CHARM/SimNIBS that maps coordinates from MNI (template) space back into the subject’s native millimeter (mm) space.

- **Note:**  recorded localite coordinates are relative to a neuronavigation planning image [1] (e.g., `T1_forneuronav`) that was loaded in the neuronavigation software. SimNIBS' nonlinear transform matrices to MNI are relative to the `final_tissues.nii.gz` segmentation image [2] however, so we first need to apply a transform from the planning to the segmentation image to the coordinates before they can be ported to MNI space.