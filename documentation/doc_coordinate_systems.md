# Coordinate Reference Systems

The function `transform_coordinates` can be used to transform coordinates between coordinate reference systems.

1. Subject Space: planning image [voxel indices] (`grid`)
	- Reference: Specific image grid (subject’s T1 MRI, segmentation)
	- Orientation: Based on image storage; mapping real anatomy depends on image header “transform”.
2. Subject Space: localite planning image [mm] (`ras_plus` with adjusted header)
    - Subject space in mm without any affine transform to world space (see [the localite planning image](doc_placement_heuristic.md#save-localite-planning-image)).
	- Use the header from the adjusted localite planning image for the transform.
2. World (Scanner) Space / RAS (Right-Anterior-Superior) [mm] (`ras_plus`)
	- Reference: Physical space of scanner (native) or MNI brain.
	- Orientation: RAS (Right-Anterior-Superior, SPM/FreeSurfer/FSL convention); origin location varies by header (e.g., AC).
4. MNI Space RAS (Right-Anterior-Superior) [mm] (`mni`)
	- Reference: Standard stereotactic space, e.g., MNI152 (often referred to as “world” space in standard template).
	- Note: Not usually 0-centered; (0,0,0) refers to AC in MNI, which can be off-center in template volume.
5. MNI Space [voxel indices] (`undefined`; see below)
	- Reference: Specific image (often 181×217×181 grid), for operations like drawing spheres or extracting image data.
	- Conversion: Use template’s affine matrix to convert MNI mm ↔ MNI voxel indices. 	

    ```
    template_affine = niftiinfo('MNI152_T1_1mm.nii.gz').Transform.T;
    inv_affine = inv(template_affine);
    mni_vox = (inv_affine * [mni_mm 1]')'; % yields [i j k 1], take [1:3]
    ```

#### Available transforms

PRESTUS reuses available nonlinear transform matrices from SimNIBS.

- CHARM’s `Conform2MNI_nonl.nii.gz` is a nonlinear deformation field that maps points from an individual SUBJECT’S SCANNER/CONFORM (native) space to MNI (template) space, expressed in millimeters.
- `MNI2Conform_nonl.nii.gz` is the inverse of `Conform2MNI_nonl.nii.gz`. It is a 4D nonlinear deformation field created by CHARM/SimNIBS that maps coordinates from MNI (template) space back into the subject’s native millimeter (mm) space.