## Transducer definition

PRESTUS does not provide default transducer calibrations. Instead, it encourages users to set up transducers and optimize settings to match their empirical setup.

See the full [parameter documentation](doc_parameters.md#transducer-specification).

#### Static parameters

The transducer supports two configurations:

- **Annular array** (`array_shape.type = annular`)
- **Matrix array** (`array_shape.type = matrix`)

---

##### General transducer parameters

`transducer.source_freq_hz`  
Central frequency of the acoustic source (in Hz).

##### Annular transducer parameters

`transducer.n_elements`  
Number of transducer elements.

`transducer.Elements_ID_mm`  
Inner diameter of each transducer element (in mm).

`transducer.Elements_OD_mm`  
Outer diameter of each transducer element (in mm).

`transducer.curv_radius_mm`  
Radius of curvature of the transducer bowl (in mm).

`transducer.dist_to_plane_mm`   
Distance from the geometric focus to the transducer plane (in mm).

---

##### Matrix transducer parameters

`transducer.matrix.steering`  
Steering mode.  
Options: `1D` (axial), `3D` (volumetric).

`transducer.matrix.element_shape`  
Shape of individual elements.  
Options: `rect`, `disc`, `bowl`.  
Element area is defined using rectangular dimensions and projected onto the selected shape.

`transducer.matrix.elem_height_mm`  
Element height (in mm).

`transducer.matrix.elem_width_mm`  
Element width (in mm).

`transducer.matrix.outer_diameter_mm`  
Outer diameter of the transducer (in mm).

---

###### Curvature (matrix only)

`transducer.matrix.is_curved`  
Enable curved (bowl-shaped) geometry.

`transducer.matrix.curved.curv_radius_mm`  
Radius of curvature (in mm).

`transducer.matrix.curved.dist_to_plane_mm`  
Distance from geometric focus to transducer plane (in mm).

---

###### Clover configuration (optional)

`transducer.matrix.is_clover_setup`  
Enable Clover (multi-aperture) configuration.

`transducer.matrix.clover.n_leaves`  
Number of leaves (max: 3).

`transducer.matrix.clover.ROC_parent`  
Radius of curvature of the combined system (in mm).

---

###### Element layout definition

`transducer.matrix.matrix_shape.type`  
Defines how element positions are generated.  
Options: `define_here`, `extract_from_file`.

---

####### Define layout manually

`transducer.matrix.matrix_shape.define_here.grid_shape.type`  
Grid type: `rect`, `fibonacci`, `fermat`.

**Note:** These options may contain bugs due to limited testing.

---

####### Load layout from file

`transducer.matrix.matrix_shape.extract_from_file.file_path`  
Path to file containing element coordinates (x, y, z).

`transducer.matrix.matrix_shape.extract_from_file.start_row`  
Row index where data starts.

`transducer.matrix.matrix_shape.extract_from_file.start_col`  
Column index where data starts.

`transducer.matrix.matrix_shape.extract_from_file.n_elements`  
Number of elements to extract.

---

###### Optional subset selection

`transducer.matrix.matrix_shape.extract_from_file.select_random_subset`  
Enable random element selection.

`transducer.matrix.matrix_shape.extract_from_file.subset.random_seed`  
Controls reproducibility.

`transducer.matrix.matrix_shape.extract_from_file.subset.subset_n_elements`  
Number of selected elements.

---

###### Optional ROC projection

`transducer.matrix.matrix_shape.extract_from_file.project_on_new_ROC`  
Project elements onto a new curvature.

`transducer.matrix.matrix_shape.extract_from_file.ROC_projection.new_ROC_mm`  
New radius of curvature (in mm).

---

#### Dynamic parameters: amplitude and focus

These parameters are traditionally calibrated in free-water simulations (see below).
**Note:** Calibration is currently implemented only for **annular arrays**.

`transducer.source_amp`     
Amplitude of the acoustic source (in Pa). [**CALIBRATED (currently annular only)**]

`transducer.source_phase_deg`   
Phase of the acoustic source (in degrees). [**CALIBRATED**]

`transducer.source_phase_rad`   
Phase of the acoustic source (in radians). [**CALIBRATED**]

For **matrix arrays**, phase delays are currently **not provided as input**.  
Instead, steering is performed automatically by computing the phase for each element based on `transducer.trans_pos` and `transducer.focus_pos`.

See [transducer calibration](doc_calibration.md).

#### Placement parameters

See also [transducer placement](doc_placement.md). 

In `water` simulations in combination with annular arrays, the placement will be automatically determined based on the edge of the PML layer.

`transducer.trans_pos`  
Position of transducer bowl (XYZ, T1 grid voxel space).
Used for placement and **phase calculation in matrix arrays**.

`transducer.focus_pos`  
Position of stimulation target (XYZ, T1 grid voxel space).
Defines the focal point and is used for **automatic steering in matrix arrays**.

#### Target distance parameters

![PRESTUS transducer distance definitions](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/transducer_distances.png)

The figure shows PRESTUS definitions for annular transducer distances.

PRESTUS expects one either of the following to be specified (or it will attempt a geometric distance calculation based on `trans_pos` and `focus_pos`, if specified).

`expected_focal_distance_ep`    
Expected distance from the transducer exit plane to the stimulation focus (mm). This often corresponds to the focal depth setting of the driving system - transducer calibration.

`expected_focal_distance_bowl`  
Expected distance from the transducer bowl to the stimulation focus (mm). Will be internally calculated based on `expected_focal_distance_ep`, `dist_to_plane_mm`, and `curv_radius_mm`.

## Transducer modeling

The `setup_source` function creates realistic transducer sources for k-Wave simulations, supporting both simple uniform elements and advanced curved arrays that produce focused beams. The function produces binary masks (active source locations) and time-varying pressure signals. Two setups are supported:

1.	Realistic kWaveArray Elements (**recommended**)
The  kWaveArray  class (`use_kWaveArray`) in the k-Wave toolbox enables realistic simulation of multi-element phased arrays by defining transducers via physical parameters (radius, element position, apodization) rather than discrete grid points. Elements render smoothly regardless of grid resolution. Independent delays/weights per element enable dynamic focusing and element overlap is automatically handled through weighted element contributions. Supports annular rings, axisymmetric 2D simulations, matrix arrays, and GPU acceleration. This is the recommended setup type.
<br>
2.	Simple Custom Elements
Creates single- or multi-element bowl-shaped (3D) or arc-shaped (2D) pressure source regions for each transducer element using geometric functions (k-Wave function: `makeBowl`). Each element gets identical continuous wave (CW) signals. Best for uniform, non-overlapping arrays where you want simple control. The discretization of grid points can introduce  staircasing artifacts in the resulting `source.p`  masks.

