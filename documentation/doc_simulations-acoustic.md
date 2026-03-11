# Acoustic simulations

PRESTUS supports 3D acoustic wave simulations using k-Wave’s pseudo-spectral method for time-domain 3D acoustic simulations. These simulations model compressional waves (P-waves), which involve pressure changes propagating through isotropic media *without* shear deformations. To describe wave propagation, they calculate pressure and velocity fields. They are faster, less memory-intensive, and require less assumptions than elastic wave simulations (e.g., using `kWaveElastic`; not implemented) because they ignore shear waves and complex material interactions. However, this also means that scattering via shear waves is not explicitly modeled (albeit contributing to the bulk attenuation coefficient resulting from absorption and scattering).

PRESTUS by default supports 3D simulations using `kspaceFirstOrder3D`.

See the available [simulation backends](doc_backend.md).

> **Turning off acoustic simulations.** PRESTUS uses the `run_acoustic_sims` flag to indicate whether to run simulations. If set to `false`, PRESTUS will first check for existing files and check whether `overwrite_files` is set to off to load the existing matrices. This is intended to make this more intuitive for workflows in which acoustic simulations are "turned off" to run follow-up thermal simulations. This behaviour is unique for this module (turning off other modules may not run the step at all incl. loading existing data).

## Axisymmetric simulations

For debugging, it can be useful to specify and model symmetric 2D tissues (e.g, 2D `phantom` or `water` media). However, [`kspaceFirstOrder2D`](http://www.k-wave.org/documentation/kspaceFirstOrder2D.php) models an infinite cylinder, not a focusing shell. It uses a 2D Cartesian grid with **X and Y axes**. These are assumed to map onto the second and first dimensions of input matrices and coordinates (i.e., [y, x]), respectively. This may at first may seem unintuitive, but it is designed to facilitate the mapping to 3D image space voxel coordinates [i,j,k], which traditionally map onto [Y,X,Z]. In the future, this could also be resolved by assuming a XY mapping for k-Wave, and introducing a remapping of voxel space coordinates.

To model a 3D axisymmetric problem using a 2D computational grid (as recommended), [`kspaceFirstOrderAS`](http://www.k-wave.org/documentation/kspaceFirstOrderAS.php) can be used (see [this example benchmark](https://github.com/ucl-bug/k-wave/blob/main/k-Wave/examples/example_at_focused_bowl_AS.m). It significantly reduces computational cost compared to a full 3D simulation. The input structures (`kgrid`, `medium`, `source`, and `sensor`) are defined similarly to 2D simulations but interpreted in the axisymmetric coordinate system. The function `kspaceFirstOrderAS` accounts for wave propagation in such axisymmetric media and is functionally similar to `kspaceFirstOrder2D` but adapted for axisymmetry. kspaceFirstOrderAS uses a 2D axisymmetric coordinate system with **axial and radial axes**. These are assumed to map onto the first and second dimensions of input matrices and coordinates, respectively.

Multiple adjustments are made if axisymmetry is requested:

- The grid should be specified as `[radial x 2, axial]`.
- It is assumed that the (2x) radial axis is shorter than the axial axis.
- Internally, `[radial x 2, axial]` is remapped onto `[axial radial]` and only values from `[half+1:end]` are retained.
- Medium, source, and sensor are set up according to [axial, radial] dimensions.
- Following acoustic simulation, values are either mirrored into a cartesian 2D grid, or radially expanded into a cartesion 3D grid (see the note below). In cartesian 2D, the output will have [radial x 2, axial] dimensions again, in line with the default specification for kspaceFirstOrder2D.

> **Axisymmetric 2D vs. 3D output**
> Axisymmetric simulations are effectively 3D simulations that are performed efficiently on 2D input grids. By default, PRESTUS assumes that users want to follow-up an axisymmetric acoustic simulation with a 3D heating simulation. When heating is requested in the simulation pipeline, axisymmetric acoustic simulation outputs will be provided in 3D. If 2D outputs are desired, heating has to be deactivated. 2D heating simulations are possible by first running an acoustic 2D simulation without heating, and then running a heating simulation with the prior output, while deactivating (or not overwriting) the acoustic simulation.


