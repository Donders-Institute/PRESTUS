# Performing acoustic simulations

## Axisymmetric simulations (3D simulations on 2D input setup)

[kspaceFirstOrder2D](http://www.k-wave.org/documentation/kspaceFirstOrder2D.php) models an infinite cylinder, not a focusing shell. It uses a 2D Cartesian grid with **X and Y axes**. These are assumed to map onto the second and first dimensions of input matrices and coordinates (i.e., [y, x]), respectively. Please note that this at first may seem unintuitive, but it is designed to facilitate the mapping to 3D image space voxel coordinates [i,j,k], which traditionally map onto [Y,X,Z]. In the future, this could also be resolved by assuming a XY mapping for k-Wave, and introducing a remapping of voxel space coordinates.

Axisymmetric simulations (similar to [this example benchmark](https://github.com/ucl-bug/k-wave/blob/main/k-Wave/examples/example_at_focused_bowl_AS.m)) can be performed using [kspaceFirstOrderAS](http://www.k-wave.org/documentation/kspaceFirstOrderAS.php). It uses a 2D axisymmetric coordinate system with **axial and radial axes**. These are assumed to map onto the first and second dimensions of input matrices and coordinates, respectively.

Multiple adjustments are made if axisymmetry is requested:
- The grid should be specified as [radial x 2, axial].
- It is assumed that the (2x) radial axis is shorter than the axial axis.
- Internally, [radial x 2, axial] is remapped onto [axial radial]; only values from [half+1:end] are retained.
- Medium, source, and sensor are set up according to [axial, radial] dimensions.
- Following acoustic simulation, values are mirrored, i.e., output has [radial x 2, axial] dimensions again, in line with the default spec for kspaceFirstOrder2D. This is used as input for heating simualtions as well.

