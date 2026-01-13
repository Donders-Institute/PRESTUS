# Heuristic transducer placement

The function ```transducer_positioning``` identifies heuristic locations for transducer placement. 

The function targets brain regions specified by MNI coordinates, converting them to subject-native space using SimNIBS tools. It aims to identify candidate transducer positions where the geometric focus aligns with a target within a user-defined focal distance range, while ensuring >50% aperture overlaps skull/skin and low variance in distances. For each position, it computes the intersection proportion between the transducer and skin, mean/variance distances to skin and skull, and geometric focus/exit plane positions to ensure efficient coupling while minimizing aberrations. 

## Coordinate table

The `tpars_sub-XXX_target.csv` table contains transducer position candidates with the following columns. This information can be used to select a suitable candidate (e.g., by selecting a position with low intersection of the transducer and head tissue (prop_intersect < 0.05) and minimal Euclidean distance to the target; often supported by visual inspection).

| Column          | Description                                    |
|-----------------|------------------------------------------------|
| idx             | index                                          |
| trans_x         | Transducer x-coordinate (voxels)               |
| trans_y         | Transducer y-coordinate (voxels)               |
| trans_z         | Transducer z-coordinate (voxels)               |
| targ_x          | Target x-coordinate (voxels)                   |
| targ_y          | Target y-coordinate (voxels)                   |
| targ_z          | Target z-coordinate (voxels)                   |
| dist_to_target  | Euclidean distance transducer-to-target (voxels)|
| prop_intersect  | Proportion of transducer volume intersecting head |
| meandistskin    | Mean distance to skin surface in aperture plane (voxels) |
| vardistskin     | Variance of distances to skin (voxels²)        |
| meandistskull   | Mean distance to skull surface in aperture plane (voxels)|
| vardistskull    | Variance of distances to skull (voxels²)       |
