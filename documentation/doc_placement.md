# Transducer Placement

Simulations in free water will place the transducer in a homogeneous water medium without the need to specify a target. For head smulations, any transducer needs to be positioned close to the scalp (possibly at some distance to allow for coupling), and a target needs to be specified.

PRESTUS expects coordinates of both the transducer (`transducer.trans_pos`) and target (`transducer.focus_pos`) to be specified as voxels in the T1w grid. The transducer coordinate describes the bowl of the transducer (i.e., not the exit plane in the case of curved transducers). Both coordinates need to be reported in the space of the planning image.

## Manual coordinate selection

Suitable coordinates (x/y/z) may be identified in preferred imaging software based on e.g. an fMRI hotspot, anatomical marker. 

## Heuristic coordinate selection

PRESTUS can identify heuristic locations for transducer placement (see [Heuristic Tranducer Placement](doc_placement_heuristic.md)). This benefits iterative approaches without manual intervention. 

Outside of PRESTUS, alternative tools such as [PlanTUS](https://github.com/mlueckel/PlanTUS) can be used to manually identify candidate transducer locations based on a broader range of criteria.

## Neuronavigation coordinate selection

PRESTUS provides helper functions to read-in coordinates acquired with  neuronavigation systems. Currently, PRESTUS supports read-in of localite positions.

An example file (`examples/demo_localite.m`) is provided. 

The example script `demo_localite.m` highlights a workflow for extracting localite transducer positions. 

- Based on the study design, localite may encode the positions from multiple stimulations in the same output file. For a requested session, `neuronav_select_and_average_localite` selects the latest available Localite trigger XML file. 
- In a localite session, triggers may be repeatedly acquired during the full pulse train repetition duration. `neuronav_compute_series_statistics` computes statistics across stimulus series (using a user-defined voxel size and user-defined trigger trains). 
- From these extracted metrics, `neuronav_create_marker_averags` generates averaged marker positions for transducer and target locations. 
- These are converted from Localite tracker space to native voxel indices and RAS coordinates aligned with the planning T1 image (`neuronav_convert_trigger_to_voxels`).
- Optionally, coordinates can be transformed to MNI space via SIMNIBS registration matrices (`neuronav_convert_native_to_MNI`). 

Results are exported to a CSV file in the Localite data folder, ready for PRESTUS simulation input.

> **Alternative localite read-in.**
> PRESTUS currently allows an alternative Localite read-in using the following parameters. This is not yet documented.
| `transducer_from_localite`        | Load transducer position from Localite files?.                                                                       |
| `reference_transducer_distance_mm` | Distance from tracker to transducer exit plane (in mm).                                                              |

