# Example scripts

This folder contains example scripts that interact with the PRESTUS pipeline

- `tutorial_ernie` (with `tutorial_ernie.yaml`) runs a head simulation on SimNIBS' Ernie example: segmentation, MNI placement, acoustic and thermal simulation (draft tutorial: `tutorial_ernie.md`).
- `example_report_ernie` (with `example_report_ernie.yaml`, an overlay on `tutorial_ernie.yaml`) adds heating, the free-water run, and ISO CEM43 to the tutorial and prints the metrics against the ITRUSST non-significant-risk limits (not yet run end to end).
- `example_fov_ernie` (with `example_fov_ernie.yaml`, an overlay on `tutorial_ernie.yaml`) runs the tutorial on the full grid and on a partial acoustic field of view and compares the metrics (not yet run end to end).
- `example_mixed_uncertainty_ernie` (with `example_mixed_uncertainty_ernie.yaml`, an overlay on `tutorial_ernie.yaml`) runs a mixed protocol in uncertainty mode: one transducer alone, two firing together, then the first stance again, whose acoustic field is reused (not yet run end to end).
- `example_sequence_ernie` (with `example_sequence_ernie.yaml`) is a hypothetical sonication sequence on Ernie: three stances one after another with pauses, simulated as one joint thermal run; open the YAML in the Nexus to see it in the *Sonication sequence* card (not yet run end to end).
- `example_sequence_pct.yaml` is a hypothetical study combining several features: a pseudo-CT skull, two transducer types, and a sonication sequence whose middle step fires two transducers together (a firing group, summed coherently within one async run), with a limited acoustic field of view and thermal simulation; open it in the Nexus (not yet run end to end).
- `calibration_standalone` initiates the calibration of transducer phases and amplitudes for equipment at the Donders institute.
- `demo_transducer_refinement` refines the shipped transducer library to your own hydrophone axial profiles (one call, writes a user library entry; runs out of the box on demonstration profiles).
- `example_uncertainty` runs the uncertainty quantification pipeline (five-variant simulation + combined report); see [Uncertainty](../doc/doc_uncertainty.md).
- `createPhantom` can be used to generate 2D phantoms
- `createPseudoCT` can be used to start pseudo-CT generation from UTE images
- `demo_localite` illustrates localite coordinate extraction for read-in in PRESTUS
- `hpc_compile_CUDA` contains instructions to set up and compile CUDA binaries for k-Wave

## Further examples

- A simplified 2D benchmarking example on phantoms (no SimNIBS needed) is in [PRESTUS_2D_demo](https://github.com/jkosciessa/PRESTUS_2D_demo).

Legacy examples. These are outdated and do not run on the current version without changes:

- [sim_ernie](https://gin.g-node.org/PRESTUS/sim_ernie/) (DataLad dataset): acoustic and heating simulations on the Ernie head, written for a 2024 PRESTUS. Use `tutorial_ernie` above instead; the dataset remains a source of the Ernie T1 and T2 images.
- [PRESTUS intro tutorial](https://github.com/jkosciessa/PRESTUS_bin/blob/main/tutorial/PRESTUS_intro_tutorial.md): transducer calibration. For current calibration, see `demo_transducer_refinement` and `calibration_standalone`.