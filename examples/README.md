# Example scripts

This folder contains example scripts that interact with the PRESTUS pipeline

- `calibration_standalone` initiates the calibration of transducer phases and amplitudes for equipment at the Donders institute.
- `createPhantom` can be used to generate 2D phantoms
- `createPseudoCT` can be used to start pseuodCT generation from UTE images
- `demo_localite` illustrates localite coordinate extraction for read-in in PRESTUS
- `hpc_compile_CUDA` contains instructions to set up and compile CUDA binaries for k-Wave

#### Further examples


- An example of a possible way to use these scrips for acoustic and heating simulations is provided as a DataLad dataset at https://gin.g-node.org/PRESTUS/sim_ernie/.
- A simplified 2D benchmarking example can be found [here](https://github.com/jkosciessa/PRESTUS_2D_demo).
- A (currently outdated) demo that focuses on transducer calibration can be found [here](https://github.com/jkosciessa/PRESTUS_bin/blob/main/tutorial/PRESTUS_intro_tutorial.md).