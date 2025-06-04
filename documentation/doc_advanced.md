## Advanced workflows

#### Iterating input parameters

The base configuration of PRESTUS allows to specify a single setup, encompassing settings such as transducer specification, entry-target points, and temporal sequence characteristics. It may also be of interest to explore the effects of a range of different parameter settings on an outcome, for instance for benchamrking, but also for choosing a suitable stimulation sequence.

A general strategy toward such parameter iteration is to work on the basus of a single study setup configuration (as described in doc_config), whose values are iteratively overwritten in a MATLAB script. 

For many parameters, this can simply be done by specifying different values for ```parameters.<<xxx>>```. One exception is the intensity and depth setting of the transducer output, as these depend on the internal transducer calibration. PRESTUS by default provides a transducer calibration ("profiling"). This should in general be used for a given depth and free water intensity. However, this may not provide sufficient flexibility for exploring different output amplitudes. For that reason, PRERSTUS also provides the function ```transducer_calibration```, which can be used within an interative MATLAB loop to find suitable transducer phase and amplitude settings that closely replicate the desired free water output profile.

#### Concatenating successive simulations

parameters.adopted_heatmap | path to nifti file
parameters.adopted_cumulative_heat | path to nifti file
options.sequential_configs