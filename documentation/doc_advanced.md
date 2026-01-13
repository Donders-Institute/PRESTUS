## Advanced workflows

#### Iterating input parameters

The base configuration of PRESTUS allows to specify a single setup, encompassing settings such as transducer specification, entry-target points, and temporal sequence characteristics. It may also be of interest to explore the effects of a range of different parameter settings on an outcome, for instance for benchamrking, but also for choosing a suitable stimulation sequence.

A general strategy toward such parameter iteration is to work on the basus of a single study setup configuration (as described in doc_config), whose values are iteratively overwritten in a MATLAB script. 

For many parameters, this can simply be done by specifying different values for ```parameters.<<xxx>>```. One exception is the intensity and depth setting of the transducer output, as these depend on the internal transducer calibration. PRESTUS by default provides a transducer calibration ("profiling"). This should in general be used for a given depth and free water intensity. However, this may not provide sufficient flexibility for exploring different output amplitudes. For that reason, PRERSTUS also provides the function ```transducer_calibration```, which can be used within an interative MATLAB loop to find suitable transducer phase and amplitude settings that closely replicate the desired free water output profile.

#### Concatenating successive simulations

parameters.adopted_heatmap | path to nifti file
parameters.adopted_cumulative_heat | path to nifti file
options.sequential_configs


#### Modeling multiple transducers (experimental support)

Multiple transducers can be specified in layered simulations. See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

```
transducers:
  - name: right
    n_elements: 10
    ...

  - name: left
    n_elements: 10
    ...
```

Limitations: 
- simulations that do not inclue a skull or layered setup do not model more than a single transducer.
- setup with kWaveArray not supported
- exit-plane related metrics refer to the first transducer
- Thermal diffusion is simulated for the COMBINED field, but ALL focal-plane time-course heating plots reflect ONLY the focal plane of the first transducer and MAY MISS HOTSPOTS near other beams. DO NOT use these 1D/2D plots as an exhaustive safety check, but ALWAYS inspect 3D maxT and CEM43 volumes (NIfTIs).