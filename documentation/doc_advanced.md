## Advanced Workflows

#### Iterating Parameters

The base configuration of PRESTUS allows to specify a single setup, encompassing settings such as transducer specification, entry-target points, and temporal sequence characteristics. It may also be of interest to explore the effects of a range of different parameter settings on an outcome, for instance for benchamrking, but also for choosing a suitable stimulation sequence.

A general strategy toward such parameter iteration is to work on the basus of a single study setup configuration (as described in the [Quick Start Guide](doc_getting-started.md)), whose values are iteratively overwritten in a MATLAB script. 

For many parameters, this can simply be done by specifying different values for `parameters.<<xxx>>`. One exception is the intensity and depth setting of the transducer output, as these depend on the internal transducer calibration. PRESTUS by default provides a transducer calibration ("profiling"). This should in general be used for a given depth and free water intensity. However, this may not provide sufficient flexibility for exploring different output amplitudes. For that reason, PRERSTUS also provides the function `transducer_calibration`, which can be used within an interative MATLAB loop to find suitable transducer phase and amplitude settings that closely replicate the desired free water output profile ([see transducer calibration](doc_calibration.md)).

#### Sequential Simulations

Applicable when you stimulate from different coordinates in sequence (for stimulation in parallel, see [multi-transducer modeling](#multi-transducer-modeling)). Instead of starting each heating simulation with the default starting temperatures, you can start your nth stimulation with the temperature and CEM43 maps from the previous simulation. To do this, you need to feed the pipeline your other configurations.

Example:
Since each config is a structure, you can easily place multiple configs in one structure without any converting.
Let's say you have config_1, config_2 and config_3 and you want to run them in sequence.

If you want the heatmaps to not carry over, you would run your pipeline like this:
`prestus_pipeline_start(subject_id, config_1)`   
`prestus_pipeline_start(subject_id, config_2)`   
`prestus_pipeline_start(subject_id, config_3)`   

But now, you will also feed it the configs for each subsequent simulation:
`sequential_configs.config_2 = config_2`  
`sequential_configs.config_3 = config_3`  
`prestus_pipeline_start(subject_id, config_1, 'sequential_configs', sequential_configs)`   

Please note that you have to use the names `config_x` in the sequential_configs, and that you have to use integers. So names like `config_-5`, `config_0`, `config_1234` and `config_007`.

```
parameters.adopted_heatmap | path to nifti file
parameters.adopted_cumulative_heat | path to nifti file
options.sequential_configs
```

#### Multi-Transducer Modeling

Status: *experimental support*

Multiple transducers can be specified in layered simulations. See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

```
transducer:
  - name: right
    elem_n: 10
    ...

  - name: left
    elem_n: 10
    ...
```

Limitations: 
- simulations that do not inclue a skull or layered setup (e.g., 'water') only model the first transducer
- setup with kWaveArray not supported
- different source frequencies not supported (see `source_sensor_setup.m`)
- exit-plane related metrics refer to the first transducer
- Thermal diffusion is simulated for the COMBINED field, but ALL focal-plane time-course heating plots reflect ONLY the focal plane of the first transducer and MAY MISS HOTSPOTS near other beams. DO NOT use these 1D/2D plots as an exhaustive safety check, but ALWAYS inspect 3D maxT and CEM43 volumes (NIfTIs).