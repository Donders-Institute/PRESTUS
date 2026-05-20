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
`options.sequential_configs = sequential_configs` 
`prestus_pipeline_start(config_1, options)`

Please note that you have to use the names `config_x` in the sequential_configs, and that you have to use integers. So names like `config_-5`, `config_0`, `config_1234` and `config_007`.

##### Cache reuse in sequential simulations

When a follow-up simulation is dispatched, PRESTUS automatically manages which results are reused and which are recomputed:

| Stage | Default behaviour | How to override |
|---|---|---|
| Grid & medium setup | **Reused** from base run | Set `io.preproc_affix` in the sequential config |
| Acoustic simulation | **Re-run** with the follow-up config | Set `io.acoustic_cache_affix` to an existing affix to reuse prior acoustics |
| Thermal simulation | **Re-run** (cache file gets a `_seq<N>` suffix automatically) | Set `io.thermal_cache_affix` explicitly to control the cache filename |

The rationale is that the head geometry does not change between sequential targets, so preproc can always be shared. Acoustics default to a fresh run because the follow-up may target a different location or use different transducer settings. The thermal cache is always given a unique name to avoid the "already done" check silently loading the base run's thermal results.

To reuse the base run's acoustic results in the follow-up (e.g. same target, different timing protocol):

```matlab
options.sequential_configs.config_2.io.acoustic_cache_affix = config_1.io.output_affix;
```

```
io.adopted_heatmap           | path to NIfTI file (set automatically by the dispatcher)
io.adopted_cem43             | path to NIfTI file (set automatically by the dispatcher)
io.preproc_affix             | affix for grid/medium cache lookup
io.acoustic_cache_affix      | affix for acoustic cache lookup
io.thermal_cache_affix       | affix for thermal cache file
options.sequential_configs
```

#### Multi-Transducer Modeling

Status: *experimental support*

Multiple transducers can be specified in layered simulations. See [this pull request](https://github.com/Donders-Institute/PRESTUS/pull/100) for examples.

```yaml
transducer:
  - name: right
    elem_n: 10
    ...

  - name: left
    elem_n: 10
    ...
```

##### Temporal relation of multi-transducer firing

In practice, multiple transducers may fire simultaneously, with some temporal offset, or fully independently with different duty cycles or pulse sequences. PRESTUS has native support only for the simultaneous coherent case. The other scenarios require workarounds described below.

**Simultaneous coherent firing (native support)**

When multiple transducers are listed in the configuration, PRESTUS combines them into a single k-Wave source matrix. Their pressure fields superpose implicitly as k-Wave solves the linear acoustic PDE — this is the physically correct treatment for transducers that fire at the same frequency with a known phase relationship. All transducers must share the same fundamental frequency; a warning is issued if they differ.

This is the default mode and requires no additional configuration.

**Asynchronous / temporally non-overlapping firing**

If the transducers fire in separate, non-overlapping time windows (asynchronous duty cycles), their pressure fields never coexist. Because only one transducer is active at any instant, two distinct combination rules apply:

- **Acoustic safety metrics (ISPPA, MI, peak pressure):** the relevant quantity is the per-voxel maximum across transducers — the peak pressure at any voxel is set by whichever transducer produces the higher instantaneous pressure there.
- **Thermal heat deposition:** time-averaged heating accumulates from all transducers; the combined heat source entering the bioheat equation is the intensity sum across all transducers weighted by their respective duty cycles.

PRESTUS provides a dedicated `transducer_coupling: async` mode that handles both combination rules automatically. See [doc_async_transducer.md](doc_async_transducer.md) for configuration details, pipeline stages, and output files.

##### Complex Pressure Field

PRESTUS extracts the complex pressure amplitude (magnitude and phase at the fundamental frequency) from the steady-state time series recorded by k-Wave and makes it available as an optional cache output. This is not needed for the two scenarios above — but it enables a third workflow described here.

To save the complex field alongside the main results, set:

```yaml
io:
  save_p_complex: true
```

This writes a separate file `sub-XXX_<medium>_p_complex<affix>.mat` containing a single variable `p_complex` — a complex-valued array of the same spatial dimensions as `p_max_all`, where magnitude encodes peak pressure amplitude and phase encodes the spatial phase of the wavefield at the fundamental frequency. This may allow more complex results integration.

##### Known Limitations

- Simulations without a skull or layered setup (e.g. `medium: water`) only model the first transducer.
- `kWaveArray` mode is incompatible with multiple transducers.
- Different carrier frequencies across transducers are not supported.
- Exit-plane metrics (focal distance, `max_Isppa_after_exit_plane`, `real_focal_distance`) are computed only with respect to the first transducer.
- Focal-plane time-course heating plots reflect only the focal plane of the first transducer and may miss hotspots near other beams. Always inspect the 3D `maxT` and `CEM43` NIfTI volumes for a complete safety assessment.
- Post-hoc water-only simulations are not implemented for multiple transducers.