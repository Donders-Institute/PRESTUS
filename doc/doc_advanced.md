# Advanced Modeling

## Iterating parameters

The base configuration of PRESTUS allows you to specify a single setup encompassing transducer specification, entry-target points, and temporal sequence characteristics. To explore the effects of different parameter settings on an outcome — for benchmarking or choosing a stimulation sequence — a general strategy is to work from a single study config (see the [Getting Started guide](doc_getting-started.md)) and iteratively overwrite values in a MATLAB script.

For most parameters this means setting different values for `parameters.<field>`. One exception is the intensity and depth setting of the transducer output, as these depend on the internal transducer calibration. Where the built-in profiling does not provide sufficient flexibility, use `transducer_calibration` within an iterative MATLAB loop to find phase and amplitude settings that replicate the desired free-water output profile (see [Transducer Calibration](doc_calibration.md)).

---

## Advanced simulation modes

PRESTUS supports several advanced simulation modes. Each is documented on its own page:

| Mode | Description |
|---|---|
| [Uncertainty](doc_uncertainty.md) | Run default / liberal / conservative medium-property variants to bracket plausible in-situ exposure |
| [Multi-Intensity](doc_multiintensity.md) | Scale one acoustic simulation to multiple target ISPPAs in parallel |
| [Multi-Transducer](doc_multitransducer.md) | Coherent or asynchronous multi-transducer simulations |
| [Sequential](doc_sequential.md) | Chain simulations so that each run inherits the thermal state of the previous one |

Modes can be combined: uncertainty + sequential, multi-intensity + uncertainty, async multi-transducer + multi-intensity sweep. Composition rules and affix conventions are described in each mode's page, with combined-affix examples in [Sequential — Output affix handling](doc_sequential.md#output-affix-handling).
