# Advanced Modelling

!!! note "At a glance"
    - **What:** Pipelines that run several related simulations from one config: uncertainty bounds, several target intensities, several transducers, or a sequence of sonications.
    - **When you need it:** To bracket a safety estimate, compare intensities, model two transducers, or carry heat from one sonication into the next.
    - **Before you start:** A working single-run study config ([Your first simulation](doc_getting-started.md)).

## Advanced simulation modes

Each mode is documented on its own page:

| Mode | Description |
|---|---|
| [Uncertainty](doc_uncertainty.md) | Run the default medium and four ITRUSST-bound variants (Temp-H/L, Trans-H/L) to bracket plausible heating and in-situ intensity |
| [Multi-Intensity](doc_multiintensity.md) | Scale one acoustic simulation to multiple target ISPPAs in parallel |
| [Multi-Transducer](doc_multitransducer.md) | Coherent or asynchronous multi-transducer simulations |
| [Sequential](doc_sequential.md) | Chain simulations so that each run inherits the thermal state of the previous one |

Modes can be combined: uncertainty + sequential, multi-intensity + uncertainty, async multi-transducer + multi-intensity sweep. Composition rules and affix conventions are described on each mode's page, with combined-affix examples in [Sequential: Output affix handling](doc_sequential.md#output-affix-handling).

## Iterating parameters

A study config specifies a single setup: transducer, entry-target points and temporal sequence. To explore the effect of a parameter, for benchmarking or for choosing a stimulation sequence, work from one study config ([Your first simulation](doc_getting-started.md)) and overwrite values in a MATLAB script, setting `parameters.<field>` per iteration.

The exception is the intensity and depth setting of the transducer output, which depend on the transducer calibration. Where the built-in profiling does not provide sufficient flexibility, use `calibration_transducer` within an iterative MATLAB loop to find phase and amplitude settings that replicate the desired free-water output profile ([Transducer Calibration](doc_calibration.md)).
