# PRESTUS Parameters

Every configuration parameter is listed with its default and what it does, grouped by the YAML struct it belongs to. The structs are spread over six pages by topic:

| Page | Structs |
|---|---|
| [Run, Paths, and Outputs](doc_parameters_run.md) | `subject_id`, `platform`, `path`, `startup`, `io`, `modules`, `report`, `study`, `hpc` |
| [Transducer](doc_parameters_transducer.md) | `transducer` |
| [Calibration](doc_parameters_calibration.md) | `calibration` (`config_calibration.yaml`) |
| [Placement and Time Reversal](doc_parameters_placement.md) | `placement`, `time_reversal` |
| [Head Model and Grid](doc_parameters_head-model.md) | `segmentation`, `headmodel`, `pct`, `layers`, `medium_properties`, `grid` |
| [Simulation and Sonication](doc_parameters_simulation.md) | `simulation`, `acoustic_backend`, `thermal_backend`, `timing`, `thermal`, `analysis` |

To find one key, use the search box or **Find setting…** (Ctrl+F) in the Nexus.

Parameters are organised in nested structs that map directly to [YAML](https://yaml.org/) keys. `config_default.yaml` defines all parameters and their default values and should not be edited; a study config is merged on top and only needs the values that differ (see [Configuring a study](doc_configuration.md) for the full layering, project defaults, and overrides). Parameters can also be set in MATLAB before calling `prestus_pipeline_start(parameters)`. Mandatory parameters have no meaningful default and must always be provided.

