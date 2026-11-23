# Concepts

This page shows the steps PRESTUS takes from an MRI to a simulated field and its safety metrics, and where each step is configured. The sequence helps to place the other pages and to see which inputs a result depends on. Terms used throughout the documentation are defined in the [Glossary](doc_glossary.md).

On this page:

- [The Simulation Pipeline](#the-simulation-pipeline): the pipeline from MRI to outputs, and what a free-water run leaves out.

## The Simulation Pipeline

PRESTUS estimates where an ultrasound transducer on the head focuses its energy, how strong the field is inside the brain, and how much it heats the tissue. It does so in a sequence of steps, each building on the previous one:

![Pipeline from MRI to outputs: Segmentation, Transducer, Placement, Grid & medium, Acoustic, Thermal (optional), Outputs](assets/figures/fig_pipeline-overview_diagram.svg)

- **MRI:** the participant's T1, optionally with a T2 or UTE scan.
- **[Segmentation](doc_segmentation.md):** SimNIBS charm labels skin, skull, brain, and other tissues.
- **[Transducer](doc_transducer.md):** geometry, drive, and calibration from the library.
- **[Placement](doc_placement.md):** where the transducer sits and where it aims.
- **[Grid & medium](doc_medium.md):** a 3D grid of tissue properties (sound speed, density, attenuation).
- **[Acoustic](doc_simulations-acoustic.md):** the pressure field of the transducer through the head (k-Wave), scaled to the intensity you asked for.
- **[Thermal](doc_simulations-thermal.md):** heating over the sonication protocol (optional).
- **[Outputs](doc_outputs.md):** metrics table, NIfTI maps, and an HTML report with safety checks.

Every step reads its settings from a config file ([Configuring a study](doc_configuration.md)), and each can be switched on or off ([Pipeline Modules](doc_modules.md)). A simulation in **free water** skips segmentation, and placement and grid set-up run in reduced form (positions on a homogeneous water grid): no MRI is needed, which makes it the quickest way to try PRESTUS ([First simulation](doc_first-simulation.md)).

Simulations are estimates. Tissue properties, especially of the skull, are uncertain, and results depend on the segmentation, the grid, and the placement. PRESTUS is a research tool; it is not approved for clinical use.
