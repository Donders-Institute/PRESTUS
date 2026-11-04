# Concepts

This page explains what PRESTUS does in one picture. Terms used throughout the documentation are defined in the [Glossary](doc_glossary.md).

## What PRESTUS does

PRESTUS estimates where an ultrasound transducer on the head focuses its energy, how strong the field is inside the brain, and how much it heats the tissue. It does so in a sequence of steps, each building on the previous one:

```
  T1 (+ T2/UTE) MRI of the participant
        │
  1. Segmentation        SimNIBS charm labels skin, skull, brain, …
        │
  2. Placement           where the transducer sits and where it aims
        │
  3. Grid & medium       a 3D grid of tissue properties (sound speed, density, attenuation)
        │
  4. Acoustic simulation the pressure field of the transducer through the head (k-Wave),
        │                scaled to the intensity you asked for
  5. Thermal simulation  heating over the sonication protocol (optional)
        │
  6. Outputs             metrics table, NIfTI maps, HTML report with safety checks
```

Every step reads its settings from a config file ([Configuring a study](doc_configuration.md)), and each can be switched on or off ([Modules](doc_modules.md)). A simulation in **free water** skips steps 1 to 3: no MRI is needed, which makes it the quickest way to try PRESTUS ([Your first simulation](doc_getting-started.md)).

Simulations are estimates. Tissue properties, especially of the skull, are uncertain, and results depend on the segmentation, the grid and the placement. PRESTUS is a research tool; it is not approved for clinical use.
