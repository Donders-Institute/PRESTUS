# PRESTUS GUI

!!! note "At a glance"
    - **What:** A point-and-click editor for PRESTUS configs that can also start a run and show its results.
    - **When you need it:** If you prefer forms to editing YAML files. Prefer the [PRESTUS Nexus](doc_nexus.md); `prestus_gui` is deprecated.
    - **Before you start:** A MATLAB desktop session. A config made here is an ordinary YAML file ([Configuring a study](doc_configuration.md)).

!!! warning "Deprecated"
    `prestus_gui` is deprecated and will be removed in a future release. It still works but is no longer developed. Use the [PRESTUS Nexus](doc_nexus.md) to build, check and run study configs. The [Localite ingest GUI](#localite-ingest-gui) (`neuronav_gui`) is not affected.

`prestus_gui` is a point-and-click front-end to the pipeline. Every field maps to a YAML parameter ([doc_parameters.md](doc_parameters.md)).

![PRESTUS GUI](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/prestus_gui.png)

## Starting the GUI

```matlab
prestus_gui
```

Use **⬆ Load** to fill the fields from an existing YAML, **⬇ Save** to write the current state, and **↺ Defaults** to reset to `config_default.yaml`.

If no telemetry decision has been recorded yet, a dialog asks whether to enable anonymous usage statistics (Allow / Decline / Ask me later); see [doc_telemetry.md](doc_telemetry.md).

### Configuration loading

The GUI loads configs the way the pipeline does: defaults first, your study YAML merged on top ([Configuring a study](doc_configuration.md#how-a-config-is-assembled)).

**⬆ Load** updates the fields present in the loaded YAML and leaves the others as shown. For a clean state, click **↺ Defaults** first, then **⬆ Load**.

**When running**, the GUI writes the current field values to a temporary YAML and passes it to `load_parameters`, which merges it on top of `config_default.yaml`. Fields left blank are omitted from the saved YAML and fall back to the default. The exception is `path.t2_pattern`: leaving it blank disables T2 entirely (`charm` is run with T1 only).

??? note "Details: how Load merges"
    Clicking **⬆ Load** calls `apply_params_to_gui` on the loaded file directly, without a prior reset to defaults. Fields present in the loaded YAML are updated; fields not present keep whatever value is currently shown.

## macOS: Full Disk Access required

On macOS, MATLAB launched from the Finder or Dock runs under stricter sandboxing than MATLAB launched from a terminal. If full disk access is required, open **System Settings → Privacy & Security → Full Disk Access** and add `MATLAB.app`, then restart MATLAB. Terminal-launched MATLAB is not affected.

## Run

Press **▶ Run Simulation** to execute the pipeline with the current settings. The log area captures all terminal output produced during the run.

## Results

After a simulation completes, the **Results** tab shows orthogonal cross-sections of the output maps.

| Sub-tab | Content |
|---|---|
| **Acoustic** | X/Y/Z slices of the acoustic pressure / ISPPA map |
| **Thermal** | X/Y/Z slices of the temperature / CEM43 map |
| **Report** | Embedded HTML summary report |

Use **Refresh** to reload results, **Open HTML Report** to open the report in a browser, and **Open Output Folder** to reveal the output directory.

## Localite ingest GUI

`neuronav_gui` turns a raw Localite session folder into ingested, named targets, with a review step before anything is written. Its single-session and batch workflows are described in [Localite ingestion](doc_placement_neuronav_workflows.md#gui-workflow); the PRESTUS Nexus offers the same steps on its *Placement ▸ Localite* page ([Localite page](doc_nexus.md#localite-page)).
