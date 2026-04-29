# PRESTUS GUI

The PRESTUS graphical user interface (`prestus_gui`) provides a point-and-click front-end to the pipeline. Every field maps directly to a YAML parameter — see [doc_parameters.md](doc_parameters.md) for the full reference.

![PRESTUS GUI](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/prestus_gui.png)

## Starting the GUI

```matlab
prestus_gui
```

Use **⬆ Load** to populate fields from an existing YAML, **⬇ Save** to export the current state, and **↺ Defaults** to reset to `config_default.yaml`.

### Configuration loading

The GUI and the pipeline share the same two-layer loading logic:

1. `config_default.yaml` is always loaded first and provides values for every parameter.
2. A study-specific YAML is merged on top via `MergeStruct` — only keys present in the study file overwrite the defaults; everything else retains its default value.

**In the GUI**, clicking **⬆ Load** calls `apply_params_to_gui` on the loaded file directly (no prior reset to defaults). This means:

- Fields present in the loaded YAML are updated.
- Fields **not** present in the loaded YAML keep whatever value is currently shown.

To get a clean state before loading, click **↺ Defaults** first, then **⬆ Load**.

**When running**, the GUI writes the current field values to a temporary YAML and passes it to `load_parameters`, which merges it on top of `config_default.yaml`. Fields left blank in the GUI are omitted from the saved YAML and therefore fall back to the default. The exception is `path.t2_pattern` — leaving it blank disables T2 entirely (`charm` is run with T1 only).

---

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
