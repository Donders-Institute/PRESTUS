# Placing from a Recorded Pose

![Pipeline steps; this page belongs to 3 Placement](assets/figures/fig_stage-strip-placement_diagram.svg){ .stage-strip }

Localite placement reads one recorded pose per transducer, either from the ingested markers of a session or from a single XML file you name. The series or marker it reads sets the position and beam axis of the simulated transducer, so a wrong `target` or `position` simulates another pose than the one used for stimulation.

On this page:

- [Automatic Mode: Ingested Markers](#automatic-mode-ingested-markers): reading the ingested pair of a subject and session.
- [File Mode: A Single XML File](#file-mode-a-single-xml-file): reading one TriggerMarkers, GUMMarkers, or InstrumentMarker file.
- [Several Transducers (Bilateral Setups)](#several-transducers-bilateral-setups): one pose per transducer.
- [Many Subjects and Sessions](#many-subjects-and-sessions): manifest rows from the Nexus, and the functions for standalone use.

!!! note "At a glance"
    - **When you need it:** You set `placement.mode: 'localite'` for a run.
    - **Before you start:** `transducer.focal_distance_ep` set (see [Placement from a Recorded Pose](doc_placement_neuronav.md#placement-from-a-recorded-pose)), and either sessions ingested as in [Localite Ingestion](doc_placement_neuronav_workflows.md) or a single Localite XML file. Folder layout and path keys: [Localite File Layout](doc_placement_neuronav_coordinates.md#localite-file-layout).

---

## Automatic Mode: Ingested Markers

Leave `placement.localite.file` empty. The ingested pair under `path.localite_post` is read. If it does not exist yet and `path.localite_raw` is set, ingestion runs first (with `neuronav.coil_map` if configured, otherwise for the requested series only, named `series_<position>` or `target`).

```yaml
path:
  localite_raw:  '/data/localite_raw'
  localite_post: '/data/derivatives/localite'
placement:
  mode: 'localite'
  localite:
    session: 2          # -> ses-02
    target: 'L_M1'      # ingested target name, or leave '' and use:
    position: 2         # detected-series index, matched against series_index in the JSON sidecar
```

To choose series and names interactively, ingest with `neuronav_gui` first; placement then only reads the result. A raw series that was not ingested is an error (`PRESTUS:localite:noSeries`), and an untracked series (position at the origin) is refused.

## File Mode: A Single XML File

Set `placement.localite.file`. TriggerMarkers, GUMMarkers, and InstrumentMarker files are accepted; the type is detected from the content.

```yaml
placement:
  mode: 'localite'
  localite:
    file: '/path/to/TriggerMarkers_Coil0_20240820122011842.xml'
    position: 2        # detected series in the file (1-based, after the pulse-count filter), or:
    target: ''         # description of an InstrumentMarker / GUMMarkers entry
```

## Several Transducers (Bilateral Setups)

`target` and `position` take one entry per transducer; transducer *i* uses entry *i*, and each gets its own pose. A single value applies to every transducer. A list of the wrong length is an error (`PRESTUS:localite:listLength`).

```yaml
placement:
  mode: 'localite'
  localite:
    session: 2
    target: [L_M1, R_M1]   # transducer(1) -> L_M1, transducer(2) -> R_M1
```

## Many Subjects and Sessions

For a study, the Nexus's Study tab turns the `neuronav_gui` outputs into manifest rows (**Import from Localite…**). It reads `neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json` files under `path.localite_post`, and makes either one run per target, or one run per session with the targets assigned to transducers 1, 2, … (left first, alphabetical, or a given order). The runs use `placement.mode: localite`, so PRESTUS reads the ingested markers itself. See [Group analysis](doc_group_study-setup.md#importing-runs-from-localite).

??? note "Entry points and standalone use"
    Placement runs in `preproc_transducer_placement`; the older `placement.localite.enabled = 1` switch in `preproc_head` uses the same code. Both go through `neuronav_localite_pose`, and the pose's column 4 (position) and column 1 (beam direction) go to `localite_matrix_to_positions`. All paths share the same XML parser and position math.

    For multi-position or multi-session workflows outside the pipeline, call `neuronav_ingest_markers` → `neuronav_select_localite` → `neuronav_compute_series_statistics` → `neuronav_convert_trigger_to_voxels` directly (see `examples/demo_localite.m` and [Target coordinates from several sources](doc_placement_neuronav_targets.md#target-coordinates-from-several-sources)). Average transducer and target locations can be provided in several coordinate spaces; native image space is the one PRESTUS uses. Standard-space definitions (experimental) allow e.g. plotting average transducer positions in standard space.
