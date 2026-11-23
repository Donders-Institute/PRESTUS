# PRESTUS Nexus (GUI)

The Nexus is a desktop app to set up, check, and run a PRESTUS study, and to look at its placements and results. You fill in forms, the Nexus checks them as you type and writes the study config, and [MATLAB](https://www.mathworks.com/products/matlab.html) runs the simulation. The config is an ordinary PRESTUS [YAML](https://yaml.org/) that holds only the values that differ from the defaults, so it also runs without the Nexus.

On this page:

- [Installation](#installation): the `prestus-nexus` conda environment and other ways to install.
- [Starting the Nexus](#starting-the-nexus): opening the Nexus from MATLAB or a terminal, with a study config and project defaults.
- [A First Simulation](#a-first-simulation): from a T1 image to results in six steps.
- [Beyond the First Simulation](#beyond-the-first-simulation): where to go for studies, protocols, calibration, Localite, and remote runs.

!!! note "At a glance"
    - **When you need it:** To build a study config without editing YAML, check it before it runs, and run it locally or on a cluster.
    - **Before you start:** The `prestus-nexus` [conda](https://docs.conda.io/) environment ([Installation](#installation)) and MATLAB, which runs the simulations.

Subpages:

- [Tabs and Pages](doc_nexus_tabs.md): every tab and page with the config keys it writes, for anyone changing settings beyond the QuickStart tab, running a study of multiple subjects, or following runs and results.

## Installation

| System | Run in the PRESTUS folder |
|---|---|
| macOS, Linux | `bash nexus/install.sh` |
| Windows | `nexus\install.bat` (or double-click it) |

The installer creates the conda environment `prestus-nexus` from `nexus/environment.yml`. If conda is missing, it offers to install [Miniforge](https://github.com/conda-forge/miniforge) into your home folder first. Run it again after a PRESTUS update: it then updates the environment. With conda already set up, `conda env create -f environment.yml` in `PRESTUS/nexus` does the same. Installing without conda, the standalone app, and the cache folder are described under [Installation options](#installation-options).

### Updating PRESTUS

**Help ▸ Check for updates…** compares your PRESTUS folder with the latest release and, on **Update to …**, moves it there. Without a window (e.g. on a cluster), use the command line:

```bash
prestus-nexus --update --check     # report the installed and the latest release
prestus-nexus --update             # move to the latest release (--to v0.7.1: to that release)
python3 nexus/prestus_nexus/update.py --check   # the same without the Nexus environment
```

The update needs a PRESTUS folder made with `git clone`. It changes nothing when PRESTUS files were edited in the folder, or when the folder holds a version with changes of its own (a development version); the message then says what to do. The bundled toolboxes that are installed (such as k-Wave) follow the release. Afterwards, restart the Nexus and MATLAB; when the result says so, run the installer again first. Let running and queued simulations finish before updating: they would switch to the new version part-way through.

??? note "What an update leaves alone"
    - Untracked files (your study configs, results) and submodules you did not install.
    - A study's own copy of `config_default.yaml`: see [Project defaults](#project-defaults).
    - A PRESTUS folder downloaded as an archive: the check names the latest release, which you download yourself.

    The check reads the list of release tags from the PRESTUS repository on GitHub. Nothing is contacted unless you ask for a check or an update.

??? developer "How the update moves the checkout"
    A release is a git tag `vX.Y.Z` (`nexus/prestus_nexus/update.py`). Changed tracked files stop the update. On a branch the update fast-forwards to the tag (`git merge --ff-only`); a checkout that sits on a tag is switched to the new one (`git checkout --detach`), and only there `--to` can also go back to an earlier release. `--check` exits with 10 when a newer release exists, 0 when not, and 1 on an error.

### Installation Options

??? note "pip, Qt, and h5py"
    Without conda, `pip install -e .` in `PRESTUS/nexus` installs the same packages into any Python ≥ 3.10 environment; the list lives in `nexus/pyproject.toml`.

    Install [Qt](https://www.qt.io/) from [pip](https://pip.pypa.io/) only (as `environment.yml` does). conda's `pyqt` together with pip's `PyQt6-WebEngine` fails to load (`Library not loaded: …QtPrintSupport`). In an older environment built with conda's `pyqt`, run `conda remove pyqt` and then `pip install -e .`.

    Reading `-v7.3` HEAT files on the *Heating* sub-tab also needs [h5py](https://www.h5py.org/) (not installed by default; PRESTUS saves `-v7` unless MATLAB's preference says otherwise).

??? note "Cache folder"
    The Nexus runs from a local copy of its own code in the user's cache folder (`~/Library/Caches/prestus_nexus` on macOS, `~/.cache/prestus_nexus` on Linux); on the Donders project drive this halves the start-up time. The parsed `config/equipment` files are cached in the same folder while they are unchanged. Set `PRESTUS_NEXUS_CACHE` to use another cache folder.

??? developer "Local copy of the code"
    On start the Nexus copies the files of `nexus/prestus_nexus` that changed since the last start into the cache folder and runs from there, so Python does not read its modules from a network drive. The checkout stays the source: an edit there is copied on the next start. Set `PRESTUS_NEXUS_LOCAL_COPY=0` to run from the checkout directly.

`nexus/packaging/prestus_nexus.spec` builds a standalone app with [PyInstaller](https://pyinstaller.org/) (`pip install -e . pyinstaller`, then `pyinstaller packaging/prestus_nexus.spec` in `nexus/`). The app does not contain PRESTUS: on first start it asks for the PRESTUS folder (or use `--prestus-root` / `PRESTUS_ROOT`) and remembers it.

## Starting the Nexus

From MATLAB:

```matlab
prestus_nexus                                  % untitled study (File > New / Open)
prestus_nexus('config_study.yaml')             % open a study config
prestus_nexus(Python='/path/to/env/bin/python')
```

From a terminal:

```bash
conda activate prestus-nexus
prestus-nexus [config_study.yaml]
```

From MATLAB, `prestus_nexus` finds the environment's Python by itself: `$PRESTUS_NEXUS_PYTHON` if set, the Python of the last start, the conda environment `prestus-nexus`, then `python3`. It uses the first that can start the Nexus and remembers it. If none can, the error lists what it tried; pass `Python=` to choose one.

### Opening a Study Config

The Nexus opens a given config, or an untitled one. **File ▸ New study…** creates a config and can copy `config_default.yaml` into its folder, which pins the defaults as `prestus_config_init` does. **Open…** and **Open recent** switch to another config. The subject, paths, T1/T2 patterns, and SimNIBS folder are set on the Setup tab. The Issues panel checks them as you type: the T1 must resolve, and the `m2m_sub-XXX` folder must exist or [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) must be set.

??? developer "Compatibility with MATLAB's YAML loader"
    PRESTUS reads configs with `yaml.loadFile` (SnakeYAML, YAML 1.1). The Nexus:

    - quotes all strings (unquoted `yes`/`no`/`04:00:00` would change type);
    - reads floats such as `2.17e-3` the way SnakeYAML does;
    - writes every element of a transducer list with the same top-level keys, which `yaml.loadFile` requires.

    A round-trip test writes configs with the Nexus, loads them with `load_parameters`, and compares every value. It covers both a serial-based and a fully inline transducer.

??? developer "Config schema and overlay"
    The schema is derived from `config_default.yaml` at startup. `nexus/prestus_nexus/overlay.yaml` adds only widget kinds, allowed values, units, visibility rules, keys the code reads that the default omits, deprecated keys, readable field labels where a key does not read as words (`display`; other labels are derived in `labels.py`), readable drop-down labels for choice values (`choice_labels`; hovering an item shows the value written to the YAML), and the schematic icons shown with drop-down choices. Quote a value that contains a comma (`{label: "a, b"}`): unquoted, YAML splits it, and the drift test reports the stray key. `nexus/tests/test_schema.py` fails when the overlay no longer matches the default config or the code. See `nexus/README.md` for running the tests.

### Project Defaults

When the study folder holds a copy of `config_default.yaml` (from *File ▸ New study…* or `prestus_config_init`), the Nexus merges it onto the toolbox default exactly as `load_parameters` does: the copy's values win, keys added to PRESTUS since the copy was made keep their toolbox defaults. The *Review YAML* page names both files. **Tools ▸ Update project defaults…** replaces the copy with the toolbox default (keeping a dated backup), optionally after writing the values the copy pinned differently into the study config so the simulation settings stay the same; `prestus_config_upgrade(folder)` does the same from MATLAB.

## A First Simulation

The steps below take one participant from a T1 image to results. The screenshots show the demo subject `sub-009`.

### 1. Fill In the QuickStart Tab

![Nexus QuickStart tab with all four steps done](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_start-filled_screenshot.png)

The Nexus opens on *QuickStart*, which asks for four things and keeps every other setting at its PRESTUS default:

1. **Head image:** the participant's T1-weighted [NIfTI](https://nifti.nimh.nih.gov/). The Nexus finds the T2, the SimNIBS segmentation, and the subject number from it.
2. **Transducer:** the device and driving system from the equipment library, the focal depth, and the free-water intensity. The calibration PRESTUS will use is shown below.
3. **Target:** where to focus, in MNI coordinates. Give the target alone and PRESTUS finds a transducer position, or give both.
4. **Check and run:** where the simulation runs, and any errors still to fix.

Each step shows *Done* once it is complete. *Continue to Run* saves the study and opens the Run tab. The other tabs hold every setting, with the same values the QuickStart tab wrote.

### 2. Check the Transducer

![Nexus Transducer tab with device panel, cross-section, and calibration coverage](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_transducer_screenshot.png)

The *Transducer* tab shows what the library knows about the device. The device panel gives its focal range and the calibration entry PRESTUS will use; green values mean your focal depth lies inside that range. The cross-section draws the transducer and the focal depth you set (✕). The *Calibration* sub-tab (*Transducer ▸ Calibration*) plots how well the calibration fits at each depth, with your study's depth as a vertical line.

### 3. Check the Placement

![Nexus placement viewer with focus, bowl, beam, and FOV on the T1](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_placement-viewer_screenshot.png)

*Placement ▸ Viewer* draws the transducer, the beam, and the simulated field of view on the participant's T1 and segmentation. The readout on the right reports the placement: how deep the target lies, how far the transducer sits from the scalp, and how much skull the beam crosses. Red values need attention. In the example, the target lies 2.2 mm short of the [focal depth setting](doc_glossary.md#focal-depth-setting), which *Set focal depth to target* corrects.

To change a manual placement, choose under *Click sets* whether a click moves the focus or the transducer, then click a slice. *Snap to scalp* moves the transducer onto the scalp along the beam. Placements from [Localite](https://www.localite.de/), the heuristic search, or [PlanTUS](https://github.com/mlueckel/PlanTUS) are shown read-only.

### 4. Clear the Issues

![Nexus Issues panel listing one warning and three notes](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_issues-popup_screenshot.png)

The error and warning counts in the header open the Issues panel. It lists what MATLAB would otherwise report only when the run starts: misspelt keys, missing files, intensities that cannot be reached, and the placement viewer's red flags. Clicking an issue opens the field it refers to. Fix the errors; read the warnings and notes.

### 5. Run the Simulation

![Nexus Run tab with pipeline stages, HPC card, and run log](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_run_screenshot.png)

The *Run* tab lists the pipeline stages, which you can switch off for partial runs, and where to run. *Run locally* runs MATLAB on this computer and streams its log. *Submit to cluster* submits from a cluster login node; *Submit over SSH* submits from a laptop with the project drive mounted. Submitted jobs appear under *Jobs*. *Review YAML* shows the exact config that will be written.

The first run on a participant also creates the segmentation, so it takes longest.

### 6. Look at the Results

![Nexus Results tab showing the intensity map on the T1](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_results-acoustic_screenshot.png)

The *Results* tab reads what the runs wrote to the subject's output folder. *Run* picks a run when the folder holds several. *Acoustic* and *Thermal* show the maps on the T1, centred on their peak; clicking a slice reads the value there. *Report* opens the HTML report, *Data table* the metrics, and *Compare runs* puts the key metrics of all runs side by side, with the [air-cavity check](doc_head-model.md#air-cavity-check) columns when a run has them.

![Nexus Results tab showing tissue temperature time courses](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_results-heating_screenshot.png)

*Heating* plots the temperature of each tissue over the sonication and the steady-state after it.

## Beyond the First Simulation

| To | Go to | Details |
|---|---|---|
| Run multiple subjects or conditions | *Study* tab | [Study tab](doc_nexus_tabs.md#study-tab), [Group analysis](doc_group_study-setup.md#study-setup) |
| Set the sonication protocol, or exchange it with the TUS calculator | *Sonication* tab | [Timing diagram and TUS calculator](doc_nexus_tabs.md#timing-diagram-and-tus-calculator) |
| Calibrate a transducer to your hydrophone profiles | *Transducer ▸ Calibration* | [Calibration page](doc_nexus_tabs.md#calibration-page) |
| Ingest Localite neuronavigation recordings | *Placement ▸ Localite* | [Localite page](doc_nexus_tabs.md#localite-page) |
| Submit from a laptop over SSH | *Run* tab | [Run tab](doc_nexus_tabs.md#run-tab) |
| Rerun with another intensity, timing, or placement, or chain a sequential run | *Results ▸ Follow-up…* | [Follow-up runs](doc_nexus_tabs.md#follow-up-runs) |
| Undo a change, find a setting, find a tab's docs, or export a figure | *Edit* and *Help* menus, a figure's download button | [Menus and shortcuts](doc_nexus_tabs.md#menus-and-shortcuts) |
| See settings hidden as advanced | **Advanced** on any tab, or the advanced icon on a card | [Advanced options](doc_nexus_tabs.md#advanced-options) |
| Pin the defaults for a project | *File ▸ New study…* | [Project defaults](#project-defaults) |

[Tabs and Pages](doc_nexus_tabs.md) describes every tab, page, and control, with the config keys each one writes, including the Calibration and Localite pages that compute in MATLAB.
