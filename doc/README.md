[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15095860.svg)](https://doi.org/10.5281/zenodo.15095860)

![PRESTUS logo](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/logo_PRESTUS.png){ .logo }

# PRESTUS: PREprocessing & Simulations for Transcranial Ultrasound Stimulation

PRESTUS is an open-source MATLAB toolbox for imaging-informed simulations of transcranial ultrasound stimulation (TUS). From a T1-weighted MRI it builds a head model, places a transducer, simulates the acoustic field and the heating it causes, and writes a report with the metrics used for planning and safety reporting. Simulations run with the [k-Wave engine](doc_backend.md) on a laptop, a GPU workstation or a cluster.

PRESTUS is intended for basic research. It is not approved for clinical use.

## Where to start

New to PRESTUS? Follow these in order:

1. [Concepts](doc_concepts.md): what PRESTUS does, in one picture. Unfamiliar terms are in the [Glossary](doc_glossary.md).
2. [Installation](doc_installation.md): MATLAB, PRESTUS and k-Wave are enough to begin.
3. [Your first simulation](doc_getting-started.md): a free-water run in a few minutes, then what a head simulation needs.
4. [Reading your results](doc_results.md): what to check in the report, and how to read the safety summary.
5. [Configuring a study](doc_configuration.md): your own config files for your own data.

Looking for something specific?

| I want to… | Go to |
|---|---|
| simulate a transducer I own | [Transducer Library](doc_transducer_library.md) |
| decide where the transducer goes | [Placement](doc_placement.md) |
| simulate heating for my protocol | [Thermal](doc_simulations-thermal.md) |
| simulate a pose recorded with Localite | [Localite placement](doc_placement_neuronav.md) |
| use subject-specific skull properties | [pseudo-CT](doc_pseudoCT.md) |
| run many subjects | [Studies & group analysis](doc_group.md) |
| run on a cluster | [Running on a cluster](doc_hpc.md) |
| build a config with a form instead of YAML | [PRESTUS Nexus](doc_nexus.md) |
| look up a setting | [Parameters](doc_parameters.md) |
| fix an error | [Troubleshooting](doc_troubleshooting.md) |

## What PRESTUS does

- [Segmentation](doc_preproc.md) of the T1 (and T2 or UTE) image with SimNIBS `charm`, and preprocessing into a simulation grid.
- [Medium mapping](doc_medium.md) of water, skin, skull (cortical and trabecular) and brain, optionally with a [(pseudo-)CT](doc_pseudoCT.md) for continuous skull properties.
- A [transducer library](doc_transducer_library.md) with calibrated phases and drive for the devices in use at the Donders Institute, and [calibration](doc_calibration.md) of your own against hydrophone profiles.
- [Placement](doc_placement.md) by typed coordinates, MNI coordinates, automatic search, or a pose recorded with neuronavigation.
- [Acoustic](doc_simulations-acoustic.md) and [thermal](doc_simulations-thermal.md) simulation of a sonication protocol, with k-Wave, [BabelBrain or jwave](doc_backend.md) as the engine.
- [Outputs](doc_outputs.md): a metrics table, NIfTI maps in subject and MNI space, and an HTML report with ITRUSST safety reference levels and a methods paragraph.
- [Cluster](doc_hpc.md) execution with SLURM or qsub, and [group](doc_group.md) reports across subjects.

A recent overview poster can be found [here](https://jkosciessa.github.io/downloads/2025-FUN25-PRESTUS.pdf).

## Versions

The `main` branch is the stable release line and the [`development` branch](https://github.com/Donders-Institute/PRESTUS/tree/development) carries new features. [Releases](https://github.com/Donders-Institute/PRESTUS/releases) are tagged with a version number; the run log prints the git hash of the checkout. Use the version selector in the page header to read the documentation of another branch.

## Citation

If you use PRESTUS for a publication, please cite:

> Kosciessa, J. Q., Chetverikov, A., Cornelissen, M., van der Zee, K., Carpino, E., Wimmers, M., Adelhöfer, N., Meijer, S., & Verhagen, L. (2026). PRESTUS (0.6.0). Zenodo. https://doi.org/10.5281/zenodo.15095860

## License and disclaimer

Released under the GNU General Public License v3.0 (see LICENSE). To help develop PRESTUS, see [Contributing](CONTRIBUTING.md).

!!! warning "Research use only"
    This software is currently under development and is provided "as is" without warranty of any kind, either express or implied, including but not limited to the implied warranties of merchantability, fitness for a particular purpose, and non-infringement. This tool is intended for research purposes only and is not designed, intended, or approved for medical or clinical use, diagnosis, or treatment of patients. The contributors accept no liability for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) arising in any way out of the use of this software. Use at your own risk.
