# BabelBrain — external submodule

This directory is the expected location for BabelBrain as a git submodule.
PRESTUS uses BabelBrain as an optional acoustic simulation backend via
`functions/babelbrain/babelbrain_wrapper.m`.

## Adding the submodule

```bash
cd <PRESTUS root>
git submodule add https://github.com/ProteusMRIgHIFU/BabelBrain external/BabelBrain
git submodule update --init --recursive
```

## Python environment

BabelBrain requires Python 3.9+ with GPU compute libraries. The recommended
setup is a dedicated conda environment:

```bash
conda create -n babelbrain python=3.10
conda activate babelbrain
cd external/BabelBrain
pip install -e .
```

For Apple Silicon (Metal backend):
```bash
pip install torch torchvision torchaudio   # Metal-accelerated via MPS
```

Set `parameters.babelbrain.python_env` in your PRESTUS config to the
environment path, or leave empty to use the system `python3`.

## PRESTUS integration

The bridge entry point is `run_babelbrain_case.py` (this directory).
It is called automatically by `babelbrain_wrapper.m`.

Activate the BabelBrain backend in a PRESTUS config:

```yaml
acoustic_backend: babelbrain

babelbrain:
  tx_system: CTX_500          # must match a BabelBrain preset name
  computing_backend: 3        # 1=CUDA  2=OpenCL  3=Metal  4=CPU
  device_name: Apple
  ppw: 6
  do_refocusing: false
  z_beyond_focus_mm: 40
  python_env: /path/to/conda/envs/babelbrain   # or leave empty for system python3
  tx_mech_adjust_z_m: 0.003
```

The BabelBrain acoustic result is then passed directly to PRESTUS thermal
simulation in the same way as k-Wave output.

## Supported transducer presets

| Preset | Description |
|---|---|
| `CTX_500` | Sonic Concepts CTX-500, 7-ring annular, 500 kHz |
| `CTX_250` | CTX-250, 7-ring annular, 250 kHz |
| `CTX_250_2ch` | CTX-250, 2-channel variant |
| `R15287` | Sonic Concepts R15287, 4-ring, 500 kHz |
| `R15473` | Sonic Concepts R15473, 8-ring, 500 kHz |
| `R15646` | Sonic Concepts R15646, 10-ring, 500 kHz |
| `DPX_500` | Dual-frequency 5-element annular, 500 kHz |
| `H317` | Single-element bowl, 500 kHz |
| `H246` | Single-element bowl, 500 kHz (smaller aperture) |
| `H301` | Single-element bowl, 1 MHz |
| `IGT64_500` | IGT 64-channel phased array, 500 kHz |
| `Single` / `BSonix` | Generic single elements |

For a PRESTUS transducer not in this list, `babelbrain_transducer_map.m`
will generate a custom element CSV. Follow the instructions it prints to
register the transducer as a BabelBrain preset.

## Version compatibility

Tested against BabelBrain commit range as documented in
`/Volumes/2425122.01/v05/doc/notes_babelbrain_prestus_integration.md`.
The `RUN_SIM().RunCases()` API is the stable programmatic interface;
internal BabelBrain changes to `CalculateFieldProcess.py` (the GUI layer)
do not affect this integration.
