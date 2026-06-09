#!/usr/bin/env python3
"""
run_babelbrain_case.py — PRESTUS → BabelBrain acoustic simulation bridge

Called by babelbrain_wrapper.m via system(). Invokes RUN_SIM().RunCases()
for the requested transducer preset and writes *DataForSim.h5 output that
PRESTUS reads back via read_babelbrain_output.m.

Usage (all arguments required unless marked optional):
    python run_babelbrain_case.py \
        --basedir /path/to/cache \
        --id sub-001 \
        --target M1_left \
        --tx_system R15473 \
        --freq 300000 \
        --ppw 6 \
        --computing_backend 3 \
        --device_name Apple \
        --tx_mech_adjust_z 0.003 \
        --z_beyond_focus_mm 40 \
        --babelbrain_root /path/to/BabelBrain \
        [--weights_file /path/to/weights.h5] \
        [--do_refocusing] \
        [--dry_run]

Exit codes:
    0  — success
    1  — import or simulation error (message printed to stderr)
"""

import argparse
import sys
import os
import yaml


def parse_args():
    p = argparse.ArgumentParser(description='PRESTUS–BabelBrain bridge')
    p.add_argument('--get_version',         action='store_true',
                   help='Print version info as KEY=VALUE lines and exit. '
                        'Only --babelbrain_root is used with this flag.')
    p.add_argument('--basedir',             default='')
    p.add_argument('--id',                  default='', dest='subject_id')
    p.add_argument('--target',              default='')
    p.add_argument('--tx_system',           default='')
    p.add_argument('--freq',                default=500e3, type=float)
    p.add_argument('--ppw',                 default=6, type=int)
    p.add_argument('--computing_backend',   default=3, type=int,
                   help='1=CUDA 2=OpenCL 3=Metal 4=CPU')
    p.add_argument('--device_name',         default='Apple')
    p.add_argument('--tx_mech_adjust_z',    default=0.0, type=float,
                   help='Transducer z offset into skin [m] (default 0: geometry already in NIfTI)')
    p.add_argument('--z_steering',          default=0.0, type=float,
                   help='ZSteering: axial offset of electronic focus from preset FocalLength [m]. '
                        'Derived automatically from focal_distance_mm vs FocalLength when no '
                        'calibrated weights file is provided.')
    p.add_argument('--z_beyond_focus_mm',   default=40.0, type=float)
    p.add_argument('--babelbrain_root',     default='',
                   help='Path to BabelBrain source root (contains TranscranialModeling/). '
                        'If empty, searches relative to this script and common v05 locations.')
    p.add_argument('--weights_file',        default='',
                   help='Path to OptimizedWeightsFile (.h5); omit to use geometric focusing')
    p.add_argument('--do_refocusing',       action='store_true',
                   help='Enable BabelBrain skull refocusing back-propagation')
    p.add_argument('--dry_run',             action='store_true',
                   help='Pass bDryRun=True to RunCases — validates paths and imports without GPU')
    return p.parse_args()


# Map tx_system string → (module_path, preset_folder_name)
# module_path is relative to BabelBrain root (TranscranialModeling package)
# preset_folder is under BabelBrain/BabelBrain/Babel_<name>/default.yaml
PRESET_MODULE_MAP = {
    'CTX_500':      ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'CTX500'),
    'CTX_250':      ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'CTX250'),
    'CTX_250_2ch':  ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'CTX250_2ch'),
    'R15287':       ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'R15287'),
    'R15473':       ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'R15473'),
    'R15646':       ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'R15646'),
    'DPX_500':      ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'DPX500'),
    'DPXPC_300':    ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'DPXPC300'),
    'R15148':       ('TranscranialModeling.BabelIntegrationANNULAR_ARRAY',       'R15148'),
    'DomeTx':       ('TranscranialModeling.BabelIntegrationDomeTx',              'DomeTx'),
    'H317':         ('TranscranialModeling.BabelIntegrationH317',                'H317'),
    'H246':         ('TranscranialModeling.BabelIntegrationH246',                'H246'),
    'H301':         ('TranscranialModeling.BabelIntegrationH301',                'H301'),
    'Single':       ('TranscranialModeling.BabelIntegrationSingle',              'SingleTx'),
    'BSonix':       ('TranscranialModeling.BabelIntegrationSingle',              'SingleTx'),
    'REMOPD':       ('TranscranialModeling.BabelIntegrationREMOPD',              'REMOPD'),
    'I12378':       ('TranscranialModeling.BabelIntegrationI12378',              'I12378'),
    'IGT64_500':    ('TranscranialModeling.BabelIntegrationIGT64_500',           'IGT64_500'),
    'ATAC':         ('TranscranialModeling.BabelIntegrationATAC',                'ATAC'),
}


def find_babelbrain_root(hint=''):
    """Locate BabelBrain source root (must contain TranscranialModeling/)."""
    candidates = []
    if hint:
        candidates.append(hint)

    # This script lives in functions/babelbrain/; the BabelBrain submodule is at
    # ../../external/BabelBrain relative to this file (i.e. PRESTUS/external/BabelBrain).
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(script_dir, '..', '..', 'external', 'BabelBrain'))

    for c in candidates:
        c = os.path.normpath(c)
        if os.path.isdir(os.path.join(c, 'TranscranialModeling')):
            return c
    raise FileNotFoundError(
        "Cannot locate BabelBrain source root (needs TranscranialModeling/ subdirectory).\n"
        "Ensure the submodule is initialised:\n"
        "  git submodule update --init external/BabelBrain\n"
        "Or pass --babelbrain_root explicitly."
    )


def load_preset_config(babelbrain_root, tx_system):
    """Load default.yaml for the given preset; return dict or {}."""
    _, folder = PRESET_MODULE_MAP.get(tx_system, (None, None))
    if not folder:
        return {}
    yaml_path = os.path.join(babelbrain_root, 'BabelBrain', f'Babel_{folder}', 'default.yaml')
    if not os.path.isfile(yaml_path):
        print(f'[run_babelbrain_case] Warning: preset config not found: {yaml_path}')
        return {}
    with open(yaml_path) as f:
        cfg = yaml.safe_load(f)
    print(f'[run_babelbrain_case] Loaded preset config: {yaml_path}')
    return cfg


def load_run_sim(babelbrain_root, tx_system):
    """Import the correct RUN_SIM class for the given preset name."""
    if tx_system not in PRESET_MODULE_MAP:
        raise ImportError(
            f"Unknown tx_system '{tx_system}'. "
            f"Known presets: {sorted(PRESET_MODULE_MAP.keys())}"
        )
    module_path, _ = PRESET_MODULE_MAP[tx_system]

    if babelbrain_root not in sys.path:
        sys.path.insert(0, babelbrain_root)

    import importlib
    mod = importlib.import_module(module_path)
    cls = getattr(mod, 'RUN_SIM', None)
    if cls is None:
        raise ImportError(f"Module '{module_path}' has no RUN_SIM class.")
    return cls


def convert_weights_if_needed(weights_path):
    """Convert PRESTUS split-format weights HDF5 to BabelBrain native complex format.

    PRESTUS export_babelbrain_weights.m writes:
        /CALIBRATION/real  [N x 1] float32
        /CALIBRATION/imag  [N x 1] float32
    (MATLAB cannot write native complex HDF5 datasets.)

    BabelBrain ReadFromH5py expects:
        /CALIBRATION       [N] complex64  (single dataset at root level)

    If the file already has the BabelBrain layout (complex dataset at /CALIBRATION),
    it is returned unchanged. Otherwise a converted copy is written alongside the
    original with suffix _bb.h5 and its path is returned.
    """
    import h5py
    import numpy as np

    with h5py.File(weights_path, 'r') as f:
        if 'CALIBRATION' not in f:
            raise KeyError(f"No 'CALIBRATION' key in {weights_path}")
        cal = f['CALIBRATION']
        # BabelBrain format: dataset (not group), complex dtype
        if isinstance(cal, h5py.Dataset):
            if np.issubdtype(cal.dtype, np.complexfloating):
                print('[run_babelbrain_case] Weights already in BabelBrain complex format.')
                return weights_path
            # Single real dataset — treat as amplitude-only (unusual, pass through)
            print('[run_babelbrain_case] Warning: CALIBRATION is real-valued; using as-is.')
            return weights_path
        # PRESTUS split format: group with 'real' and 'imag' datasets
        if isinstance(cal, h5py.Group) and 'real' in cal and 'imag' in cal:
            w_real = cal['real'][:]
            w_imag = cal['imag'][:]
        else:
            raise ValueError(
                f"Unrecognised CALIBRATION layout in {weights_path}. "
                "Expected complex dataset or group with 'real'/'imag' datasets."
            )

    weights_complex = (w_real + 1j * w_imag).astype(np.complex64).ravel()

    out_path = weights_path.replace('.h5', '_bb.h5')
    if os.path.isfile(out_path):
        os.remove(out_path)

    with h5py.File(out_path, 'w') as f:
        f.create_dataset('CALIBRATION', data=weights_complex)
        # Copy metadata datasets if present
        try:
            with h5py.File(weights_path, 'r') as fsrc:
                for key in fsrc.keys():
                    if key != 'CALIBRATION':
                        fsrc.copy(key, f)
        except Exception:
            pass  # metadata is optional

    print(f'[run_babelbrain_case] Converted PRESTUS weights → BabelBrain format: {out_path}')
    return out_path


def get_babelbrain_version(bb_root):
    """Return version info dict. Prints KEY=VALUE lines to stdout for MATLAB to parse."""
    import subprocess

    version_txt = ''
    version_file = os.path.join(bb_root, 'BabelBrain', 'version.txt')
    if os.path.isfile(version_file):
        with open(version_file) as f:
            version_txt = f.read().strip()

    git_describe = ''
    git_commit   = ''
    try:
        git_describe = subprocess.check_output(
            ['git', '-C', bb_root, 'describe', '--tags', '--always'],
            stderr=subprocess.DEVNULL, text=True).strip()
        git_commit = subprocess.check_output(
            ['git', '-C', bb_root, 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        pass

    return {
        'BABELBRAIN_VERSION_TXT':      version_txt,
        'BABELBRAIN_GIT_DESCRIBE':     git_describe,
        'BABELBRAIN_GIT_COMMIT':       git_commit,
        'BABELBRAIN_ROOT':             bb_root,
    }


def main():
    args = parse_args()

    try:
        bb_root = find_babelbrain_root(args.babelbrain_root)
    except FileNotFoundError as e:
        print(f'[run_babelbrain_case] {e}', file=sys.stderr)
        sys.exit(1)
    print(f'[run_babelbrain_case] BabelBrain root: {bb_root}')

    # --get_version: emit KEY=VALUE lines and exit (used by babelbrain_wrapper.m)
    if args.get_version:
        info = get_babelbrain_version(bb_root)
        for k, v in info.items():
            print(f'{k}={v}')
        sys.exit(0)

    try:
        RunSimClass = load_run_sim(bb_root, args.tx_system)
    except ImportError as e:
        print(f'[run_babelbrain_case] Import error: {e}', file=sys.stderr)
        sys.exit(1)

    # Load preset geometry from default.yaml
    preset_cfg = load_preset_config(bb_root, args.tx_system)

    # Build geometry kwargs for RunCases (passed through to BabelFTD_Simulations)
    geom_kw = {}
    if preset_cfg:
        if 'TxDiam' in preset_cfg:
            geom_kw['Aperture']     = float(preset_cfg['TxDiam'])
        if 'FocalLength' in preset_cfg:
            geom_kw['FocalLength']  = float(preset_cfg['FocalLength'])
        if 'InDiameters' in preset_cfg:
            geom_kw['InDiameters']  = [float(x) for x in preset_cfg['InDiameters']]
        if 'OutDiameters' in preset_cfg:
            geom_kw['OutDiameters'] = [float(x) for x in preset_cfg['OutDiameters']]

    if geom_kw:
        print(f'[run_babelbrain_case] Preset geometry: '
              f"Aperture={geom_kw.get('Aperture','?')} m  "
              f"FocalLength={geom_kw.get('FocalLength','?')} m  "
              f"n_rings={len(geom_kw.get('InDiameters',[]))}")

    try:
        R = RunSimClass()
    except TypeError:
        R = RunSimClass(TxSystem=args.tx_system)

    # Build RunCases kwargs
    kw = dict(
        targets        = [args.target],
        ID             = args.subject_id,
        basedir        = args.basedir,
        Frequencies    = [args.freq],
        basePPW        = [args.ppw],
        COMPUTING_BACKEND           = args.computing_backend,
        deviceName                  = args.device_name,
        TxMechanicalAdjustmentZ     = args.tx_mech_adjust_z,
        ZSteering                   = args.z_steering,
        bUseRayleighForWater        = True,
        bTightNarrowBeamDomain      = True,
        zLengthBeyonFocalPointWhenNarrow = args.z_beyond_focus_mm / 1e3,  # mm → m
        bDoRefocusing               = args.do_refocusing,
        bDisplay                    = False,
        bDryRun                     = args.dry_run,
        **geom_kw,
    )

    if args.weights_file and os.path.isfile(args.weights_file):
        kw['OptimizedWeightsFile'] = convert_weights_if_needed(args.weights_file)
        print(f'[run_babelbrain_case] Using weights: {kw["OptimizedWeightsFile"]}')
    else:
        print('[run_babelbrain_case] No weights file — using geometric focusing.')

    if args.dry_run:
        print('[run_babelbrain_case] DRY RUN — no GPU computation will occur.')

    print(f'[run_babelbrain_case] Starting RunCases for '
          f'{args.subject_id}/{args.target} '
          f'({args.tx_system}, {args.freq/1e3:.0f} kHz, PPW={args.ppw})')

    try:
        output_files = R.RunCases(**kw)
    except Exception as e:
        print(f'[run_babelbrain_case] RunCases error: {e}', file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if output_files:
        for f in output_files:
            print(f'[run_babelbrain_case] Output: {f}')
    elif args.dry_run:
        print('[run_babelbrain_case] Dry run completed (no output files expected).')
    else:
        print('[run_babelbrain_case] Warning: RunCases returned no output files.',
              file=sys.stderr)
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
