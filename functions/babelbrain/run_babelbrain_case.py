#!/usr/bin/env python3
"""
run_babelbrain_case.py — PRESTUS → BabelBrain acoustic simulation bridge

Called by babelbrain_wrapper.m via run_command() (java.lang.ProcessBuilder,
no shell). Invokes RUN_SIM().RunCases()
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
import importlib
import sys
import os
import yaml
import numpy as np


def _read_proc_status_kb(key):
    with open('/proc/self/status') as f:
        for line in f:
            if line.startswith(key + ':'):
                return int(line.split()[1])
    return None


def _cgroup_memory_files():
    """Return (current, peak) file paths of this process's memory cgroup, or (None, None).

    Under SLURM the cgroup also holds the parent MATLAB process, so its usage
    is the figure the OOM handler acts on.
    """
    try:
        with open('/proc/self/cgroup') as f:
            lines = f.read().splitlines()
    except OSError:
        return None, None
    for line in lines:                        # cgroup v1 memory controller
        _, ctrls, rel = line.split(':', 2)
        if 'memory' in ctrls.split(',') and os.path.isdir('/sys/fs/cgroup/memory'):
            d = '/sys/fs/cgroup/memory' + rel
            return (os.path.join(d, 'memory.usage_in_bytes'),
                    os.path.join(d, 'memory.max_usage_in_bytes'))
    for line in lines:                        # cgroup v2 unified hierarchy
        if line.startswith('0::'):
            d = '/sys/fs/cgroup' + line[3:]
            return os.path.join(d, 'memory.current'), os.path.join(d, 'memory.peak')
    return None, None


def _read_bytes_gb(path):
    try:
        with open(path) as f:
            return int(f.read().strip()) / 1024**3
    except (OSError, ValueError, TypeError):
        return float('nan')


def _start_memory_monitor(interval_s=5.0, step_gb=0.5):
    """Log this process's peak RSS and the cgroup usage whenever the peak grows.

    An OOM kill (SIGKILL) skips atexit handlers, so the monitor prints
    incrementally: the last '[memory_bb]' line before a kill is the high-water
    mark. The tag keeps these lines distinct from MATLAB's '[memory]' lines in
    the pipeline/SLURM logs (parsed by parse_memory_tags.m). Linux only;
    silently does nothing elsewhere.
    """
    if not os.path.isfile('/proc/self/status'):
        return
    import atexit
    import threading
    import time
    cg_cur, cg_peak = _cgroup_memory_files()
    state = {'last_logged_gb': 0.0}

    def _line(tag):
        rss = (_read_proc_status_kb('VmRSS') or 0) / 1024**2
        hwm = (_read_proc_status_kb('VmHWM') or 0) / 1024**2
        return hwm, (f'[memory_bb] {tag}: python rss={rss:.1f} GB (peak {hwm:.1f} GB) | '
                     f'cgroup current={_read_bytes_gb(cg_cur):.1f} GB, '
                     f'peak={_read_bytes_gb(cg_peak):.1f} GB')

    def _loop():
        t0 = time.time()
        while True:
            try:
                hwm, msg = _line(f't={time.time() - t0:.0f}s')
                if hwm - state['last_logged_gb'] >= step_gb:
                    state['last_logged_gb'] = hwm
                    print(msg, flush=True)
            except Exception:
                return                       # diagnostics must never break the run
            time.sleep(interval_s)

    threading.Thread(target=_loop, name='memory-monitor', daemon=True).start()
    atexit.register(lambda: print(_line('exit')[1], flush=True))


class _TeeStream:
    """Write to multiple streams simultaneously (console + log file)."""
    def __init__(self, *streams):
        self.streams = streams
    def write(self, data):
        for s in self.streams:
            s.write(data)
    def flush(self):
        for s in self.streams:
            s.flush()
    def fileno(self):
        return self.streams[0].fileno()


def resolve_device_name(computing_backend, device_name):
    """Return the first matching device name for the given backend.

    When device_name is the sentinel default ('CPU' or 'auto') but the backend
    is CUDA (1) or OpenCL (2), auto-detect the first available device so the
    user does not need to know the exact hardware string.
    """
    _is_sentinel = device_name in ('CPU', 'auto')
    if computing_backend == 1 and _is_sentinel:
        # Use nvidia-smi for device detection — avoids CUDA library version mismatches
        try:
            import subprocess
            out = subprocess.check_output(
                ['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'],
                stderr=subprocess.STDOUT
            ).decode().strip()
            lines = [l.strip() for l in out.splitlines() if l.strip()]
            if lines:
                name = lines[0]
                print(f'[run_babelbrain_case] Auto-detected CUDA device: {name}')
                return name
            else:
                print('[run_babelbrain_case] ERROR: nvidia-smi found no GPUs. '
                      'Ensure the SLURM job requests a GPU (--gres=gpu:1) or set device_name explicitly.',
                      file=sys.stderr)
                sys.exit(1)
        except Exception as e:
            print(f'[run_babelbrain_case] nvidia-smi device detection failed: {e}', file=sys.stderr)
            sys.exit(1)
    elif computing_backend == 2 and _is_sentinel:
        try:
            import pyopencl as cl
            platforms = cl.get_platforms()
            for pf in platforms:
                devices = pf.get_devices(cl.device_type.GPU)
                if devices:
                    name = devices[0].name.strip()
                    print(f'[run_babelbrain_case] Auto-detected OpenCL device: {name}')
                    return name
            print('[run_babelbrain_case] ERROR: no OpenCL GPU found. '
                  'Set --device_name explicitly or use a different --computing_backend.',
                  file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f'[run_babelbrain_case] OpenCL device detection failed: {e}', file=sys.stderr)
            sys.exit(1)
    return device_name


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
    p.add_argument('--computing_backend',   default=4, type=int,
                   help='1=CUDA 2=OpenCL 3=Metal 4=CPU (default)')
    p.add_argument('--device_name',         default='CPU',
                   help='Hardware device string passed to BabelBrain. '
                        'For CUDA/OpenCL, omit to auto-detect via nvidia-smi / pyopencl.')
    p.add_argument('--tx_mech_adjust_z',    default=0.0, type=float,
                   help='Transducer z offset into skin [m] (default 0: geometry already in NIfTI)')
    p.add_argument('--aperture',            default=0.0, type=float,
                   help='Transducer aperture (outer diameter) [m]. Required with --tx_system PRESTUS '
                        'and by the single-element presets (Single, BSonix), whose default.yaml has '
                        'no geometry; overrides the preset TxDiam when > 0.')
    p.add_argument('--focal_length',        default=0.0, type=float,
                   help='Transducer focal length (radius of curvature) [m]. Required with --tx_system '
                        'PRESTUS and by the single-element presets; overrides the preset FocalLength '
                        'when > 0.')
    p.add_argument('--in_diameters',        default='',
                   help='--tx_system PRESTUS: comma-separated inner ring diameters [m].')
    p.add_argument('--out_diameters',       default='',
                   help='--tx_system PRESTUS: comma-separated outer ring diameters [m].')
    p.add_argument('--z_steering',          default=0.0, type=float,
                   help='ZSteering: axial offset of electronic focus from preset FocalLength [m] '
                        '(default 0; babelbrain_wrapper.m derives it).')
    p.add_argument('--z_beyond_focus_mm',      default=40.0,  type=float,
                   help='Axial domain extent beyond focal point [mm] (default 40). '
                        'Increase to reduce boundary reflections at the focus.')
    p.add_argument('--no_tight_narrow_beam',   action='store_true',
                   help='Disable bTightNarrowBeamDomain (use full BabelBrain domain). '
                        'Slower but eliminates lateral truncation artefacts.')
    p.add_argument('--babelbrain_root',     default='',
                   help='Path to BabelBrain source root (contains TranscranialModeling/). '
                        'If empty, uses external/BabelBrain of this PRESTUS copy.')
    p.add_argument('--weights_file',        default='',
                   help='Path to OptimizedWeightsFile (.h5); omit to use geometric focusing')
    p.add_argument('--do_refocusing',       action='store_true',
                   help='Enable BabelBrain skull refocusing back-propagation')
    p.add_argument('--dry_run',             action='store_true',
                   help='Pass bDryRun=True to RunCases — validates paths and imports without GPU')
    p.add_argument('--water_only',          action='store_true',
                   help='Run water-only baseline simulation (bWaterOnly=True). '
                        'Produces *_Water_DataForSim.h5 alongside the tissue DataForSim.h5. '
                        'Required by BabelBrain thermal for ISPPA normalization.')
    p.add_argument('--no_rayleigh_for_water', action='store_true',
                   help='Compute the water-only field with the full-wave FDTD solver instead of '
                        'the Rayleigh–Sommerfeld integral (sets bUseRayleighForWater=False). Use '
                        'this to compare against other full-wave backends (k-Wave/jwave) on a like-'
                        'for-like basis; the Rayleigh default focuses slightly tighter/deeper and '
                        'so reads as a rigid axial offset + higher peak vs full-wave solvers.')
    p.add_argument('--use_ct',              action='store_true',
                   help='Activate BabelBrain CT path (bUseCT=True). Requires a CT.nii.gz '
                        'alongside the BabelViscoInput mask. Use when passing raw HU values '
                        'from a UTE/PETRA scan for BabelBrain to map to acoustic properties.')
    p.add_argument('--ct_mapping',          default='Webb-Marsac',
                   help='HU-to-property mapping method (MappingMethod in __init__). '
                        'Valid: Webb-Marsac (default), Aubry, Pichardo, McDannold, '
                        'Marsac-Aubry, Pichardo-Marsac, McDannold-Marsac. '
                        'Only used when --use_ct is set.')
    p.add_argument('--air_regions',         action='store_true',
                   help='CT path only: write <target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz from the raw-HU '
                        'CT.nii.gz (BabelBrain recipe: HU window, 3-voxel median filter, exterior air '
                        'removed) and enable bExtractAirRegions, so BabelBrain treats internal air '
                        '(sinuses, mastoid cells) as a pressure-release reflector. Ignored with --water_only.')
    p.add_argument('--air_hu_window',       nargs=2, type=float, default=[-1200.0, -400.0],
                   metavar=('LO', 'HI'),
                   help='HU range classified as air for --air_regions (default -1200 -400, as BabelBrain).')
    p.add_argument('--air_keep_boundary',   action='store_true',
                   help='--air_regions: keep air components that touch the domain boundary. By default '
                        'they are removed with the largest component, since cropping can split the '
                        'exterior air (and the coupling gap in front of the transducer) into pieces.')
    p.add_argument('--use_petra',           action='store_true',
                   help="Pass bPETRA=True to __init__ -- uses Bjorn's density formula for "
                        'PETRA-sequence-derived pseudo-CT. Leave unset for ZTE/UTE pCT.')
    p.add_argument('--no_skull_atten_fix',  action='store_true',
                   help='Disable the non-CT skull-attenuation fix. By default, on the non-CT '
                        '(binary mask) path PRESTUS reproduces BabelBrain CT-mode bone treatment '
                        '(ShearSoS=0) and sets the bone Q correction (default 1.0, see --q_correction) '
                        'so the injected bone attenuation is realised as specified. '
                        'Set this flag only to reproduce the un-patched BabelBrain behaviour.')
    p.add_argument('--no_zdomain_fix',      action='store_true',
                   help='Disable the ZDomainStart neutralization fix. By default '
                        'CalculateDomainZReference is patched to return 0 so the transducer is '
                        'placed at exactly the requested TxMechanicalAdjustmentZ. Without the fix, '
                        'BabelBrain silently adds a skin-surface-dependent domain correction '
                        '(GlobalStartSkin - LineOfSightStartSkin) on top of the requested '
                        'adjustment, moving the source millimetres off the PRESTUS-intended '
                        'position (e.g. -5.5 mm). Set this flag only to reproduce the '
                        'un-patched BabelBrain behaviour.')
    p.add_argument('--medium_properties_json', default='',
                   help='Path to a JSON dump of PRESTUS parameters.medium_properties. When given, '
                        "BabelBrain's built-in non-CT material database (MatFreq) is overwritten "
                        'with these values (sound_speed, density, attenuation converted from '
                        'dB/cm/MHz to Np/m; shear set to 0) so k-Wave / jwave / BabelBrain assume '
                        'identical tissue acoustics. No effect on the CT path bone mapping.')
    p.add_argument('--bone_shear',          action='store_true',
                   help='Elastic bone (babelbrain.shear: prestus_elastic). Keeps BabelBrain shear on for '
                        'the non-CT path and sets cortical/trabecular shear speed, shear attenuation and '
                        'longitudinal attenuation from the babelbrain_elastic entries of the skull layers '
                        'in --medium_properties_json (the kwave_elastic values at the drive frequency). '
                        'Requires --medium_properties_json; not available with --use_ct.')
    p.add_argument('--log_file',            default='',
                   help='Path to write a copy of all stdout/stderr output for post-hoc debugging.')
    p.add_argument('--source_amp_pa',       default=0.0, type=float,
                   help='Source pressure amplitude [Pa] derived from calibrated particle velocity '
                        '(velocity_m_s * rho_water * c_water). When provided, overrides '
                        "BabelBrain's default of 50 kPa so simulated ISPPA matches calibration.")
    p.add_argument('--q_correction',        default=0.0, type=float,
                   help='Bone Q-factor correction for the relaxation attenuation model '
                        '(default 0 = 1.0 on the non-CT path, BabelBrain default 3.0 on the CT path). '
                        'Realised bone attenuation scales ~1/Q: the bone-slab phantom (6 mm cortical, '
                        '300 kHz, PPW 6) realises 91%% of the injected alpha at Q=1 (k-Wave: 91%%) and '
                        "29%% at BabelBrain's CT-mode 3.0. Values > 0 override both paths.")
    p.add_argument('--debug',               action='store_true',
                   help='Keep BabelBrain debug-only outputs that PRESTUS itself does not read, e.g. '
                        '*_Sub_NORM.nii.gz (brain-masked, normalized FullElasticSolution — useful to '
                        'see where the brain/skull mask boundary falls relative to the solved field).')
    return p.parse_args()


# Preset names BabelBrain has renamed; old configs keep working.
TX_ALIASES = {'REMOPD': 'TheraPhas'}

# Reserved tx_system: the PRESTUS transducer's own ring geometry (--aperture,
# --focal_length, --in_diameters, --out_diameters) on BabelBrain's generic
# focused annular template, which BabelBrain also uses for custom transducers.
PRESTUS_GEOMETRY = 'PRESTUS'

# BabelBrain package modules, relative to its source root.
BB_BASE    = 'TranscranialModeling.babel_integration.integration_templates.babel_integration_base'
BB_HELPERS = 'TranscranialModeling.babel_integration.babel_integration_helpers'
BB_ANNULAR = 'TranscranialModeling.babel_integration.integration_templates.babel_integration_focused_annular_array'
BB_TXGEOM  = 'TranscranialModeling.tx_geometries'


def add_babelbrain_paths(bb_root):
    """Put the BabelBrain source root first on sys.path and its BabelBrain/ app
    folder last (for Utils, ArtifactIO), as BabelBrain's own offline examples do.
    The root must come first so `import BabelBrain` finds the package."""
    if bb_root not in sys.path:
        sys.path.insert(0, bb_root)
    app_dir = os.path.join(bb_root, 'BabelBrain')
    if app_dir not in sys.path:
        sys.path.append(app_dir)


def resolve_preset(bb_root, tx_system):
    """Look a preset up in BabelBrain's transducer registry.

    Returns a dict with the BabelBrain preset name, the RUN_SIM module path, the
    preset default.yaml path and the transducer_type (geometry family), derived
    the same way BabelBrain's CalculateFieldProcess does.
    """
    if tx_system == PRESTUS_GEOMETRY:
        return {'name': PRESTUS_GEOMETRY, 'type': 'focused_annular_array',
                'module': BB_ANNULAR, 'yaml': ''}
    name = TX_ALIASES.get(tx_system, tx_system)
    reg_path = os.path.join(bb_root, 'BabelBrain', 'SelFiles', 'transducer_list.yaml')
    with open(reg_path) as f:
        registry = {e['name']: e for e in yaml.safe_load(f) if not e.get('custom', False)}
    if name not in registry:
        raise ImportError(f"Unknown tx_system '{tx_system}'. Known presets: {sorted(registry)}")
    entry = registry[name]
    ttype = entry['transducer_type']
    return {
        'name':   name,
        'type':   ttype,
        'module': f'TranscranialModeling.babel_integration.{ttype}.babel_integration_{name}',
        'yaml':   os.path.join(bb_root, 'BabelBrain', 'babel_transducers', ttype,
                               entry['module_name'], 'default.yaml'),
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


def load_preset_config(preset):
    """Load default.yaml for the given preset (resolve_preset dict); return dict or {}."""
    yaml_path = preset['yaml']
    tx_system = preset['name']
    if not yaml_path:
        return {}   # PRESTUS geometry: everything comes from the command line
    if not os.path.isfile(yaml_path):
        print(f'[run_babelbrain_case] WARNING: preset config not found: {yaml_path}', file=sys.stderr)
        print(f'[run_babelbrain_case] WARNING: BabelBrain will use hardcoded geometry defaults '
              f'for {tx_system}. Verify that transducer aperture, focal length, and ring '
              f'diameters are correct for your hardware.', file=sys.stderr)
        return {}
    with open(yaml_path) as f:
        cfg = yaml.safe_load(f)
    print(f'[run_babelbrain_case] Loaded preset config: {yaml_path}')
    return cfg


def load_run_sim(preset):
    """Import the RUN_SIM class for the given preset (resolve_preset dict)."""
    import importlib
    mod = importlib.import_module(preset['module'])
    cls = getattr(mod, 'RUN_SIM', None)
    if cls is None:
        raise ImportError(f"Module '{preset['module']}' has no RUN_SIM class.")
    return cls


def convert_weights_if_needed(weights_path):
    """Convert legacy PRESTUS split-format weights HDF5 to BabelBrain native complex format.

    export_babelbrain_weights.m writes the native layout (compound r/i
    /CALIBRATION with 'type' attributes), which passes through unchanged.
    The older split layout has:
        /CALIBRATION/real  [N x 1] float32
        /CALIBRATION/imag  [N x 1] float32

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
        # Carry the raw drive (for the AllPhi-referenced correction) through the
        # conversion — SaveToH5py below rewrites the file from scratch, so these
        # would otherwise be lost and the correction would silently fall back.
        drive_phases = np.asarray(f['DRIVE/phases_rad']).ravel() if ('DRIVE' in f and 'phases_rad' in f['DRIVE']) else None
        drive_amp    = np.asarray(f['DRIVE/amp_ratio']).ravel()  if ('DRIVE' in f and 'amp_ratio'  in f['DRIVE']) else None

    weights_complex = (w_real + 1j * w_imag).astype(np.complex64).ravel()

    out_path = weights_path.replace('.h5', '_bb.h5')
    if os.path.isfile(out_path):
        os.remove(out_path)

    # Use SaveToH5py so each dataset gets the 'type' attribute that ReadFromH5py requires.
    from BabelViscoFDTD.H5pySimple import SaveToH5py
    data_dict = {'CALIBRATION': weights_complex}
    if drive_phases is not None and drive_amp is not None:
        data_dict['DRIVE_phases_rad'] = drive_phases.astype(np.float64)
        data_dict['DRIVE_amp_ratio']  = drive_amp.astype(np.float64)
    SaveToH5py(data_dict, out_path)

    print(f'[run_babelbrain_case] Converted PRESTUS weights → BabelBrain format: {out_path}')

    # Read-back verification: confirm /CALIBRATION is a complex dataset that BabelBrain can open
    with h5py.File(out_path, 'r') as _f:
        _cal = _f.get('CALIBRATION')
        if _cal is None:
            raise ValueError(
                f'Weights verification failed: /CALIBRATION missing from converted file {out_path}.'
            )
        if not isinstance(_cal, h5py.Dataset):
            raise ValueError(
                f'Weights verification failed: /CALIBRATION in {out_path} is a group, not a dataset.'
            )
        if not np.issubdtype(_cal.dtype, np.complexfloating):
            raise ValueError(
                f'Weights verification failed: /CALIBRATION dtype={_cal.dtype} in {out_path}; '
                'expected complex64. BabelBrain ReadFromH5py will reject this file.'
            )
        _n = _cal.size
    print(f'[run_babelbrain_case] Weights verified: {_n} elements, dtype=complex64')
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


def _patch_annular_array_weights(RunSimClass, dump_path=None, anchor_path=None,
                                 tx_mech_adjust_z=None):
    """Monkey-patch AdjustWeightAmplitudes on the SimulationConditions class of
    annular-array modules to handle per-ring calibration weights (N_rings values)
    by expanding each ring weight to all its mesh subelements.  Per-subelement
    weights (N_subelems values) pass through unchanged.

    When dump_path is given (debug mode), the actual source BabelBrain built —
    sub-element positions, area weights, per-ring geometric steer (AllPhi),
    applied weights, and the preset geometry — is written to an HDF5 file the
    first time AdjustWeightAmplitudes runs, for offline comparison against the
    ideal (analytical) source geometry."""
    import numpy as np
    import importlib

    # Resolve the preset's SimulationConditions class (RUN_SIM → BabelFTD_Simulations
    # → SimulationConditions) and patch only focused annular arrays.
    sim_cls = getattr(RunSimClass, '_BabelFTDSimClass', None)
    SimCond = getattr(sim_cls, '_SimConditionsClass', None)
    if SimCond is None:
        print('[run_babelbrain_case] Warning: could not find SimulationConditions — '
              'per-ring weight expansion patch not applied.')
        return
    if not issubclass(SimCond, importlib.import_module(BB_ANNULAR).SimulationConditions):
        return

    def _dump_source(self, ret):
        if not dump_path or getattr(self, '_prestus_source_dumped', False):
            return
        try:
            import h5py, os
            os.makedirs(os.path.dirname(dump_path), exist_ok=True)
            bpp = getattr(self, 'BasePhasedArrayProgramming', None)
            allphi = np.angle(np.asarray(bpp).ravel()) if bpp is not None else np.array([])
            with h5py.File(dump_path, 'w') as f:
                # Sub-element source geometry BabelBrain actually integrates over.
                f['center']   = np.asarray(self._Tx['center'], dtype=np.float64)   # [N_sub x 3] positions [m]
                f['ds']       = np.asarray(self._Tx['ds'], dtype=np.float64).ravel()  # [N_sub] surface-area weights [m^2]
                f['elemdims'] = np.asarray([int(e[0]) for e in self._Tx['elemdims']], dtype=np.int64)  # subelems/ring
                if 'elemcenter' in self._Tx:
                    f['elemcenter'] = np.asarray(self._Tx['elemcenter'], dtype=np.float64)  # [N_rings x 3]
                f['AllPhi_rad']  = allphi.astype(np.float64)                          # per-ring geometric steer
                f['weights_ret'] = np.asarray(ret, dtype=np.complex64).ravel()        # applied per-subelem weight
                for attr, key in [('_SourceAmpPa', 'SourceAmpPa'), ('_SpatialStep', 'SpatialStep'),
                                  ('_Frequency', 'Frequency'), ('_FocalLength', 'FocalLength'),
                                  ('_Aperture', 'Aperture')]:
                    v = getattr(self, attr, None)
                    if v is not None:
                        try: f[key] = float(v)
                        except Exception: pass
                for attr, key in [('_InDiameters', 'InDiameters'), ('_OutDiameters', 'OutDiameters')]:
                    v = getattr(self, attr, None)
                    if v is not None:
                        f[key] = np.asarray(v, dtype=np.float64).ravel()
            self._prestus_source_dumped = True
            print(f'[run_babelbrain_case] Source dump written: {dump_path} '
                  f'(N_sub={np.asarray(self._Tx["center"]).shape[0]}, N_rings={self._Tx["NumberElems"]})')
        except Exception as e:
            print(f'[run_babelbrain_case] Source dump failed (non-fatal): {e}')

    def _write_anchor_sidecar(self):
        # Source-placement DIAGNOSTIC sidecar. The output resample in
        # read_babelbrain_output.m anchors on the /TargetLocation tissue marker
        # (the label-5 voxel PRESTUS wrote at focus_pos), which is the exact
        # inverse of the forward mask transform — so tissue lands back on the
        # original medium_masks by construction and no anchor correction is
        # applied from this file. What still matters is whether BabelBrain
        # actually placed the SOURCE where PRESTUS asked: with the ZDomainStart
        # neutralization active the only remaining discrepancy is the whole-voxel
        # snap-back ("Correction") that keeps the bowl behind the FDTD source-
        # injection plane (typically ≤ 1 native voxel). We record
        #     residual = realized CoC z - (marker z + TxMechanicalAdjustmentZ)
        # i.e. realized-minus-requested source position. read_babelbrain_output.m
        # reads it back purely as a diagnostic and warns when it exceeds a voxel.
        if not anchor_path:
            return
        try:
            import json, os
            C = np.asarray(self._Tx['center'], dtype=np.float64)   # sub-elem positions [m], domain frame
            Rc = float(self._FocalLength)                            # bowl radius = focal length [m]
            fsl = getattr(self, '_FocalSpotLocation', None)
            zdim = getattr(self, '_ZDim', None)
            if fsl is None or zdim is None:
                print('[run_babelbrain_case] Anchor sidecar skipped: _FocalSpotLocation/_ZDim unavailable.')
                return
            r2 = C[:, 0]**2 + C[:, 1]**2
            coc_z    = float(np.mean(C[:, 2] + np.sqrt(np.maximum(Rc**2 - r2, 0.0))))  # realized source centre of curvature z [m]
            marker_z = float(np.asarray(zdim).ravel()[int(fsl[2])])                    # focal-spot marker z [m]
            # Requested CoC position: marker sits at preset FocalLength from the
            # unadjusted CoC, so requested CoC z = marker z + TxMechanicalAdjustmentZ
            # + ZSteering-independent geometry (CoC-at-marker when unadjusted).
            req_adjust_m = float(tx_mech_adjust_z) if tx_mech_adjust_z is not None else 0.0
            requested_coc_z = marker_z + req_adjust_m
            residual_mm = (coc_z - requested_coc_z) * 1e3
            offset_mm   = (coc_z - marker_z) * 1e3   # kept for backward compat / inspection
            neutralized = getattr(self, '_prestus_zdomain_neutralized_m', None)
            os.makedirs(os.path.dirname(anchor_path), exist_ok=True)
            with open(anchor_path, 'w') as fh:
                json.dump({'source_residual_mm': residual_mm,
                           'requested_tx_adjust_mm': req_adjust_m * 1e3,
                           'zdomain_neutralized_mm': (neutralized * 1e3) if neutralized is not None else None,
                           'anchor_offset_mm': offset_mm,
                           'source_coc_z_mm': coc_z * 1e3,
                           'focalspot_marker_z_mm': marker_z * 1e3,
                           'curv_radius_mm': Rc * 1e3}, fh)
            print(f'[run_babelbrain_case] Source placement: realized CoC {coc_z*1e3:.3f} mm, '
                  f'requested {requested_coc_z*1e3:.3f} mm (marker {marker_z*1e3:.3f} '
                  f'{req_adjust_m*1e3:+.3f} adj) -> residual {residual_mm:+.3f} mm '
                  f'-> {anchor_path}')
        except Exception as e:
            print(f'[run_babelbrain_case] Anchor sidecar write failed (non-fatal): {e}')

    def _adjust_weight_amplitudes(self):
        ret = _compute_weights(self)
        _dump_source(self, ret)
        _write_anchor_sidecar(self)
        return ret

    def _compute_weights(self):
        if self._OptimizedWeights is None:
            return np.ones((self._Tx['center'].shape[0], 1), np.complex64)

        n_sub = self._Tx['center'].shape[0]
        n_rings = self._Tx['NumberElems']

        # AllPhi-referenced correction: if the weights file carries the raw
        # intended per-ring drive (/DRIVE/phases_rad, /DRIVE/amp_ratio), rebuild
        # the per-ring weight against BabelBrain's OWN internally-computed
        # AllPhi so the realised per-ring phase equals phases_rad exactly:
        #   effective = AllPhi + angle(weight) = AllPhi + (phases_rad - AllPhi)
        #             = phases_rad.
        # BasePhasedArrayProgramming (= exp(1j*AllPhi), per ring) is set in
        # CalculateRayleighFieldsForward BEFORE AdjustWeightAmplitudes is called,
        # so AllPhi is available here. This sidesteps any analytic reconstruction
        # of BabelBrain's geometric steer (which does not reproduce the voxelised-
        # grid AllPhi).
        drive = None
        try:
            import h5py
            with h5py.File(self._OptimizedWeightsFile, 'r') as _f:
                # Converted (_bb.h5) layout uses flat keys; original PRESTUS file
                # uses a /DRIVE group.
                if 'DRIVE_phases_rad' in _f and 'DRIVE_amp_ratio' in _f:
                    drive = (np.asarray(_f['DRIVE_phases_rad']).ravel().astype(np.float64),
                             np.asarray(_f['DRIVE_amp_ratio']).ravel().astype(np.float64))
                elif 'DRIVE' in _f and 'phases_rad' in _f['DRIVE']:
                    drive = (np.asarray(_f['DRIVE/phases_rad']).ravel().astype(np.float64),
                             np.asarray(_f['DRIVE/amp_ratio']).ravel().astype(np.float64))
        except Exception as _e:
            print(f'[run_babelbrain_case] Could not read DRIVE from weights file '
                  f'({_e}); using stored CALIBRATION weights.')

        bpp = getattr(self, 'BasePhasedArrayProgramming', None)
        if drive is not None and bpp is None:
            print('[run_babelbrain_case] /DRIVE present but AllPhi (BasePhasedArrayProgramming) '
                  'not yet computed; using stored CALIBRATION weights for this call.')
        if drive is not None and bpp is not None:
            phases_rad, amp_ratio = drive
            allphi = np.angle(np.asarray(bpp).ravel())
            if phases_rad.size != n_rings or allphi.size != n_rings:
                raise ValueError(
                    f'AllPhi-referenced correction: /DRIVE has {phases_rad.size} rings and '
                    f'AllPhi has {allphi.size}, but the Tx has {n_rings} rings.')
            w_ring = amp_ratio * np.exp(1j * (phases_rad - allphi))
            print('[run_babelbrain_case] AllPhi-referenced weights: realised per-ring phase '
                  '= intended drive (deg): '
                  f'{np.round(np.rad2deg(phases_rad - phases_rad[0]), 1).tolist()}')
            w_expanded = np.empty(n_sub, dtype=np.complex64)
            base = 0
            for n in range(n_rings):
                count = self._Tx['elemdims'][n][0]
                w_expanded[base:base + count] = w_ring[n]
                base += count
            return w_expanded.reshape((n_sub, 1))

        if self._OptimizedWeights.size == n_sub:
            # Already per-subelement — use directly (normal base-class path)
            return self._OptimizedWeights.reshape((n_sub, 1)).astype(np.complex64)

        if self._OptimizedWeights.size == n_rings:
            # Per-ring calibration weights — expand to per-subelement
            print(f'[run_babelbrain_case] Expanding {n_rings} per-ring weights '
                  f'to {n_sub} subelements.')
            w_expanded = np.empty(n_sub, dtype=np.complex64)
            base = 0
            for n in range(n_rings):
                count = self._Tx['elemdims'][n][0]
                w_expanded[base:base + count] = self._OptimizedWeights[n]
                base += count
            return w_expanded.reshape((n_sub, 1))

        raise ValueError(
            f"OptimizedWeightsFile {self._OptimizedWeightsFile} has "
            f"{self._OptimizedWeights.size} subelements, but the Tx has "
            f"{n_sub} elements ({n_rings} rings)."
        )

    SimCond.AdjustWeightAmplitudes = _adjust_weight_amplitudes
    print(f'[run_babelbrain_case] Patched AdjustWeightAmplitudes on SimulationConditions '
          f'for per-ring weight expansion.')


def _patch_noncat_skull_attenuation(bb_root, bone_q=1.0, keep_shear=False):
    """Make the non-CT (binary-mask) path attenuate through the skull correctly
    by reproducing BabelBrain's *CT-mode* bone treatment.

    THE BUG
    -------
    On the non-CT path, babel_integration_base.Step1_InitializeConditions builds the
    per-material Q-correction array as ``QCorrArr = np.ones(5)`` (or ``np.ones(8)``
    with brain segmentation) — i.e. it applies *no* Q correction to cortical or
    trabecular bone. The transcranial field then comes out almost identical to
    the water-only field (ISPPA ~30 W/cm² vs ~11 W/cm² for k-Wave / jwave with
    identical settings).

    The CT path adds bone with ``ShearSoS = 0`` and ``QCorrArr[bone] =
    self._QCorrection`` (default 3.0). This routine sets on the non-CT path:

      * ``self._Shear = 0``  → bone modelled as a fluid (no shear / mode
        conversion), matching the CT path and the fluid solvers (k-Wave, jwave).
        Set *before* Step1 runs so the bone AddMaterial() calls and the CFL /
        temporal-step calculation in UpdateConditions both see shear = 0.
      * ``QCorrArr[bone] = bone_q`` (default 1.0; ``--q_correction`` overrides).
        Realised bone attenuation scales ~1/Q. Bone-slab phantom (6 mm cortical,
        78.5 Np/m, 300 kHz, PPW 6; abs/noabs ratio): Q=1 realises 91% of the
        injected alpha (k-Wave 91%), Q=3 29%, Q=6 14%, Q=10 8.7%, Q=14 6.3%. The
        CT-mode 3.0 is BabelBrain's calibration for its own CT-derived bone
        attenuation and must not be applied to PRESTUS-specified alpha; the
        shear-speed setting is the operative part of the fix. Set
        *after* Step1 runs: the QCorr array is only consumed by the FDTD
        propagation (Step2+), never by the temporal-step calc, so overwriting it
        post-setup is safe.

    WHY A PATCH AND NOT A BenchmarkTestFile
    ---------------------------------------
    BabelBrain's BenchmarkTestFile mechanism could inject materials + QCorrArr, but
    it also requires a MaterialMap whose shape matches BabelBrain's *internal*
    domain exactly. With ``bTightNarrowBeamDomain`` (babelbrain.tight_narrow_beam) the domain is
    cropped to follow the acoustic beam, so that shape is emergent from the
    Rayleigh-field setup and is not known until BabelBrain has already run. A
    BenchmarkTestFile would therefore need a throw-away run to discover the grid.
    Patching the two offending values lets BabelBrain build the domain and material
    map itself, needs no pre-run, and touches nothing under external/.

    The patch is self-guarding: it only alters the non-CT, non-water, non-homogeneous,
    non-benchmark branch. It is a no-op for CT / water-only / benchmark runs.

    With ``keep_shear`` (``--bone_shear``) ``self._Shear`` is left at BabelBrain's 1.0,
    so the bone shear speed and attenuation injected into MatFreq are used; only the
    bone Q correction is set.
    """
    import inspect  # noqa: F401  (kept for parity with sibling patch helper)

    add_babelbrain_paths(bb_root)
    try:
        BabelFTD_Simulations_BASE = importlib.import_module(BB_BASE).BabelFTD_Simulations_BASE
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not import BabelFTD_Simulations_BASE; '
              f'non-CT skull-attenuation fix NOT applied ({e}).', file=sys.stderr)
        return

    if getattr(BabelFTD_Simulations_BASE, '_prestus_skull_fix', False):
        return  # already patched (idempotent)

    orig_step1 = BabelFTD_Simulations_BASE.Step1_InitializeConditions

    # Bone material indices in the non-CT layered material map:
    #   0=water 1=skin 2=cortical 3=trabecular 4=brain [5=WM 6=GM 7=CSF]
    BONE_IDX = (2, 3)

    def patched_step1(self):
        non_ct_layered = (self._CTFNAME is None
                          and not self._bWaterOnly
                          and not getattr(self, '_bForceHomogenousMedium', False)
                          and len(getattr(self, '_BenchmarkTestFile', '') or '') == 0)

        if non_ct_layered and not keep_shear:
            # Bone → fluid. Must be set before Step1's AddMaterial() calls and the
            # CFL calculation, both of which read self._Shear.
            self._Shear = 0.0

        orig_step1(self)

        if non_ct_layered:
            ss = self._SIM_SETTINGS
            qc = np.array(ss._QCorrection, dtype=float)
            if qc.ndim == 1 and qc.size >= 4:
                qc = qc.copy()
                for idx in BONE_IDX:
                    if idx < qc.size:
                        qc[idx] = bone_q
                ss._QCorrection = qc
                shear_str = 'ShearSoS(bone) kept' if keep_shear else 'ShearSoS(bone)=0'
                print('[run_babelbrain_case] Non-CT skull-attenuation fix applied: '
                      f'{shear_str}, QCorrArr[bone]→{bone_q:g} '
                      f'(QCorrArr now {qc.tolist()}).')
            else:
                print('[run_babelbrain_case] Warning: unexpected non-CT QCorrArr '
                      f'(shape={np.shape(ss._QCorrection)}); skull QCorr fix skipped. '
                      'Skull attenuation may be too low.', file=sys.stderr)

    patched_step1.__name__ = 'Step1_InitializeConditions'
    BabelFTD_Simulations_BASE.Step1_InitializeConditions = patched_step1
    BabelFTD_Simulations_BASE._prestus_skull_fix = True
    print('[run_babelbrain_case] Patched Step1_InitializeConditions: non-CT skull attenuation '
          'now matches CT-mode bone treatment (ShearSoS=0, QCorr for bone).')


def _patch_zdomain_reference(bb_root):
    """Neutralize BabelBrain's ZDomainStart so the transducer is placed at exactly
    the requested TxMechanicalAdjustmentZ.

    THE PROBLEM
    -----------
    CalculateRayleighFieldsForward (every geometry template)
    positions the transducer at
        Tx_z += TxMechanicalAdjustmentZ + ZDomainStart
    where ZDomainStart = CalculateDomainZReference() = GlobalStartSkin -
    LineOfSightStartSkin (SimulationConditionsBASE): the z-difference between the
    globally shallowest skin voxel anywhere in the domain and the skin voxel on
    the beam axis. In BabelBrain's own GUI workflow this compensates for oblique
    trajectories through its full-head domain. PRESTUS, however, computes the
    intended source position itself (babelbrain_placement_correction.m derives
    TxMechanicalAdjustmentZ = preset_FL - |focus_pos - trans_pos| so the apex sits
    at its true distance from the target) and hands BabelBrain a beam-aligned,
    possibly FOV-cropped mask — the skin-surface heuristic then silently moves
    the source off the requested position by GlobalStartSkin - LineOfSightStartSkin
    (e.g. -5.5 mm on a pFOV grid whose lateral corners contain skin voxels
    shallower than the on-axis skin). Every millimetre of that shift
    (a) changes the physical source-to-skull distance the FDTD actually simulates,
    and (b) desynchronizes the returned field from the tissue mask it was solved
    on once the output resample re-anchors on the tissue marker.

    THE PATCH
    ---------
    Override SimulationConditionsBASE.CalculateDomainZReference to return 0.0, so
    Tx placement honours TxMechanicalAdjustmentZ exactly. The original per-run
    value is recorded on the instance (self._prestus_zdomain_neutralized_m) and
    printed, so the magnitude of the neutralized correction stays visible in the
    log. Disable with --no_zdomain_fix to reproduce the un-patched behaviour.
    """
    add_babelbrain_paths(bb_root)
    try:
        SimulationConditionsBASE = importlib.import_module(BB_BASE).SimulationConditionsBASE
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not import SimulationConditionsBASE; '
              f'ZDomainStart neutralization NOT applied ({e}).', file=sys.stderr)
        return

    if getattr(SimulationConditionsBASE, '_prestus_zdomain_fix', False):
        return  # already patched (idempotent)

    orig_fn = SimulationConditionsBASE.CalculateDomainZReference

    def patched_zref(self):
        try:
            orig_val = float(orig_fn(self))
        except Exception as e:
            print(f'[run_babelbrain_case] ZDomainStart: original CalculateDomainZReference '
                  f'failed ({e}); using 0.')
            orig_val = float('nan')
        self._prestus_zdomain_neutralized_m = orig_val
        print(f'[run_babelbrain_case] ZDomainStart neutralized: BabelBrain wanted '
              f'{orig_val*1e3:.3f} mm domain correction on top of TxMechanicalAdjustmentZ; '
              'using 0 so the source is placed exactly as requested.')
        return 0.0

    patched_zref.__name__ = 'CalculateDomainZReference'
    SimulationConditionsBASE.CalculateDomainZReference = patched_zref
    SimulationConditionsBASE._prestus_zdomain_fix = True
    print('[run_babelbrain_case] Patched CalculateDomainZReference: ZDomainStart forced to 0 '
          '(source placed at exactly the requested TxMechanicalAdjustmentZ).')


def _patch_resave_normalized(bb_root, write_output=False):
    """Monkeypatch ResaveNormalized: fix its negative-index bug, and only
    actually write *_Sub_NORM.nii.gz when write_output is True.

    Step10_GetResults calls this solely to write *_Sub_NORM.nii.gz — a
    brain-masked, max-normalized copy of *_Sub.nii.gz used by BabelBrain's own
    GUI viewer. PRESTUS never reads it (read_babelbrain_output.m /
    read_babelbrain_thermal.m only read DataForSim.h5), so it's skipped by
    default; pass write_output=True (wired from simulation.debug) to keep it
    as a debug aid — it visualizes where the brain/skull mask boundary falls
    relative to the solved field.

    Without lower-bounds clamping, voxel coordinates that fall slightly
    outside the mask volume (affine rounding) produce negative indices, which
    numpy fancy-indexing wraps around to the end of the array rather than
    clamping to zero. The upper-bounds clamp was already in BabelBrain; this
    adds the matching lower-bounds clamp.

    The function lives in babel_integration_helpers; the base module takes it
    with `import *`, so it is rebound in the base module, where Step10 looks it
    up. bApplyOnlyMask (phase maps) masks without normalising, as upstream.
    """
    add_babelbrain_paths(bb_root)
    try:
        BASE = importlib.import_module(BB_BASE)
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not patch ResaveNormalized: {e}',
              file=sys.stderr)
        return

    if getattr(BASE, '_prestus_resave_fix', False):
        return

    import numpy as _np
    import nibabel as _nib

    def _patched_resave_normalized(rpath, mask, bApplyOnlyMask=False):
        assert '_Sub.nii.gz' in rpath
        if not write_output:
            return
        NRPath = rpath.replace('_Sub.nii.gz', '_Sub_NORM.nii.gz')
        Results = _nib.load(rpath)
        ResultsData = Results.get_fdata()
        MaskData = mask.get_fdata()
        ii, jj, kk = _np.mgrid[0:ResultsData.shape[0],
                                0:ResultsData.shape[1],
                                0:ResultsData.shape[2]]
        Indexes = _np.c_[(ii.flatten().T, jj.flatten().T, kk.flatten().T,
                          _np.ones((kk.size, 1)))].T
        PosResults = Results.affine.dot(Indexes)
        IndexesMask = _np.round(_np.linalg.inv(mask.affine).dot(PosResults)).astype(int)
        IndexesMask[0, IndexesMask[0, :] < 0] = 0
        IndexesMask[1, IndexesMask[1, :] < 0] = 0
        IndexesMask[2, IndexesMask[2, :] < 0] = 0
        IndexesMask[0, IndexesMask[0, :] >= MaskData.shape[0]] = MaskData.shape[0] - 1
        IndexesMask[1, IndexesMask[1, :] >= MaskData.shape[1]] = MaskData.shape[1] - 1
        IndexesMask[2, IndexesMask[2, :] >= MaskData.shape[2]] = MaskData.shape[2] - 1
        SubMask = MaskData[IndexesMask[0, :],
                           IndexesMask[1, :],
                           IndexesMask[2, :]].reshape(ResultsData.shape)
        ResultsData[SubMask < 4] = 0
        if not bApplyOnlyMask:
            ResultsData /= ResultsData.max()
        NormalizedNifti = _nib.Nifti1Image(ResultsData, Results.affine, header=Results.header)
        NormalizedNifti.to_filename(NRPath)

    _patched_resave_normalized.__name__ = 'ResaveNormalized'
    BASE.ResaveNormalized = _patched_resave_normalized
    BASE._prestus_resave_fix = True
    state = 'writing *_Sub_NORM.nii.gz (debug mode)' if write_output else 'skipped (unused *_Sub_NORM.nii.gz GUI output)'
    print(f'[run_babelbrain_case] Patched ResaveNormalized: negative-index bounds clamp added; {state}.')


def _patch_calculate_temperature_effects(bb_root):
    """Monkeypatch CalculateTemperatureEffects to guard against empty tissue masks.

    BabelBrain calls .max() on boolean-indexed slices of ResTemp and FinalDose
    for brain, skin, and skull regions.  If any tissue class is absent from the
    simulation domain the mask is all-False and .max() raises ValueError on an
    empty array.

    This patch rewrites the six bare .max() calls (TI / TIS / TIC and CEMBrain /
    CEMSkin / CEMSkull) to guard with ``if <mask>.any() else 0.0``.

    Implementation: uses inspect.getsource + string replacement + exec to redefine
    just the function in the module namespace — the 1185-line function body is NOT
    copied into PRESTUS.  The patch is a safe no-op (with a warning) if the source
    does not match (different BabelBrain version) or is already guarded.
    """
    import inspect
    import textwrap

    add_babelbrain_paths(bb_root)
    try:
        import ThermalModeling.CalculateTemperatureEffects as CTE_mod
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not patch CalculateTemperatureEffects: {e}',
              file=sys.stderr)
        return

    if getattr(CTE_mod, '_prestus_thermal_fix', False):
        return

    orig_fn = CTE_mod.CalculateTemperatureEffects
    try:
        src = textwrap.dedent(inspect.getsource(orig_fn))
    except (OSError, TypeError) as e:
        print(f'[run_babelbrain_case] Warning: could not get source of '
              f'CalculateTemperatureEffects: {e}', file=sys.stderr)
        return

    if 'SelBrain.any()' in src:
        # Already guarded (either by external edit or a prior run of this patch).
        CTE_mod._prestus_thermal_fix = True
        return

    _REPLACEMENTS = [
        ('ResTemp[SelBrain].max()',
         'ResTemp[SelBrain].max() if SelBrain.any() else 0.0'),
        ('ResTemp[SelSkin].max()',
         'ResTemp[SelSkin].max() if SelSkin.any() else 0.0'),
        ('ResTemp[SelSkull].max()',
         'ResTemp[SelSkull].max() if SelSkull.any() else 0.0'),
        ('FinalDose[SelBrain].max()/60',
         'FinalDose[SelBrain].max()/60 if SelBrain.any() else 0.0'),
        ('FinalDose[SelSkin].max()/60',
         'FinalDose[SelSkin].max()/60 if SelSkin.any() else 0.0'),
        ('FinalDose[SelSkull].max()/60',
         'FinalDose[SelSkull].max()/60 if SelSkull.any() else 0.0'),
    ]
    if not all(old in src for old, _ in _REPLACEMENTS):
        print('[run_babelbrain_case] Warning: CalculateTemperatureEffects source does not '
              'match expected pattern; empty-mask guard NOT applied (BabelBrain version '
              'mismatch?).', file=sys.stderr)
        return

    src_patched = src
    for old, new in _REPLACEMENTS:
        src_patched = src_patched.replace(old, new, 1)

    try:
        code = compile(src_patched, inspect.getfile(orig_fn), 'exec')
        exec(code, CTE_mod.__dict__)  # noqa: S102 — intentional monkeypatch
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: failed to apply '
              f'CalculateTemperatureEffects patch: {e}', file=sys.stderr)
        return

    CTE_mod._prestus_thermal_fix = True
    print('[run_babelbrain_case] Patched CalculateTemperatureEffects: '
          'empty-mask guards added (TI/TIS/TIC/CEM*).')


def _patch_annular_surface(bb_root):
    """Monkeypatch generate_annular_curved_surface: deterministic ring-band tessellation.

    tx_geometries.generate_annular_curved_surface lays out each ring's surface
    elements with two `np.arange(step/2, step*(0.5+n), step)` calls, one over the
    polar-angle bands (betas) and one azimuthal (alphas). Because the stop is a
    floating-point multiple of the step, np.arange intermittently emits n+1
    samples instead of n, appending a band beyond the ring's outer radius.

    Observed on the 10-ring R15473 at 300 kHz: ring 9 got 4 bands instead of 3,
    extending its outer radius from 30.75 mm to 31.215 mm (into ring 10, erasing
    the kerf). That inflated ring 9's area ~41% and the free-water focal
    intensity ~13% above k-Wave / jwave / analytical Rayleigh.

    Fix: count-based sampling, `start + step*(np.arange(int(n)) + 0.5)`, which
    gives exactly n samples with the same values. The function is shared by the
    focused annular, single-element and flat annular templates, which resolve it
    from tx_geometries' globals at call time, so rebinding it there covers all of
    them. Implementation: inspect.getsource + targeted replacement + exec into
    the module namespace. Safe no-op (with a warning) if the source does not match.
    """
    import inspect
    import textwrap

    add_babelbrain_paths(bb_root)
    try:
        TXG = importlib.import_module(BB_TXGEOM)
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not patch ring tessellation: {e}',
              file=sys.stderr)
        return

    if getattr(TXG, '_prestus_surface_fix', False):
        return

    orig_fn = getattr(TXG, 'generate_annular_curved_surface', None)
    if orig_fn is None:
        print('[run_babelbrain_case] Warning: generate_annular_curved_surface not found in '
              'tx_geometries; ring-tessellation fix NOT applied.', file=sys.stderr)
        return
    try:
        src = textwrap.dedent(inspect.getsource(orig_fn))
    except (OSError, TypeError) as e:
        print(f'[run_babelbrain_case] Warning: could not get source of '
              f'generate_annular_curved_surface: {e}', file=sys.stderr)
        return

    _REPLACEMENTS = [
        ('np.arange(inner_beta + beta_step / 2, inner_beta + beta_step * (1 / 2 + num_steps), beta_step)',
         '(inner_beta + beta_step * (np.arange(int(num_steps)) + 0.5))'),
        ('np.arange(alpha_step / 2, alpha_step * (1 / 2 + num_alphas), alpha_step)',
         '(alpha_step * (np.arange(int(num_alphas)) + 0.5))'),
    ]
    if not all(old in src for old, _ in _REPLACEMENTS):
        print('[run_babelbrain_case] Warning: generate_annular_curved_surface source does not '
              'match expected pattern; ring-tessellation fix NOT applied (BabelBrain '
              'version mismatch?).', file=sys.stderr)
        return

    src_patched = src
    for old, new in _REPLACEMENTS:
        src_patched = src_patched.replace(old, new, 1)

    try:
        code = compile(src_patched, inspect.getfile(orig_fn), 'exec')
        exec(code, TXG.__dict__)  # noqa: S102 — intentional monkeypatch
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: failed to apply ring-tessellation '
              f'patch: {e}', file=sys.stderr)
        return

    TXG._prestus_surface_fix = True
    print('[run_babelbrain_case] Patched generate_annular_curved_surface: deterministic '
          'ring-band tessellation (fixes np.arange endpoint overshoot inflating ring area).')


# dB/cm → Np/m : 1 Np = 8.685889638 dB, and 100 cm per m.
_DB_CM_TO_NP_M = 100.0 / 8.685889638

# PRESTUS medium_properties field → BabelBrain MatFreq material key. The generic
# 'skull' is listed first so that, when a config also defines the split
# cortical/trabecular tissues, those specific entries overwrite it.
# 'Water' here only updates MatFreq[f]['Water'], which BabelBrain reads in
# GetSmallestSOS (grid spacing). The simulated water material comes from the
# module-level Material (1120 kHz loop leftover) and stays at 1000 kg/m³ /
# 1500 m/s. Deliberately not patched: water also sets FactorConvPtoU, so
# changing it would shift the source amplitude (see doc_babelbrain.md).
_PRESTUS_TO_BB_MATERIAL = [
    ('water',            'Water'),
    ('skin',             'Skin'),
    ('brain',            'Brain'),
    ('skull',            'Cortical'),     # undifferentiated skull → both bone slots …
    ('skull',            'Trabecular'),
    ('skull_cortical',   'Cortical'),     # … overridden by split tissues when present
    ('skull_trabecular', 'Trabecular'),
]


def _write_air_regions(basedir, subject_id, target, freq, ppw, hu_window=(-1200.0, -400.0),
                       keep_boundary=False):
    """Write the AirRegions.nii.gz BabelBrain reads as its reflector mask (AIRMASK).

    Mirrors BabelDatasetPreps.py's air extraction (bExtractAirRegions): voxels
    of the raw-HU CT inside hu_window, 3-voxel median filter, connected
    components, the largest component (air around the head) removed. In
    addition, components touching the domain boundary are removed unless
    keep_boundary is set: PRESTUS hands BabelBrain an already-cropped grid, in
    which the exterior air -- including the coupling gap between transducer and
    scalp, which is air (-1000 HU) in the pseudoCT but water in the
    simulation -- can be split into several pieces.

    Must run before _quantize_ct_for_babelbrain, which overwrites CT.nii.gz
    with the index map. BabelBrain (babel_integration_base.RunCases) loads the
    file from basedir/ID/<target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz with the
    same axis-2 flip as the CT, so writing it on the CT grid keeps them aligned.
    In the solver, stress and particle velocity are never updated inside the
    mask (ReflectorMask), i.e. the region acts as a pressure-release reflector.
    """
    import nibabel as nib
    from scipy import ndimage

    freq_kHz = round(freq / 1e3)
    stem = os.path.join(basedir, subject_id, f'{target}_{freq_kHz}kHz_{ppw}PPW_')
    ct_path = stem + 'CT.nii.gz'
    if not os.path.isfile(ct_path):
        raise FileNotFoundError(f'[run_babelbrain_case] air regions: {ct_path} not found.')
    ct_nii = nib.load(ct_path)
    ct = np.asarray(ct_nii.get_fdata(), dtype=np.float64)

    lo, hi = float(hu_window[0]), float(hu_window[1])
    air = ((ct >= lo) & (ct <= hi)).astype(np.uint8)     # NaN padding compares False
    air = ndimage.median_filter(air, 3)
    labels, n = ndimage.label(air == 1, structure=np.ones((3, 3, 3)))
    n_raw = int(air.sum())
    removed = 0
    if n > 0:
        sizes = np.bincount(labels.ravel())[1:]
        drop = {int(np.argmax(sizes)) + 1}
        if not keep_boundary:
            faces = [labels[0], labels[-1], labels[:, 0], labels[:, -1], labels[:, :, 0], labels[:, :, -1]]
            for f in faces:
                drop.update(int(v) for v in np.unique(f) if v > 0)
        drop_mask = np.isin(labels, list(drop))
        removed = int(drop_mask.sum())
        air[drop_mask] = 0
    n_kept = int(air.sum())
    vox_mm3 = float(np.prod(ct_nii.header.get_zooms()[:3]))
    print(f'[run_babelbrain_case] Air regions: HU [{lo:.0f}, {hi:.0f}], {n} component(s); '
          f'{removed} exterior/boundary voxels removed, {n_kept} internal air voxels kept '
          f'({n_kept * vox_mm3 / 1e3:.2f} cm³).')
    if n_kept == 0:
        print('[run_babelbrain_case] Air regions: no internal air in the simulation domain '
              '(the mask is empty; results equal a run without --air_regions).')

    out = stem + 'AirRegions.nii.gz'
    img = nib.Nifti1Image(air, ct_nii.affine)
    img.header.set_data_dtype(np.uint8)
    img.to_filename(out)
    print(f'[run_babelbrain_case] Air regions written: {out}')
    return out, n_kept


def _quantize_ct_for_babelbrain(basedir, subject_id, target, freq, ppw):
    """Convert the raw-HU CT.nii.gz written by write_babelbrain_input.m into the
    quantized index-map format + companion CT-cal.npz that BabelBrain's CT path
    actually requires.

    write_babelbrain_input.m writes CT.nii.gz as a raw Hounsfield-equivalent
    volume, resampled/gridded identically to BabelViscoInput.nii.gz -- but
    BabelBrain (BabelDatasetPreps.py Step1) never reads a raw-HU CT.nii.gz
    directly. It expects CTFNAME to already be a quantized index map (small
    integers indexing into a sorted table of unique HU levels), with that
    table saved separately as `<CT.nii.gz prefix>CT-cal.npz` (key
    'UniqueHU') -- babel_integration_base.py loads both:
        DensityCTMap = flip(nibabel.load(CTFNAME).get_fdata(), axis=2).astype(uint32)
        AllBoneHU = np.load(CTFNAME.split('CT.nii.gz')[0]+'CT-cal.npz')['UniqueHU']
    and uses DensityCTMap values directly as indices into AllBoneHU / the
    material table, so an unquantized float/HU CT.nii.gz is either a crash
    (index far out of range) or silently wrong material properties throughout
    the skull.

    Quantization mirrors BabelDatasetPreps.py's own scheme exactly (uniform
    round-to-M-levels quantization, M=2**CT_quantification-1, default
    CT_quantification=10 bits), restricted to the bone region -- taken here
    from the label-2/3 (cortical/trabecular) voxels of the tissue mask
    write_babelbrain_input.m writes alongside it, since PRESTUS has no
    separate bone-segmentation mask of its own to reuse. Voxels outside the
    bone region map to index 0, matching BabelDatasetPreps.py's
    `ndataCTMap[nfct==False]=0`.
    """
    import nibabel as nib

    freq_kHz = round(freq / 1e3)
    ct_path = os.path.join(basedir, subject_id, f'{target}_{freq_kHz}kHz_{ppw}PPW_CT.nii.gz')
    mask_path = os.path.join(basedir, subject_id, f'{target}_{freq_kHz}kHz_{ppw}PPW_BabelViscoInput.nii.gz')

    if not os.path.isfile(ct_path):
        raise FileNotFoundError(
            f'[run_babelbrain_case] CT quantization: {ct_path} not found. '
            'write_babelbrain_input.m should have written it alongside BabelViscoInput.nii.gz.')
    if not os.path.isfile(mask_path):
        raise FileNotFoundError(f'[run_babelbrain_case] CT quantization: {mask_path} not found.')

    ct_nii = nib.load(ct_path)
    ct_data = np.asarray(ct_nii.get_fdata(), dtype=np.float64)
    mask_data = np.asarray(nib.load(mask_path).get_fdata())

    bone = (mask_data == 2) | (mask_data == 3)
    if not np.any(bone):
        raise ValueError(
            '[run_babelbrain_case] CT quantization: no cortical/trabecular (label 2/3) voxels found '
            f'in {mask_path} -- cannot determine the bone region to quantize CT values over.')

    CT_quantification = 10   # bits; matches BabelDatasetPreps.py's own default
    bone_vals = ct_data[bone]
    min_data = bone_vals.min()
    max_data = bone_vals.max()
    a_range = max_data - min_data
    if a_range <= 0:
        raise ValueError(
            f'[run_babelbrain_case] CT quantization: bone-region HU range is degenerate '
            f'(min={min_data}, max={max_data}) in {ct_path}.')
    m_levels = 2 ** CT_quantification - 1
    res_step = a_range / m_levels
    qx = res_step * np.round((m_levels / a_range) * (bone_vals - min_data)) + min_data
    ct_data[bone] = qx

    unique_hu = np.unique(ct_data[bone])
    print(f'[run_babelbrain_case] CT quantization: {len(unique_hu)} unique bone HU levels '
          f'(range [{min_data:.1f}, {max_data:.1f}]).')

    ct_map = np.zeros(ct_data.shape, dtype=np.uint32)
    for n, d in enumerate(unique_hu):
        ct_map[ct_data == d] = n
    ct_map[~bone] = 0

    cal_path = ct_path[:ct_path.rindex('CT.nii.gz')] + 'CT-cal.npz'
    np.savez_compressed(cal_path, UniqueHU=unique_hu)
    print(f'[run_babelbrain_case] CT calibration table written: {cal_path}')

    ct_map_img = nib.Nifti1Image(ct_map, ct_nii.affine)
    ct_map_img.header.set_data_dtype(np.uint32)
    ct_map_img.to_filename(ct_path)
    print(f'[run_babelbrain_case] CT index map written (overwriting raw-HU CT.nii.gz): {ct_path}')


def _inject_prestus_materials(bb_root, freq, props_json_path, bone_shear=False):
    """Overwrite BabelBrain's built-in non-CT material database (``MatFreq``) with
    PRESTUS's configured medium properties, so k-Wave / jwave / BabelBrain all
    assume identical tissue acoustics.

    PRESTUS stores attenuation as ``alpha_coeff`` in dB/(cm·MHz^alpha_power);
    BabelBrain's ``MatFreq`` stores a per-frequency longitudinal attenuation in
    Np/m. This converts ``alpha_coeff * f_MHz**alpha_power`` (dB/cm) to Np/m.

    Shear speed / attenuation are set to 0 (fluid model), consistent with the CT
    path and with ``_patch_noncat_skull_attenuation``. The realised bone
    attenuation still depends on the Q-correction applied there; a validation run
    against k-Wave is the check that the two agree.

    With ``bone_shear`` (babelbrain.shear: prestus_elastic), a skull layer with a
    ``babelbrain_elastic`` entry gets its shear speed [m/s], longitudinal and shear
    attenuation [Np/m at the drive frequency] from there. PRESTUS
    (babelbrain_skull_elastic.m) has already converted them: kwave_elastic's y = 2
    coefficient a2 [dB/(cm MHz^2)], exact at f0, times f_MHz^2 gives dB/cm, times
    100 / 8.686 gives Np/m, the unit of the MatFreq / AddMaterial columns
    [Density, LongSoS, ShearSoS, LongAtt, ShearAtt].

    Returns the set of BabelBrain material keys that received elastic values.
    """
    import json

    add_babelbrain_paths(bb_root)
    try:
        bbase = importlib.import_module(BB_HELPERS)
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not import babel_integration_helpers; '
              f'PRESTUS material injection skipped ({e}).', file=sys.stderr)
        return set()

    try:
        with open(props_json_path) as f:
            props = json.load(f)
    except Exception as e:
        print(f'[run_babelbrain_case] Warning: could not read medium properties JSON '
              f'{props_json_path}; material injection skipped ({e}).', file=sys.stderr)
        return set()

    MatFreq = getattr(bbase, 'MatFreq', None)
    if not MatFreq:
        print('[run_babelbrain_case] Warning: MatFreq not found in babel_integration_helpers; '
              'PRESTUS material injection skipped.', file=sys.stderr)
        return set()

    # Locate the frequency key BabelBrain will actually index (MatFreq[self._Frequency]).
    keys = np.array(sorted(MatFreq.keys()), dtype=float)
    k = float(keys[np.argmin(np.abs(keys - freq))])
    if abs(k - freq) > 1.0:
        print(f'[run_babelbrain_case] Warning: no MatFreq entry within 1 Hz of {freq:.0f} Hz '
              f'(closest {k:.0f}); PRESTUS material injection skipped.', file=sys.stderr)
        return set()

    f_MHz = freq / 1e6
    applied = {}
    elastic = set()
    for pk, bk in _PRESTUS_TO_BB_MATERIAL:
        m = props.get(pk)
        if not isinstance(m, dict):
            continue
        try:
            density = float(m['density'])
            c_long  = float(m['sound_speed'])
        except (KeyError, TypeError, ValueError):
            continue
        a_coeff = float(m.get('alpha_coeff', 0.0) or 0.0)
        a_power = float(m.get('alpha_power', 1.0) or 1.0)
        long_att_npm = a_coeff * (f_MHz ** a_power) * _DB_CM_TO_NP_M
        shear_sos, shear_att_npm = 0.0, 0.0
        el = m.get('babelbrain_elastic') if bone_shear else None
        if isinstance(el, dict):
            shear_sos     = float(el['shear_sound_speed'])
            shear_att_npm = float(el['shear_att_np_m'])
            long_att_npm  = float(el['long_att_np_m'])     # already Np/m at f0
            elastic.add(bk)
        else:
            elastic.discard(bk)    # a later fluid entry overrides an earlier elastic one
        MatFreq[k][bk] = np.array([density, c_long, shear_sos, long_att_npm, shear_att_npm])
        applied[bk] = (density, c_long, long_att_npm, shear_sos, shear_att_npm)

    if applied:
        print(f'[run_babelbrain_case] Injected PRESTUS materials into MatFreq[{k:.0f} Hz]'
              + (':' if bone_shear else ' (shear set to 0):'))
        for bk, (d, c, a, cs, sa) in applied.items():
            print(f'    {bk:11s} density={d:7.1f} kg/m³  LongSoS={c:7.1f} m/s  '
                  f'LongAtt={a:8.2f} Np/m  ShearSoS={cs:7.1f} m/s  ShearAtt={sa:8.2f} Np/m')
    else:
        print('[run_babelbrain_case] Warning: no PRESTUS materials matched known tissues; '
              'BabelBrain built-in database left unchanged.', file=sys.stderr)
    return elastic


def main():
    args = parse_args()

    # Set up log file tee (stdout + stderr → console and file simultaneously)
    _log_fh = None
    if args.log_file:
        os.makedirs(os.path.dirname(os.path.abspath(args.log_file)), exist_ok=True)
        _log_fh = open(args.log_file, 'w', buffering=1)
        sys.stdout = _TeeStream(sys.__stdout__, _log_fh)
        sys.stderr = _TeeStream(sys.__stderr__, _log_fh)

    args.device_name = resolve_device_name(args.computing_backend, args.device_name)

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

    add_babelbrain_paths(bb_root)
    try:
        preset = resolve_preset(bb_root, args.tx_system)
        if preset['name'] != args.tx_system:
            print(f"[run_babelbrain_case] Preset '{args.tx_system}' is now '{preset['name']}' in BabelBrain.")
        RunSimClass = load_run_sim(preset)
    except (ImportError, OSError) as e:
        print(f'[run_babelbrain_case] Import error: {e}', file=sys.stderr)
        sys.exit(1)

    # Fix the non-CT path so the skull actually attenuates (see function docstring).
    # No-op for CT / water-only runs; disabled with --no_skull_atten_fix.
    if args.no_skull_atten_fix:
        print('[run_babelbrain_case] --no_skull_atten_fix set: leaving BabelBrain non-CT '
              'skull attenuation UN-patched (expect inflated ISPPA).')
    else:
        _patch_noncat_skull_attenuation(
            bb_root, bone_q=args.q_correction if args.q_correction > 0 else 1.0,
            keep_shear=args.bone_shear)

    # Place the source at exactly the requested TxMechanicalAdjustmentZ (see
    # function docstring): BabelBrain otherwise adds a skin-surface-dependent
    # ZDomainStart on top, silently moving the source millimetres off the
    # PRESTUS-intended position and desynchronizing field vs tissue after the
    # marker-anchored output resample. Disabled with --no_zdomain_fix.
    if args.no_zdomain_fix:
        print('[run_babelbrain_case] --no_zdomain_fix set: leaving BabelBrain ZDomainStart '
              'UN-patched (source placement will include the skin-surface domain correction; '
              'tissue-anchored output alignment will be off by that amount).')
    else:
        _patch_zdomain_reference(bb_root)

    # Replace BabelBrain's built-in material DB with PRESTUS's configured values so
    # all solver backends assume identical tissue acoustics. Requires the skull
    # fix above to be active (otherwise the injected bone attenuation is not
    # realised). Water-only runs do not read it (BabelBrain's own water). On the
    # CT path bone comes from the HU mapping, but Skin and Brain still come from
    # MatFreq (babel_integration_helpers.py), so their injected values apply there.
    bone_shear = args.bone_shear and not args.water_only
    if bone_shear and (args.use_ct or not args.medium_properties_json):
        print('[run_babelbrain_case] --bone_shear needs --medium_properties_json and the non-CT '
              'path (BabelBrain sets bone shear to 0 on the CT path).', file=sys.stderr)
        sys.exit(1)
    if args.medium_properties_json:
        if args.no_skull_atten_fix:
            print('[run_babelbrain_case] Warning: --medium_properties_json given together with '
                  '--no_skull_atten_fix; injected bone attenuation will NOT be realised (QCorr '
                  'left at 1.0 for bone).', file=sys.stderr)
        elastic = _inject_prestus_materials(bb_root, args.freq, args.medium_properties_json,
                                            bone_shear=bone_shear)
        if bone_shear and not {'Cortical', 'Trabecular'} <= elastic:
            print('[run_babelbrain_case] --bone_shear: elastic properties were not injected for '
                  f'both bone materials (got {sorted(elastic)}); stopping instead of running '
                  'with BabelBrain\'s built-in shear values.', file=sys.stderr)
            sys.exit(1)

    # Defensive crash-prevention patches (no result change for normal inputs).
    _patch_resave_normalized(bb_root, write_output=args.debug)
    _patch_calculate_temperature_effects(bb_root)

    # Correctness patch: deterministic annular ring-band tessellation. Fixes a
    # np.arange endpoint overshoot in generate_annular_curved_surface that intermittently adds a
    # spurious outer sub-ring (inflating that ring's area / the focal intensity).
    _patch_annular_surface(bb_root)

    # Load preset geometry from default.yaml
    preset_cfg = load_preset_config(preset)

    # Build geometry kwargs for RunCases (passed through to BabelFTD_Simulations).
    # The YAML → kwarg names follow BabelBrain's own Step 2 caller; phased arrays
    # take their element table from here (the device classes hold no geometry).
    geom_kw = {}
    for yaml_key, kw_key, conv in (
            ('TxDiam',                     'Aperture',                       float),
            ('FocalLength',                'FocalLength',                    float),
            ('InDiameters',                'InDiameters',                    lambda v: np.array([float(x) for x in v])),
            ('OutDiameters',               'OutDiameters',                   lambda v: np.array([float(x) for x in v])),
            ('Elements',                   'elements',                       lambda v: v),
            ('NumElements',                'num_elements',                   int),
            ('ElementSize',                'element_size',                   float),
            ('CoordinateSystem',           'coordinate_system',              str),
            ('DistanceTxBottomToOutplane', 'distance_tx_bottom_to_outplane', float)):
        if yaml_key in preset_cfg:
            geom_kw[kw_key] = conv(preset_cfg[yaml_key])

    # Geometry passed by PRESTUS (the PRESTUS-geometry preset and single-element
    # presets have none in default.yaml).
    for key, val in (('Aperture', args.aperture), ('FocalLength', args.focal_length)):
        if val > 0:
            if key in geom_kw and abs(geom_kw[key] - val) > 1e-6:
                print(f'[run_babelbrain_case] {key}: command line {val} m overrides preset {geom_kw[key]} m')
            geom_kw[key] = val
    for key, val in (('InDiameters', args.in_diameters), ('OutDiameters', args.out_diameters)):
        if val:
            geom_kw[key] = np.array([float(x) for x in val.split(',')])

    # Fail here with the fix, not a KeyError/ValueError deep in BabelBrain: the
    # simple_focused RUN_SIM.RunCases reads Aperture/FocalLength with no default,
    # and the PRESTUS-geometry preset has no other source for its rings.
    needed = {'simple_focused': ('Aperture', 'FocalLength'),
              PRESTUS_GEOMETRY: ('Aperture', 'FocalLength', 'InDiameters', 'OutDiameters')}
    key = PRESTUS_GEOMETRY if preset['name'] == PRESTUS_GEOMETRY else preset['type']
    if key in needed:
        missing = [k for k in needed[key] if k not in geom_kw]
        if missing:
            raise ValueError(f'{args.tx_system} needs {", ".join(missing)}: pass --aperture, '
                             f'--focal_length [m] (and --in_diameters/--out_diameters for PRESTUS).')
        if key == PRESTUS_GEOMETRY and len(geom_kw['InDiameters']) != len(geom_kw['OutDiameters']):
            raise ValueError('--in_diameters and --out_diameters need one value per ring.')

    # Only forward geometry the preset's sim class accepts by name: each geometry
    # family takes a different subset (e.g. flat annular H246 has no FocalLength),
    # and BabelFTD_Simulations_BASE rejects the rest.
    sim_cls = getattr(RunSimClass, '_BabelFTDSimClass', None)
    if geom_kw and sim_cls is not None:
        import inspect
        accepted = set()
        for cls in sim_cls.__mro__:
            if '__init__' in cls.__dict__:
                accepted |= {n for n, p in inspect.signature(cls.__init__).parameters.items()
                             if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)}
        dropped = sorted(k for k in geom_kw if k not in accepted)
        if dropped:
            print(f'[run_babelbrain_case] {args.tx_system}: preset class takes no '
                  f'{", ".join(dropped)}; using its built-in geometry for these.')
            geom_kw = {k: v for k, v in geom_kw.items() if k in accepted}

    if geom_kw:
        print(f'[run_babelbrain_case] Preset geometry: '
              f"Aperture={geom_kw.get('Aperture','?')} m  "
              f"FocalLength={geom_kw.get('FocalLength','?')} m  "
              f"n_rings={len(geom_kw.get('InDiameters',[]))}")

    R = RunSimClass()

    # Build RunCases kwargs
    kw = dict(
        TxSystem       = preset['name'],   # stored in DataForSim; read by the thermal step
        targets        = [args.target],
        ID             = args.subject_id,
        basedir        = args.basedir,
        Frequencies    = [args.freq],
        basePPW        = [args.ppw],
        COMPUTING_BACKEND           = args.computing_backend,
        deviceName                  = args.device_name,
        TxMechanicalAdjustmentZ     = args.tx_mech_adjust_z,
        ZSteering                   = args.z_steering,
        bUseRayleighForWater        = (not args.no_rayleigh_for_water),
        bWaterOnly                  = args.water_only,
        bTightNarrowBeamDomain      = not args.no_tight_narrow_beam,
        zLengthBeyonFocalPointWhenNarrow = args.z_beyond_focus_mm / 1e3,  # mm → m
        bDoRefocusing               = args.do_refocusing,
        bDisplay                    = False,
        bDryRun                     = args.dry_run,
        **geom_kw,
    )

    if args.use_ct:
        kw['bUseCT'] = True
        kw['MappingMethod'] = args.ct_mapping
        kw['bPETRA'] = args.use_petra
        # BabelBrain only reads an air mask when bExtractAirRegions=True
        # (babel_integration_base.py RunCases sets AIRMASK to
        # <target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz, else None). With
        # --air_regions PRESTUS writes that file from the raw-HU CT first
        # (before quantization overwrites CT.nii.gz). Never for water-only runs,
        # where the reflector would sit in free water.
        use_air = args.air_regions and not args.water_only
        if use_air:
            _write_air_regions(args.basedir, args.subject_id, args.target, args.freq, args.ppw,
                               args.air_hu_window, args.air_keep_boundary)
        elif args.air_regions:
            print('[run_babelbrain_case] --air_regions ignored for the water-only run.')
        kw['bExtractAirRegions'] = use_air
        print(f'[run_babelbrain_case] CT path enabled: mapping={args.ct_mapping}'
              f'{"  bPETRA=True" if args.use_petra else ""}'
              f'{"  air regions=on" if use_air else ""}')
        _quantize_ct_for_babelbrain(args.basedir, args.subject_id, args.target,
                                     args.freq, args.ppw)
    elif args.air_regions:
        print('[run_babelbrain_case] Warning: --air_regions needs --use_ct (pseudoCT path); ignored.',
              file=sys.stderr)

    if args.weights_file and os.path.isfile(args.weights_file):
        kw['OptimizedWeightsFile'] = convert_weights_if_needed(args.weights_file)
        print(f'[run_babelbrain_case] Using weights: {kw["OptimizedWeightsFile"]}')
        src_dump = None
        if args.debug:
            src_dump = os.path.join(args.basedir, args.subject_id,
                                    f'{args.target}_babelbrain_source_dump.h5')
            print(f'[run_babelbrain_case] Debug: will dump BabelBrain source geometry to {src_dump}')
        # Written on every run with a weights file: source-placement diagnostic
        # (realized vs requested Tx position). read_babelbrain_output.m reads it back and warns when the
        # residual exceeds a native voxel (see _write_anchor_sidecar).
        anchor_sidecar = os.path.join(args.basedir, args.subject_id,
                                      f'{args.target}_prestus_anchor.json')
        _patch_annular_array_weights(RunSimClass, dump_path=src_dump, anchor_path=anchor_sidecar,
                                     tx_mech_adjust_z=args.tx_mech_adjust_z)
    else:
        print('[run_babelbrain_case] No weights file — using geometric focusing.')

    if args.q_correction > 0:
        # BabelFTD_Simulations_BASE.__init__ reads this (and only this exact key)
        # to override its QCorrection default (CT path); the non-CT skull patch
        # receives the same value directly.
        os.environ['BABEL_PYTEST_QFACTOR'] = repr(float(args.q_correction))
        print(f'[run_babelbrain_case] Bone QCorrection override: {args.q_correction:g}')

    if args.source_amp_pa > 0:
        kw['pressure'] = args.source_amp_pa
        print(f'[run_babelbrain_case] Source pressure amplitude: {args.source_amp_pa:.1f} Pa '
              f'({args.source_amp_pa/1e3:.2f} kPa)')

    if args.dry_run:
        print('[run_babelbrain_case] DRY RUN — no GPU computation will occur.')

    print(f'[run_babelbrain_case] Starting RunCases for '
          f'{args.subject_id}/{args.target} '
          f'({args.tx_system}, {args.freq/1e3:.0f} kHz, PPW={args.ppw})')
    print(f'[run_babelbrain_case] Water field solver: '
          f'{"Rayleigh-Sommerfeld (analytic)" if not args.no_rayleigh_for_water else "full-wave FDTD"} '
          f'(bUseRayleighForWater={not args.no_rayleigh_for_water})')

    _start_memory_monitor()
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
    elif args.water_only:
        # Some BabelBrain versions return None/empty for water-only runs even on success.
        print('[run_babelbrain_case] Water-only run completed (output file list may be empty).')
    else:
        print('[run_babelbrain_case] Warning: RunCases returned no output files.',
              file=sys.stderr)
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
