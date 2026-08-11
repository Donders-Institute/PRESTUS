#!/usr/bin/env python3
"""
run_babelbrain_thermal.py — PRESTUS → BabelBrain thermal simulation bridge

Called by babelbrain_thermal_wrapper.m via run_command() (java.lang.ProcessBuilder,
no shell). Invokes BabelBrain's
CalculateThermalProcess on a *DataForSim.h5 acoustic output file and writes
a *ThermalField*.h5 result that PRESTUS reads back via read_babelbrain_thermal.m.

Usage:
    python run_babelbrain_thermal.py \
        --acoustic_h5 /path/to/*DataForSim.h5 \
        --isppa 5.0 \
        --duty_cycle 0.3 \
        --prf 1000 \
        --duration_s 80 \
        --computing_backend 3 \
        --device_name Apple \
        [--babelbrain_root /path/to/BabelBrain] \
        [--baseline_temp 37.0] \
        [--sel_p p_amp]

Exit codes:
    0  — success
    1  — error (message printed to stderr)
"""

import argparse
import sys
import os


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


def parse_args():
    p = argparse.ArgumentParser(description='PRESTUS–BabelBrain thermal bridge')
    p.add_argument('--acoustic_h5',        required=True,
                   help='Path to *DataForSim.h5 from acoustic simulation')
    p.add_argument('--isppa',              required=True, type=float,
                   help='ISPPA [W/cm²] at the focus')
    p.add_argument('--duty_cycle',         required=True, type=float,
                   help='Duty cycle [0–1]')
    p.add_argument('--prf',                required=True, type=float,
                   help='Pulse repetition frequency [Hz]')
    p.add_argument('--duration_s',         required=True, type=float,
                   help='Pulse train duration [s]')
    p.add_argument('--computing_backend',  required=True, type=int,
                   help='1=CUDA 2=OpenCL 3=Metal 4=CPU')
    p.add_argument('--device_name',        required=True)
    p.add_argument('--babelbrain_root',    default='',
                   help='BabelBrain source root; auto-detected if empty')
    p.add_argument('--baseline_temp',      default=37.0, type=float,
                   help='Baseline temperature [°C] (default 37.0)')
    p.add_argument('--sel_p',             default='p_amp',
                   help='Pressure key in DataForSim.h5 (default: p_amp; '
                        'use p_amp_refocus if refocusing was used)')
    p.add_argument('--duration_off_s',     default=0.0, type=float,
                   help='Off duration per cycle [s] (default 0)')
    p.add_argument('--repetitions',        default=1, type=int,
                   help='Number of sonication repetitions (default 1)')
    p.add_argument('--frequency',          default=7e5, type=float,
                   help='Transducer centre frequency [Hz] (default 700000)')
    p.add_argument('--tx_system',          default='',
                   help='BabelBrain TxSystem identifier (default empty string)')
    p.add_argument('--limit_bhte_iter',    default=100, type=int,
                   help='Max BHTE iterations per subprocess (default 100)')
    p.add_argument('--log_file',            default='',
                   help='Path to write a copy of all stdout/stderr output for post-hoc debugging.')
    return p.parse_args()


def find_babelbrain_root(hint=''):
    candidates = []
    if hint:
        candidates.append(hint)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(script_dir, '..', '..', 'external', 'BabelBrain'))
    for c in candidates:
        c = os.path.normpath(c)
        if os.path.isdir(os.path.join(c, 'TranscranialModeling')):
            return c
    raise FileNotFoundError(
        "Cannot locate BabelBrain source root.\n"
        "Ensure the submodule is initialised: git submodule update --init external/BabelBrain\n"
        "Or pass --babelbrain_root explicitly."
    )


def main():
    args = parse_args()

    # Set up log file tee (stdout + stderr → console and file simultaneously)
    _log_fh = None
    if args.log_file:
        os.makedirs(os.path.dirname(os.path.abspath(args.log_file)), exist_ok=True)
        _log_fh = open(args.log_file, 'w', buffering=1)
        sys.stdout = _TeeStream(sys.__stdout__, _log_fh)
        sys.stderr = _TeeStream(sys.__stderr__, _log_fh)

    if not os.path.isfile(args.acoustic_h5):
        print(f'[run_babelbrain_thermal] Acoustic HDF5 not found: {args.acoustic_h5}',
              file=sys.stderr)
        sys.exit(1)

    try:
        bb_root = find_babelbrain_root(args.babelbrain_root)
    except FileNotFoundError as e:
        print(f'[run_babelbrain_thermal] {e}', file=sys.stderr)
        sys.exit(1)

    print(f'[run_babelbrain_thermal] BabelBrain root: {bb_root}')

    if bb_root not in sys.path:
        sys.path.insert(0, bb_root)

    try:
        from BabelBrain.Babel_Thermal.CalculateThermalProcess import CalculateThermalProcess
    except ImportError as e:
        print(f'[run_babelbrain_thermal] Import error: {e}', file=sys.stderr)
        sys.exit(1)

    # BabelBrain expects a dummy queue for GUI messaging; use a simple substitute
    import queue
    q = queue.Queue()

    AllDC_PRF_Duration = [{
        'DC':                           args.duty_cycle,
        'PRF':                          args.prf,
        'Duration':                     args.duration_s,
        'DurationOff':                  args.duration_off_s,
        'Repetitions':                  args.repetitions,
        'NumberGroupedSonications':     1,
        'PauseBetweenGroupedSonications': 0.0,
    }]

    print(f'[run_babelbrain_thermal] Acoustic: {args.acoustic_h5}')
    print(f'[run_babelbrain_thermal] ISPPA={args.isppa} W/cm²  '
          f'DC={args.duty_cycle:.3f}  PRF={args.prf:.0f} Hz  '
          f'Duration={args.duration_s} s')

    try:
        CalculateThermalProcess(
            queueMsg                    = q,
            case                        = args.acoustic_h5,
            AllDC_PRF_Duration          = AllDC_PRF_Duration,
            ExtraData                   = {},
            Isppa                       = args.isppa,
            sel_p                       = args.sel_p,
            BaselineTemperature         = args.baseline_temp,
            COMPUTING_BACKEND           = args.computing_backend,
            deviceName                  = args.device_name,
            Frequency                   = args.frequency,
            TxSystem                    = args.tx_system,
            LimitBHTEIterationsPerProcess = args.limit_bhte_iter,
            bConcatenateSimulations     = False,
        )
    except Exception as e:
        print(f'[run_babelbrain_thermal] CalculateThermalProcess error: {e}',
              file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Verify thermal output was written before returning success
    import glob
    acoustic_dir   = os.path.dirname(args.acoustic_h5)
    thermal_files  = glob.glob(os.path.join(acoustic_dir, '*ThermalField*.h5'))
    if not thermal_files:
        print(f'[run_babelbrain_thermal] ERROR: No *ThermalField*.h5 found in {acoustic_dir} '
              'after CalculateThermalProcess completed.', file=sys.stderr)
        sys.exit(1)
    print(f'[run_babelbrain_thermal] Thermal output: {thermal_files[-1]}')
    print('[run_babelbrain_thermal] Thermal simulation complete.')
    sys.exit(0)


if __name__ == '__main__':
    main()
