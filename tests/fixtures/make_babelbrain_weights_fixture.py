"""Write BabelBrain-layout CALIBRATION fixtures for test_babelbrain_weights.m.

BabelBrain's TxCalibration.py saves its fit with
    SaveToH5py({'CALIBRATION': res.x, 'deviceName': ..., ...}, 'CALIBRATION.h5')
(BabelViscoFDTD.H5pySimple). res.x is a 1-D complex numpy array with one
weight per transducer SUB-element. h5py stores it as a compound {r, i}; every
root item gets a string attribute 'type'. Arrays above 1 KiB are Blosc
compressed, so the fixtures stay at or below 1 KiB (uncompressed, readable by
MATLAB).

The real SaveToH5py is used when BabelViscoFDTD is importable; otherwise the
branches of its ProcType used here are mirrored.

Usage:
    python3 make_babelbrain_weights_fixture.py [out_dir]

Writes:
    bb_calibration_ring.h5     - 4 weights (one per ring), complex128
    bb_calibration_subelem.h5  - 64 weights (4 rings x 16 sub-elements), complex64

Ring phases [deg]: 0, 30, -45, 120; ring amplitudes: 1, 0.9, 1.1, 1.
Sub-element phases scatter symmetrically (+/-5 deg) about the ring phase, so
the phase of each ring's complex sum equals the ring phase.
"""

import os
import sys

import h5py
import numpy as np

RING_PHASE_DEG = np.array([0.0, 30.0, -45.0, 120.0])
RING_AMP = np.array([1.0, 0.9, 1.1, 1.0])
N_SUB_PER_RING = 16

try:
    from BabelViscoFDTD.H5pySimple import SaveToH5py
except ImportError:
    def _proc(group, k, v):
        if isinstance(v, dict):
            g = group.create_group(k)
            g.attrs['type'] = 'dict'
            g.attrs['type_key'] = 'str'
            for kk, vv in v.items():
                _proc(g, kk, vv)
        elif isinstance(v, list):
            g = group.create_group(k)
            g.attrs['type'] = 'list'
            for n, vv in enumerate(v):
                _proc(g, 'item_%d' % n, vv)
        elif isinstance(v, np.ndarray):
            assert v.nbytes <= 2**10, 'mirror writes uncompressed arrays only'
            ds = group.create_dataset(k, data=v)
            ds.attrs['type'] = 'ndarray'
        elif isinstance(v, str):
            ds = group.create_dataset(k, data=v)
            ds.attrs['type'] = str(type(v))
        elif np.isscalar(v):
            ds = group.create_dataset(k, data=v)
            ds.attrs['type'] = 'scalar'
        else:
            raise TypeError(type(v).__name__)

    def SaveToH5py(d, fname):
        with h5py.File(fname, 'w') as f:
            for k, v in d.items():
                _proc(f, k, v)


def cal_results(x):
    # Subset of TxCalibration.py's CalResults with the same Python types;
    # the values are illustrative.
    return {'CALIBRATION': x,
            'deviceName': 'test',
            'Frequency': 500e3,
            'FocalLength': 63.2e-3,
            'InDiameters': np.array([0.0, 32.9, 46.0, 55.7]) * 1e-3,
            'OutDiameters': np.array([32.3, 45.4, 55.2, 64.0]) * 1e-3,
            'Locations': [30, 50, 70],
            'Regularizer': 'None',
            'bUseRayleighPhase': True}


def main(out_dir):
    ring = (RING_AMP * np.exp(1j * np.deg2rad(RING_PHASE_DEG))).astype(np.complex128)
    jitter = np.deg2rad(np.linspace(-5, 5, N_SUB_PER_RING))
    sub = np.concatenate([RING_AMP[n] * np.exp(1j * (np.deg2rad(RING_PHASE_DEG[n]) + jitter))
                          for n in range(RING_PHASE_DEG.size)]).astype(np.complex64)
    SaveToH5py(cal_results(ring), os.path.join(out_dir, 'bb_calibration_ring.h5'))
    SaveToH5py(cal_results(sub), os.path.join(out_dir, 'bb_calibration_subelem.h5'))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
