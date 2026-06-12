"""Rewrite a BabelBrain DataForSim.h5 file replacing blosc compression with gzip.

MATLAB's built-in HDF5 library cannot read blosc-compressed datasets unless the
blosc filter plugin is registered at library load time.  This script reads the
file using h5py (which has blosc available in the BabelBrain Python env) and
rewrites every dataset with gzip-4 compression, which MATLAB supports natively.

Usage:
    python rewrite_hdf5_no_blosc.py <h5_path>

The original file is replaced in-place (written to a temp file then moved).
"""

import sys
import os
import tempfile
import shutil
import hdf5plugin  # registers blosc and other filters with h5py
import h5py
import numpy as np


def rewrite(src_path):
    dir_name = os.path.dirname(src_path)
    fd, tmp_path = tempfile.mkstemp(suffix='.h5', dir=dir_name)
    os.close(fd)
    try:
        with h5py.File(src_path, 'r') as src, h5py.File(tmp_path, 'w') as dst:
            _copy_group(src, dst)
        shutil.move(tmp_path, src_path)
        print(f'[rewrite_hdf5_no_blosc] Rewrote {src_path} (blosc → gzip-4)', flush=True)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def _copy_group(src_grp, dst_grp):
    for name, item in src_grp.items():
        if isinstance(item, h5py.Group):
            child = dst_grp.require_group(name)
            _copy_group(item, child)
        elif isinstance(item, h5py.Dataset):
            data = item[()]
            kwargs = {
                'data': data,
                'compression': 'gzip',
                'compression_opts': 4,
            }
            if data.ndim == 0 or data.size == 0:
                # Scalar / empty — no chunking/compression
                kwargs = {'data': data}
            dst_grp.create_dataset(name, **kwargs)
            # Copy attributes
            for k, v in item.attrs.items():
                dst_grp[name].attrs[k] = v
    # Copy group attributes
    for k, v in src_grp.attrs.items():
        dst_grp.attrs[k] = v


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(f'Usage: {sys.argv[0]} <h5_path>')
    rewrite(sys.argv[1])
