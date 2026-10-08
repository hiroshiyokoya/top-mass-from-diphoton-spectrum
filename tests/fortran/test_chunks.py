"""Chunked template generation reproduces a single run (NSKIP).

With NCALL = 5000 VEGAS makes exactly 5000 calls per iteration, and at this
low statistics it never stops early, so skipping NSKIP points' worth of
random numbers puts the generator in the same state as a single run: the
last point computed alone equals the last point of a run over all points.
"""
import os
import subprocess
import sys

import numpy as np

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT_DIR, 'scripts'))
import make_templates as mt  # noqa: E402

CFG = {'rs': 13000., 'ncall': 5000, 'itmx': 3, 'maa_step': 1.0}


def run(tmp_path, name, maa_range, nskip):
    out = str(tmp_path / name)
    subprocess.run([mt.DEFAULT_EXE],
                   input=mt.namelist(CFG, 173.0, 1.5, out, maa_range, nskip),
                   capture_output=True, text=True, check=True)
    return np.loadtxt(out, ndmin=2)


def test_chunk_equals_single_run(tmp_path):
    single = run(tmp_path, 'single.dat', (340.0, 342.0), 0)
    assert single.shape == (3, 3)
    alone = run(tmp_path, 'alone.dat', (342.0, 342.0), 0)
    chunk = run(tmp_path, 'chunk.dat', (342.0, 342.0), 2)
    assert np.array_equal(chunk[0], single[2])
    # without skipping, the point gets another random sequence
    assert not np.array_equal(alone[0], single[2])
