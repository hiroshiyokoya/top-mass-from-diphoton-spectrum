"""Spot checks of the signal of arXiv:1607.00990 (a few seconds per point).

Each point is computed with the full VEGAS statistics of the paper setup
(yaml/repro/1607.00990.yml) and compared with the curve extracted from the
paper's figure.  The full comparison of all curves is
scripts/reproduce_1607_00990.py.
"""
import os
import subprocess
import sys

import numpy as np
import pytest
import yaml

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT_DIR, 'scripts'))
import reproduce_1607_00990 as repro  # noqa: E402

CFG = yaml.safe_load(open(os.path.join(ROOT_DIR, 'yaml', 'repro',
                                       '1607.00990.yml')))
TOL = 0.005

POINTS = [  # (figure, curve, m_aa)
    ('fig1R', 'GNLO_mu40', 330.0),
    ('fig1R', 'GNLO_mu40', 342.9),   # dip
    ('fig1R', 'GNLO_mu40', 344.7),   # bump
    ('fig1R', 'one_loop', 350.0),
    ('fig1L', 'GLO_mu20', 345.0),
    ('fig4R', 'mt173', 300.0),       # FCC, pT > 0.4 m_aa
]


@pytest.mark.parametrize('fig,curve,maa', POINTS)
def test_point(fig, curve, maa, tmp_path):
    fcfg = CFG['figures'][fig]
    values = dict(CFG['common'])
    values.update(fcfg.get('set', {}))
    values.update(fcfg['curves'][curve])
    out = tmp_path / 'o.dat'
    values.update({'NCALL': CFG['ncall'], 'ITMX': CFG['itmx'],
                   'MAAMIN': maa, 'MAAMAX': maa, 'OUTFILE': str(out)})
    subprocess.run([repro.DEFAULT_EXE], input=repro.namelist(values),
                   capture_output=True, text=True, check=True)
    ours = np.loadtxt(out)[1]
    ref = np.loadtxt(os.path.join(ROOT_DIR, CFG['reference_dir'],
                                  '{}_{}.csv'.format(fig, curve)),
                     delimiter=',')
    paper = np.interp(maa, ref[:, 0], ref[:, 1])
    assert ours == pytest.approx(paper, rel=TOL)
