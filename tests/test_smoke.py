"""End-to-end smoke test: Fortran template generator -> TMDP fit.

Uses the low-statistics templates of yaml/templates/test.yml (generated
on first use, ~1 min) and *synthetic* exponential backgrounds.  It checks
that the chain runs, not that the physics is right.

    python3 -m pytest -q tests        # inside the pytmdp Docker image
"""
import os
import pathlib
import subprocess
import sys

import numpy as np
import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_TEMPLATES = os.path.join(ROOT_DIR, 'Template', 'test')
MASSES = ['172.0', '173.0']


@pytest.fixture(scope='module')
def templates():
    names = ['Tab_{}_1.50.dat'.format(m) for m in MASSES]
    if not all(os.path.exists(os.path.join(TEST_TEMPLATES, n))
               for n in names):
        subprocess.run([sys.executable,
                        os.path.join(ROOT_DIR, 'scripts', 'make_templates.py'),
                        os.path.join(ROOT_DIR, 'yaml', 'templates', 'test.yml')],
                       check=True)
    return names


def test_template_format(templates):
    for name in templates:
        data = np.loadtxt(os.path.join(TEST_TEMPLATES, name))
        assert data.shape == (1001, 3)
        assert np.allclose(data[:, 0], 300.0 + 0.1 * np.arange(1001))
        assert np.all(data[:, 1] > 0)


def make_tmdp(templates):
    # TMDP.get_mass_width_from_filename splits the *full path* on '_'
    # (docs/REVIEW.md P6), so the work directory must not contain '_'.
    tmp_path = pathlib.Path('/tmp/tmdpsmoke{}'.format(os.getpid()))
    tmp_path.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    for name in ['Direct.dat', 'OneF.dat', 'TwoF.dat']:
        np.savetxt(tmp_path / name, 300 + rng.exponential(60, 20000) % 100)
    for name in templates:
        link = tmp_path / name
        if not link.exists():
            os.symlink(os.path.join(TEST_TEMPLATES, name), link)
    yml = tmp_path / 'input.yml'
    yml.write_text('\n'.join([
        'dir: {}/'.format(tmp_path), 'rs: 13000', 'lum: 3000', 'corr: 1.2',
        'kgg: 0.1', 'Nevnt: 515000', 'hbin: 100', 'fitopt: LMNQR',
        'fmin: 330', 'fmax: 360', 'sig_dir: 0.0717', 'sig_one: 0.0105',
        'sig_two: 0.0001', 'files_dir: [Direct.dat]',
        'files_one: [OneF.dat]', 'files_two: [TwoF.dat]',
        'files_sig: [{}]'.format(templates[1]),
        'files_template: [{}]'.format(', '.join(templates))]) + '\n')
    sys.path.insert(0, ROOT_DIR)
    import TMDP
    tmdp = TMDP.TMDP(str(yml))
    tmdp.read_template()
    return tmdp


def test_fit_chain(templates):
    tmdp = make_tmdp(templates)
    assert tmdp.list_mass == [172.0, 173.0]
    tmdp.genEvents()
    for pfnc in tmdp.list_pfnc:
        tmdp.hGen.Fit(pfnc, tmdp.fit_options)
        assert np.isfinite(pfnc.GetChisquare())


@pytest.mark.xfail(strict=True, reason='docs/REVIEW.md P1: with ROOT 6.34 '
                   'TH1::FillRandom(TH1*, n) silently fills nothing when the '
                   'binnings differ and n > ~10 x nbins')
def test_pseudo_data_not_empty(templates):
    tmdp = make_tmdp(templates)
    tmdp.genEvents()
    assert tmdp.hGen.Integral() == pytest.approx(
        int(tmdp.Nsig) + int(tmdp.Nbg))
