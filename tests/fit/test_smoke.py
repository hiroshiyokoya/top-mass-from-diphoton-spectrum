"""End-to-end smoke test: Fortran template generator -> TMDP fit.

Uses the low-statistics templates of config/templates/test.yml (generated
on first use, ~1 min) and *synthetic* exponential backgrounds.  It checks
that the chain runs, not that the physics is right.

    python3 -m pytest -q tests        # inside the tmdp Docker image
"""
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
TEST_TEMPLATES = os.path.join(ROOT_DIR, 'Template', 'test')
MASSES = ['172.0', '173.0']


@pytest.fixture(scope='module')
def templates():
    names = ['Tab_{}_1.50.dat'.format(m) for m in MASSES]
    if not all(os.path.exists(os.path.join(TEST_TEMPLATES, n))
               for n in names):
        subprocess.run([sys.executable,
                        os.path.join(ROOT_DIR, 'scripts', 'make_templates.py'),
                        os.path.join(ROOT_DIR, 'config', 'templates', 'test.yml')],
                       check=True)
    return names


def test_template_format(templates):
    for name in templates:
        data = np.loadtxt(os.path.join(TEST_TEMPLATES, name))
        assert data.shape == (1001, 3)
        assert np.allclose(data[:, 0], 300.0 + 0.1 * np.arange(1001))
        assert np.all(data[:, 1] > 0)


def make_tmdp(templates, tmp_path):
    # tmp_path from pytest contains '_' (e.g. .../test_fit_chain0), which
    # broke the template-name parsing before docs/REVIEW.md P6 was fixed.
    rng = np.random.default_rng(1)
    for name in ['Direct.dat', 'OneF.dat', 'TwoF.dat']:
        np.savetxt(tmp_path / name, 300 + rng.exponential(60, 20000) % 100)
    for name in templates:
        os.symlink(os.path.join(TEST_TEMPLATES, name), tmp_path / name)
    yml = tmp_path / 'input.yml'
    yml.write_text('\n'.join([
        'dir: {}/'.format(tmp_path), 'rs: 13000', 'lum: 3000', 'corr: 1.2',
        'kgg: 0.1', 'Nevnt: 515000', 'hbin: 100', 'fitopt: LMNQR',
        'fmin: 330', 'fmax: 360', 'sig_dir: 0.0717', 'sig_one: 0.0105',
        'sig_two: 0.0001', 'seed: 12345', 'files_dir: [Direct.dat]',
        'files_one: [OneF.dat]', 'files_two: [TwoF.dat]',
        'files_sig: [{}]'.format(templates[1]),
        'files_template: [{}]'.format(', '.join(templates))]) + '\n')
    sys.path.insert(0, os.path.join(ROOT_DIR, 'fit'))
    import TMDP
    tmdp = TMDP.TMDP(str(yml))
    tmdp.read_template()
    return tmdp


def test_fit_chain(templates, tmp_path):
    assert '_' in str(tmp_path)
    tmdp = make_tmdp(templates, tmp_path)
    assert tmdp.list_mass == [172.0, 173.0]
    tmdp.genEvents()
    for pfnc in tmdp.list_pfnc:
        tmdp.hGen.Fit(pfnc, tmdp.fit_options)
        assert np.isfinite(pfnc.GetChisquare())


def test_template_name_parsing():
    sys.path.insert(0, os.path.join(ROOT_DIR, 'fit'))
    import TMDP
    assert TMDP.get_mass_width_from_filename(
        '/a_b/Template/LHC13T_x/Tab_173.0_1.875.dat') == (173.0, 1.875)
    with pytest.raises(ValueError):
        TMDP.get_mass_width_from_filename('/tmp/Tab_173.0.dat')


def test_pseudo_data(templates, tmp_path):
    """docs/REVIEW.md P1/P4: the pseudo-data are Poisson(mu_i) around the
    expected counts, which add up to Nsig + Nbg."""
    tmdp = make_tmdp(templates, tmp_path)
    tmdp.genEvents()
    n_exp = tmdp.Nsig + tmdp.Nbg
    assert tmdp.mu_sig.sum() + tmdp.mu_bg.sum() == pytest.approx(n_exp)
    assert len(tmdp.mu_sig) == tmdp.hbin and np.all(tmdp.mu_sig > 0)
    counts = np.array([tmdp.hGen.GetBinContent(i)
                       for i in range(1, tmdp.hbin + 1)])
    assert abs(counts.sum() - n_exp) < 5 * np.sqrt(n_exp)
    mu = tmdp.mu_sig + tmdp.mu_bg
    chi2 = np.sum((counts - mu) ** 2 / mu)
    assert abs(chi2 - tmdp.hbin) < 5 * np.sqrt(2 * tmdp.hbin)
    # the signal shape is that of the template: compare with direct binning
    data = np.loadtxt(os.path.join(TEST_TEMPLATES, templates[1]))
    edges = np.linspace(300., 400., tmdp.hbin + 1)
    direct = np.array([np.trapz(data[(data[:, 0] >= lo - 1e-9)
                                     & (data[:, 0] <= hi + 1e-9), 1], dx=0.1)
                       for lo, hi in zip(edges[:-1], edges[1:])])
    assert tmdp.mu_sig / tmdp.mu_sig.sum() == pytest.approx(
        direct / direct.sum(), rel=0.02)
