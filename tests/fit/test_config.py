"""The fit inputs (fit/config/) must find their templates.

Every fit input reads its templates from a directory `dir`; exactly one
template config in config/templates/ writes that directory, with the same
sqrt(s), and its points include every template the fit input reads.
No Fortran is run.
"""
import glob
import os
import sys

import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT_DIR, 'scripts'))
import make_templates as mt  # noqa: E402


def rel(paths):
    return sorted(os.path.relpath(p, ROOT_DIR).replace(os.sep, '/')
                  for p in paths)


FIT_INPUTS = rel(glob.glob(os.path.join(ROOT_DIR, 'fit', 'config', '*.yml')))
TEMPLATE_CONFIGS = rel(glob.glob(os.path.join(ROOT_DIR, 'config', 'templates',
                                              '*.yml')))


def test_layout():
    assert len(FIT_INPUTS) == 8
    assert not glob.glob(os.path.join(ROOT_DIR, '*.yml'))


@pytest.mark.parametrize('fit', FIT_INPUTS)
def test_fit_input_templates_available(fit):
    fcfg = mt.load_yaml(fit)
    owners = [c for c in TEMPLATE_CONFIGS
              if mt.norm_dir(mt.load_yaml(c)['outdir'])
              == mt.norm_dir(fcfg['dir'])]
    assert len(owners) == 1, owners
    cfg = mt.load_yaml(owners[0])
    assert float(cfg['rs']) == float(fcfg['rs'])
    names = {mt.template_name(*p) for p in mt.collect_points(cfg)}
    wanted = set(fcfg['files_sig']) | set(fcfg['files_template'])
    assert wanted <= names, sorted(wanted - names)


@pytest.mark.parametrize('name', ['Tab_173.0_1.50.dat', 'Tab_173.0_1.875.dat',
                                  'Tab_172.5_0.25.dat', 'Tab_160.0_4.00.dat'])
def test_template_name_roundtrip(name):
    assert mt.template_name(*mt.parse_template_name(name)) == name
