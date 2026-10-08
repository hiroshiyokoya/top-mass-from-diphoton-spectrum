"""yaml/fit and yaml/templates must agree.

Every fit input is listed by exactly one template config, whose outdir and
sqrt(s) match the fit input and whose points include all templates the fit
input reads.  No Fortran is run.
"""
import glob
import os
import sys

import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT_DIR, 'scripts'))
import make_templates as mt  # noqa: E402

FIT_INPUTS = sorted(os.path.relpath(p, ROOT_DIR).replace(os.sep, '/')
                    for p in glob.glob(os.path.join(ROOT_DIR, 'yaml', 'fit',
                                                    '*.yml')))
TEMPLATE_CONFIGS = sorted(os.path.relpath(p, ROOT_DIR).replace(os.sep, '/')
                          for p in glob.glob(os.path.join(
                              ROOT_DIR, 'yaml', 'templates', '*.yml')))


def owners(fit):
    return [c for c in TEMPLATE_CONFIGS
            if fit in (mt.load_yaml(c).get('fit_inputs') or [])]


def test_layout():
    assert len(FIT_INPUTS) == 8
    assert not glob.glob(os.path.join(ROOT_DIR, '*.yml'))


@pytest.mark.parametrize('fit', FIT_INPUTS)
def test_fit_input_has_one_template_config(fit):
    assert len(owners(fit)) == 1


@pytest.mark.parametrize('config', TEMPLATE_CONFIGS)
def test_template_config_consistent(config):
    cfg = mt.load_yaml(config)
    assert mt.check_config(cfg) == []
    names = {mt.template_name(*p) for p in mt.collect_points(cfg)}
    for fit in cfg.get('fit_inputs') or []:
        fcfg = mt.load_yaml(fit)
        wanted = list(fcfg['files_sig']) + list(fcfg['files_template'])
        assert set(wanted) <= names, sorted(set(wanted) - names)


@pytest.mark.parametrize('name', ['Tab_173.0_1.50.dat', 'Tab_173.0_1.875.dat',
                                  'Tab_172.5_0.25.dat', 'Tab_160.0_4.00.dat'])
def test_template_name_roundtrip(name):
    assert mt.template_name(*mt.parse_template_name(name)) == name
