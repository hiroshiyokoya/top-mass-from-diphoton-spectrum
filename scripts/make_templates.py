#!/usr/bin/env python3
"""Generate gg -> gamma gamma M_aa templates Tab_<mt>_<gt>.dat.

Runs fortran/build/mktemplate.exe once per (m_t, Gamma_t) point, in
parallel, and writes one file per point (the format read by fit/TMDP.py).
Intended to run inside the tmdp Docker image:

    python3 scripts/make_templates.py config/templates/LHC13T.yml

The points to compute are the union of
  * the grid `masses` x `widths`,
  * every grid in `grids` (a list of {masses, widths}), and
  * the explicit pairs in `points` ([[m_t, Gamma_t], ...]).
`--check` only lists the points, without running anything.  Whether the
fit inputs find all their templates is checked on the fit side
(tests/fit/test_config.py).

A manifest.json with the full configuration is written next to the
templates so that each template set records how it was made.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import yaml

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_EXE = os.path.join(ROOT_DIR, 'fortran', 'build', 'mktemplate.exe')
TEMPLATE_RE = re.compile(r'^Tab_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.dat$')

# namelist key -> (yml key, default); defaults follow MKD_gg2aa.f
PARAMS = {
    'RS': ('rs', 100000.0),
    'MUG': ('mu_green', 40.0),
    'ETMAX': ('etamax', 2.5),
    'PTMIN': ('ptmin', 40.0),
    'PTRATIO': ('ptratio', 0.4),
    'MAAMIN': ('maa_min', 300.0),
    'MAAMAX': ('maa_max', 400.0),
    'DMAA': ('maa_step', 0.1),
    'NCALL': ('ncall', 50000),
    'ITMX': ('itmx', 6),
}
# optional integer switches, passed only when set in the config
# (defaults in fortran/src/mktemplate.f: 3, 1, 0, 0)
SWITCHES = {'MODE': 'mode', 'IORDG': 'iordg', 'NQCDIN': 'nqcd',
            'LEGACYCUT': 'legacycut'}


def repo_path(path):
    return path if os.path.isabs(path) else os.path.join(ROOT_DIR, path)


def norm_dir(path):
    """Directory as written in a yml ('./Template/X/') -> 'Template/X'."""
    return os.path.normpath(path).replace(os.sep, '/')


def expand_grid(spec):
    """A list of values, or {start, stop, step} (stop inclusive)."""
    if isinstance(spec, dict):
        n = int(round((spec['stop'] - spec['start']) / spec['step'])) + 1
        return [round(spec['start'] + i * spec['step'], 6) for i in range(n)]
    if isinstance(spec, (int, float)):
        return [float(spec)]
    return [float(v) for v in spec]


def format_width(gt):
    """1.5 -> '1.50', 1.875 -> '1.875' (as in the 2018 file names)."""
    s = '{:.3f}'.format(gt)
    return s[:-1] if s.endswith('0') else s


def template_name(mt, gt):
    return 'Tab_{:.1f}_{}.dat'.format(mt, format_width(gt))


def parse_template_name(name):
    m = TEMPLATE_RE.match(os.path.basename(name))
    if not m:
        raise ValueError('not a template file name: {}'.format(name))
    return float(m.group(1)), float(m.group(2))


def load_yaml(path):
    with open(repo_path(path)) as f:
        return yaml.safe_load(f)


def collect_points(cfg):
    """Sorted, de-duplicated (mt, gt) points requested by a template config."""
    points = set()
    grids = list(cfg.get('grids') or [])
    if 'masses' in cfg or 'widths' in cfg:
        grids.append({'masses': cfg['masses'], 'widths': cfg['widths']})
    for g in grids:
        points.update((mt, gt) for gt in expand_grid(g['widths'])
                      for mt in expand_grid(g['masses']))
    points.update((float(mt), float(gt)) for mt, gt in cfg.get('points') or [])
    # one file per name: equal names mean the same template
    by_name = {template_name(mt, gt): (mt, gt) for mt, gt in points}
    return sorted(by_name.values(), key=lambda p: (p[1], p[0]))


def namelist(cfg, mt, gt, outfile):
    lines = ['&TMPL']
    for key, (ykey, default) in PARAMS.items():
        val = cfg.get(ykey, default)
        if isinstance(default, int):
            lines.append(' {}={:d},'.format(key, int(val)))
        elif key == 'PTRATIO' and ykey not in cfg:
            continue  # keep the Fortran default (single-precision 0.4)
        else:
            lines.append(' {}={!r},'.format(key, float(val)))
    for key, ykey in SWITCHES.items():
        if ykey in cfg:
            lines.append(' {}={:d},'.format(key, int(cfg[ykey])))
    lines.append(" PDFSET='{}',".format(cfg.get('pdfset', 'CT14lo')))
    lines.append(' MT={!r}, GT={!r},'.format(float(mt), float(gt)))
    lines.append(" OUTFILE='{}',".format(outfile))
    lines.append('/')
    return '\n'.join(lines) + '\n'


def run_one(exe, cfg, mt, gt, outdir):
    outfile = os.path.join(outdir, template_name(mt, gt))
    nml = namelist(cfg, mt, gt, outfile)
    t0 = time.time()
    proc = subprocess.run([exe], input=nml, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('mt={} gt={} failed:\n{}'.format(
            mt, gt, proc.stdout[-2000:] + proc.stderr[-2000:]))
    data = np.loadtxt(outfile)
    expected = int(round((cfg.get('maa_max', 400.0)
                          - cfg.get('maa_min', 300.0))
                         / cfg.get('maa_step', 0.1))) + 1
    if data.shape != (expected, 3) or not np.all(np.isfinite(data)):
        raise RuntimeError('{}: unexpected content, shape {}'.format(
            outfile, data.shape))
    return outfile, time.time() - t0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('config', help='template config (config/templates/)')
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count(),
                        help='parallel runs (default: all CPUs)')
    parser.add_argument('--exe', default=DEFAULT_EXE)
    parser.add_argument('--force', action='store_true',
                        help='recompute templates that already exist')
    parser.add_argument('--check', action='store_true',
                        help='only list the points')
    args = parser.parse_args()

    cfg = load_yaml(args.config)
    points = collect_points(cfg)
    outdir = repo_path(cfg['outdir'])
    if args.check:
        for mt, gt in points:
            print(template_name(mt, gt))
        print('{} templates -> {}'.format(len(points), cfg['outdir']))
        return
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(args.exe):
        sys.exit('{} not found: run "make -C fortran" first'.format(args.exe))

    todo = [(mt, gt) for mt, gt in points if args.force or not
            os.path.exists(os.path.join(outdir, template_name(mt, gt)))]
    print('{} templates requested, {} to compute, {} jobs -> {}'.format(
        len(points), len(todo), args.jobs, outdir), flush=True)

    manifest = {
        'config': cfg,
        'config_file': os.path.relpath(repo_path(args.config), ROOT_DIR),
        'templates': [template_name(mt, gt) for mt, gt in points],
        'git_commit': subprocess.run(
            ['git', '-C', ROOT_DIR, 'rev-parse', 'HEAD'],
            capture_output=True, text=True).stdout.strip(),
        'started': time.strftime('%Y-%m-%dT%H:%M:%S'),
    }

    failed = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, args.exe, cfg, mt, gt, outdir):
                   (mt, gt) for mt, gt in todo}
        for fut in as_completed(futures):
            mt, gt = futures[fut]
            try:
                outfile, dt = fut.result()
                done += 1
                print('[{}/{}] {} ({:.0f} s)'.format(
                    done, len(todo), os.path.basename(outfile), dt),
                    flush=True)
            except Exception as err:  # keep the other runs going
                failed.append((mt, gt))
                print('FAILED', err, file=sys.stderr, flush=True)

    manifest['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    manifest['failed'] = failed
    with open(os.path.join(outdir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)
    missing = [p for p in points if not os.path.exists(
        os.path.join(outdir, template_name(*p)))]
    print('done: {} computed, {} failed, {} of {} present'.format(
        done, len(failed), len(points) - len(missing), len(points)))
    sys.exit(1 if failed or missing else 0)


if __name__ == '__main__':
    main()
