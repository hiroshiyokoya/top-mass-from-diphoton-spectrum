#!/usr/bin/env python3
"""Generate gg -> gamma gamma M_aa templates Tab_<mt>_<gt>.dat.

Runs fortran/build/mktemplate.exe for every (m_t, Gamma_t) point and writes
one file per point (the format read by fit/TMDP.py).  Each template is
split into chunks of m_aa points (--chunk, default 25) that run as
separate processes in parallel and are joined afterwards; chunks of an
interrupted run are reused.  Each chunk skips the random numbers of the
points before it, so chunks do not repeat each other's VEGAS sequence.
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
import shutil
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


def maa_grid(cfg):
    lo = float(cfg.get('maa_min', PARAMS['MAAMIN'][1]))
    hi = float(cfg.get('maa_max', PARAMS['MAAMAX'][1]))
    step = float(cfg.get('maa_step', PARAMS['DMAA'][1]))
    n = int(round((hi - lo) / step)) + 1
    return np.round(lo + step * np.arange(n), 6)


def namelist(cfg, mt, gt, outfile, maa_range=None, nskip=0):
    lines = ['&TMPL']
    if nskip:
        lines.append(' NSKIP={:d},'.format(int(nskip)))
    for key, (ykey, default) in PARAMS.items():
        val = cfg.get(ykey, default)
        if maa_range is not None and key == 'MAAMIN':
            val = maa_range[0]
        if maa_range is not None and key == 'MAAMAX':
            val = maa_range[1]
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


def chunk_dir(outdir, mt, gt):
    return os.path.join(outdir, '.chunks', template_name(mt, gt)[:-4])


def chunk_file(outdir, mt, gt, maa0):
    return os.path.join(chunk_dir(outdir, mt, gt), '{:09.4f}.dat'.format(maa0))


def run_chunk(exe, cfg, mt, gt, maa, outdir, nskip):
    """One mktemplate run for the m_aa points `maa` of one template, which
    start at index `nskip` of the full m_aa grid (mktemplate skips the
    random numbers of the preceding points, so chunks are independent)."""
    chunk = chunk_file(outdir, mt, gt, maa[0])
    nml = namelist(cfg, mt, gt, chunk, (maa[0], maa[-1]), nskip)
    proc = subprocess.run([exe], input=nml, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('mt={} gt={} m_aa={}-{} failed:\n{}'.format(
            mt, gt, maa[0], maa[-1],
            proc.stdout[-2000:] + proc.stderr[-2000:]))
    data = np.loadtxt(chunk, ndmin=2)
    if (data.shape != (len(maa), 3) or not np.allclose(data[:, 0], maa)
            or not np.all(np.isfinite(data))):
        os.remove(chunk)
        raise RuntimeError('{}: unexpected content'.format(chunk))


def assemble(cfg, mt, gt, outdir, chunks):
    """Join the chunk files of one template into Tab_<mt>_<gt>.dat."""
    data = np.vstack([np.loadtxt(c, ndmin=2) for c in chunks])
    if not np.allclose(data[:, 0], maa_grid(cfg)):
        raise RuntimeError('{}: chunks do not cover the m_aa grid'.format(
            template_name(mt, gt)))
    outfile = os.path.join(outdir, template_name(mt, gt))
    with open(outfile, 'w') as f:
        # the line format of mktemplate.f: (F12.4,1X,2(1PE15.5))
        for maa, value, error in data:
            f.write('{:12.4f} {:15.5E}{:15.5E}\n'.format(maa, value, error))
    shutil.rmtree(chunk_dir(outdir, mt, gt))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('config', help='template config (config/templates/)')
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count(),
                        help='parallel runs (default: all CPUs)')
    parser.add_argument('--chunk', type=int, default=25,
                        help='m_aa points per mktemplate run (default 25)')
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
    maa = maa_grid(cfg)
    size = max(1, args.chunk)
    parts = [(i, maa[i:i + size]) for i in range(0, len(maa), size)]
    jobs, chunks_of = [], {}
    for mt, gt in todo:
        if args.force and os.path.isdir(chunk_dir(outdir, mt, gt)):
            shutil.rmtree(chunk_dir(outdir, mt, gt))
        os.makedirs(chunk_dir(outdir, mt, gt), exist_ok=True)
        chunks_of[(mt, gt)] = [chunk_file(outdir, mt, gt, p[0])
                               for _, p in parts]
        jobs += [(mt, gt, i, p) for i, p in parts  # reuse finished chunks
                 if not os.path.exists(chunk_file(outdir, mt, gt, p[0]))]
    print('{} templates requested, {} to compute: {} runs of <= {} m_aa '
          'points, {} parallel -> {}'.format(
              len(points), len(todo), len(jobs), size, args.jobs, outdir),
          flush=True)

    manifest = {
        'config': cfg,
        'config_file': os.path.relpath(repo_path(args.config), ROOT_DIR),
        'templates': [template_name(mt, gt) for mt, gt in points],
        'chunk': size,
        'git_commit': subprocess.run(
            ['git', '-C', ROOT_DIR, 'rev-parse', 'HEAD'],
            capture_output=True, text=True).stdout.strip(),
        'started': time.strftime('%Y-%m-%dT%H:%M:%S'),
    }

    failed = set()
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_chunk, args.exe, cfg, mt, gt, p, outdir,
                               i): (mt, gt) for mt, gt, i, p in jobs}
        step = max(1, len(futures) // 20)
        for n, fut in enumerate(as_completed(futures), 1):
            try:
                fut.result()
            except Exception as err:  # keep the other runs going
                failed.add(futures[fut])
                print('FAILED', err, file=sys.stderr, flush=True)
            if n % step == 0 or n == len(futures):
                print('[{}/{}] runs done, {:.0f} s'.format(
                    n, len(futures), time.time() - t0), flush=True)

    done = 0
    for mt, gt in todo:
        if (mt, gt) in failed:
            continue
        try:
            assemble(cfg, mt, gt, outdir, chunks_of[(mt, gt)])
            done += 1
        except Exception as err:
            failed.add((mt, gt))
            print('FAILED', err, file=sys.stderr, flush=True)
    chunks_root = os.path.join(outdir, '.chunks')
    if os.path.isdir(chunks_root) and not os.listdir(chunks_root):
        os.rmdir(chunks_root)

    manifest['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    manifest['failed'] = sorted(failed)
    with open(os.path.join(outdir, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)
    missing = [p for p in points if not os.path.exists(
        os.path.join(outdir, template_name(*p)))]
    print('done: {} computed, {} failed, {} of {} present'.format(
        done, len(failed), len(points) - len(missing), len(points)))
    sys.exit(1 if failed or missing else 0)


if __name__ == '__main__':
    main()
