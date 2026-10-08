#!/usr/bin/env python3
"""Generate gg -> gamma gamma M_aa templates Tab_<mt>_<gt>.dat.

Runs fortran/build/mktemplate.exe once per (m_t, Gamma_t) point, in
parallel, and writes the files read by TMDP.GG2AA.  Intended to run
inside the pytmdp Docker image:

    python3 scripts/make_templates.py config/templates_LHC13T.yml

A manifest.json with the full configuration is written next to the
templates so that each template set records how it was made.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import yaml

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_EXE = os.path.join(ROOT_DIR, 'fortran', 'build', 'mktemplate.exe')

# namelist key -> (yml key, default); defaults follow MKD_gg2aa.f
PARAMS = {
    'RS': ('rs', 100000.0),
    'MUG': ('mu_green', 40.0),
    'ETMAX': ('etamax', 2.5),
    'PTMIN': ('ptmin', 40.0),
    'MAAMIN': ('maa_min', 300.0),
    'MAAMAX': ('maa_max', 400.0),
    'DMAA': ('maa_step', 0.1),
    'NCALL': ('ncall', 50000),
    'ITMX': ('itmx', 6),
}


def expand_grid(spec):
    """A list of values, or {start, stop, step} (stop inclusive)."""
    if isinstance(spec, dict):
        n = int(round((spec['stop'] - spec['start']) / spec['step'])) + 1
        return [round(spec['start'] + i * spec['step'], 6) for i in range(n)]
    if isinstance(spec, (int, float)):
        return [float(spec)]
    return [float(v) for v in spec]


def template_name(mt, gt):
    # Same as the Fortran format '("Tab_",F5.1,"_",F4.2,".dat")'
    return 'Tab_{:5.1f}_{:4.2f}.dat'.format(mt, gt)


def namelist(cfg, mt, gt, outfile):
    lines = ['&TMPL']
    for key, (ykey, default) in PARAMS.items():
        val = cfg.get(ykey, default)
        if isinstance(default, int):
            lines.append(' {}={:d},'.format(key, int(val)))
        else:
            lines.append(' {}={!r},'.format(key, float(val)))
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
    parser.add_argument('config', help='yml file (see config/)')
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count(),
                        help='parallel runs (default: all CPUs)')
    parser.add_argument('--exe', default=DEFAULT_EXE)
    parser.add_argument('--force', action='store_true',
                        help='recompute templates that already exist')
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    outdir = cfg['outdir']
    if not os.path.isabs(outdir):
        outdir = os.path.join(ROOT_DIR, outdir)
    os.makedirs(outdir, exist_ok=True)
    if not os.path.exists(args.exe):
        sys.exit('{} not found: run "make -C fortran" first'.format(args.exe))

    points = [(mt, gt) for gt in expand_grid(cfg['widths'])
              for mt in expand_grid(cfg['masses'])]
    todo = [(mt, gt) for mt, gt in points if args.force or not
            os.path.exists(os.path.join(outdir, template_name(mt, gt)))]
    print('{} templates requested, {} to compute, {} jobs -> {}'.format(
        len(points), len(todo), args.jobs, outdir), flush=True)

    manifest = {
        'config': cfg,
        'config_file': os.path.abspath(args.config),
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
