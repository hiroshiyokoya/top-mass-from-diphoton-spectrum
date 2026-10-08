#!/usr/bin/env python3
"""Reproduce the gg -> gamma gamma signal of arXiv:1607.00990.

Computes every curve of Figs. 1 and 4 of the paper with
fortran/build/mktemplate.exe (in parallel, split into m_aa chunks),
compares them with the curves extracted from the paper's figures
(reference/1607.00990/), and writes

  <outdir>/<figure>_<curve>.dat   m_aa, dsigma/dm_aa [fb/GeV], error
  <outdir>/<figure>.png           our curves (lines) over the paper's (dots)
  <outdir>/summary.md, summary.json

    python3 scripts/reproduce_1607_00990.py yaml/repro/1607.00990.yml -j 20

Exit status 0 if every curve agrees with the paper within `tolerance`.
Existing chunk results are reused unless --force is given.
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
INT_KEYS = {'MODE', 'IORDG', 'NQCDIN', 'LEGACYCUT', 'NCALL', 'ITMX'}
STR_KEYS = {'PDFSET', 'OUTFILE'}


def repo_path(path):
    return path if os.path.isabs(path) else os.path.join(ROOT_DIR, path)


def namelist(values):
    lines = ['&TMPL']
    for key, val in values.items():
        if key in STR_KEYS:
            lines.append(" {}='{}',".format(key, val))
        elif key in INT_KEYS:
            lines.append(' {}={:d},'.format(key, int(val)))
        else:
            lines.append(' {}={!r},'.format(key, float(val)))
    lines.append('/')
    return '\n'.join(lines) + '\n'


def grid(spec):
    n = int(round((spec['max'] - spec['min']) / spec['step'])) + 1
    return np.round(spec['min'] + spec['step'] * np.arange(n), 6)


def jobs_for(cfg, outdir):
    """(curve key, chunk file, namelist values) for every chunk."""
    jobs = []
    for fig, fcfg in cfg['figures'].items():
        maa = grid(fcfg['grid'])
        for curve, ccfg in fcfg['curves'].items():
            key = '{}_{}'.format(fig, curve)
            for i0 in range(0, len(maa), cfg['chunk']):
                part = maa[i0:i0 + cfg['chunk']]
                values = dict(cfg['common'])
                values.update(fcfg.get('set', {}))
                values.update(ccfg)
                values.update({'NCALL': cfg['ncall'], 'ITMX': cfg['itmx'],
                               'MAAMIN': part[0], 'MAAMAX': part[-1],
                               'DMAA': fcfg['grid']['step']})
                chunk = os.path.join(outdir, 'chunks',
                                     '{}_{:04d}.dat'.format(key, i0))
                values['OUTFILE'] = chunk
                jobs.append((key, chunk, values, len(part)))
    return jobs


def run_job(exe, chunk, values, npt):
    proc = subprocess.run([exe], input=namelist(values),
                          capture_output=True, text=True)
    data = np.loadtxt(chunk, ndmin=2) if os.path.exists(chunk) else None
    if proc.returncode != 0 or data is None or data.shape != (npt, 3):
        raise RuntimeError('{} failed:\n{}'.format(
            chunk, (proc.stdout + proc.stderr)[-2000:]))
    return chunk


def extrema(x, y, lo, hi, reach=4.0):
    """Dip and bump: the local minimum in [lo, hi] followed by the largest
    rise to a local maximum within `reach` GeV.  None if there is no dip."""
    sel = (x >= lo) & (x <= hi)
    xs, ys = x[sel], y[sel]
    best = None
    for i in range(1, len(ys) - 1):
        if not (ys[i] <= ys[i - 1] and ys[i] < ys[i + 1]):
            continue
        k = np.nonzero((xs > xs[i]) & (xs <= xs[i] + reach))[0]
        j = k[np.argmax(ys[k])]
        rise = ys[j] - ys[i]
        if j < len(ys) - 1 and (best is None or rise > best[0]):
            best = (rise, xs[i], xs[j])
    return None if best is None else best[1:]


def compare(cfg, outdir):
    rows = []
    for fig, fcfg in cfg['figures'].items():
        for curve in fcfg['curves']:
            key = '{}_{}'.format(fig, curve)
            ours = np.loadtxt(os.path.join(outdir, key + '.dat'))
            ref = np.loadtxt(os.path.join(repo_path(cfg['reference_dir']),
                                          key + '.csv'), delimiter=',')
            ref = ref[(ref[:, 0] >= ours[0, 0]) & (ref[:, 0] <= ours[-1, 0])]
            pred = np.interp(ref[:, 0], ours[:, 0], ours[:, 1])
            rel = pred / ref[:, 1] - 1.0
            row = {'curve': key, 'n_ref': int(len(ref)),
                   'max_abs_rel_dev': float(np.max(np.abs(rel))),
                   'mean_rel_dev': float(np.mean(rel)),
                   'at': float(ref[np.argmax(np.abs(rel)), 0])}
            mt = fcfg.get('set', {}).get('MT', cfg['common'].get('MT'))
            mt = fcfg['curves'][curve].get('MT', mt)
            if fcfg['curves'][curve].get('MODE', 3) != 0:
                lo, hi = 2 * mt - 8, 2 * mt + 4
                for name, arr in (('ours', ours), ('paper', ref)):
                    db = extrema(arr[:, 0], arr[:, 1], lo, hi)
                    row['dip_bump_' + name] = (None if db is None else
                                               [float(v) for v in db])
            rows.append(row)
    return rows


def plot(cfg, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for fig, fcfg in cfg['figures'].items():
        f, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                    gridspec_kw={'height_ratios': [3, 1]})
        for idx, curve in enumerate(fcfg['curves']):
            key = '{}_{}'.format(fig, curve)
            ours = np.loadtxt(os.path.join(outdir, key + '.dat'))
            ref = np.loadtxt(os.path.join(repo_path(cfg['reference_dir']),
                                          key + '.csv'), delimiter=',')
            color = 'C{}'.format(idx)
            ax.plot(ours[:, 0], ours[:, 1], color=color, lw=1.2,
                    label=curve)
            ax.plot(ref[::4, 0], ref[::4, 1], 'o', color=color, ms=2.2,
                    alpha=0.6)
            sel = (ref[:, 0] >= ours[0, 0]) & (ref[:, 0] <= ours[-1, 0])
            axr.plot(ref[sel, 0], np.interp(ref[sel, 0], ours[:, 0],
                                            ours[:, 1]) / ref[sel, 1] - 1,
                     color=color, lw=1)
        ax.set_title(fcfg['title'] + '\nlines: this repository, '
                     'dots: arXiv:1607.00990', fontsize=10)
        ax.set_ylabel(r'$d\sigma/dm_{\gamma\gamma}$ [fb/GeV]')
        ax.legend(fontsize=8)
        tol = cfg['tolerance']
        axr.axhspan(-tol, tol, color='0.9')
        axr.axhline(0, color='0.5', lw=0.5)
        axr.set_ylim(-2 * tol, 2 * tol)
        axr.set_ylabel('ours/paper - 1')
        axr.set_xlabel(r'$m_{\gamma\gamma}$ [GeV]')
        f.tight_layout()
        f.savefig(os.path.join(outdir, fig + '.png'), dpi=130)
        plt.close(f)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('config')
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    parser.add_argument('--exe', default=DEFAULT_EXE)
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()

    with open(repo_path(args.config)) as f:
        cfg = yaml.safe_load(f)
    outdir = repo_path(cfg['outdir'])
    os.makedirs(os.path.join(outdir, 'chunks'), exist_ok=True)
    jobs = jobs_for(cfg, outdir)
    todo = [j for j in jobs if args.force or not os.path.exists(j[1])]
    print('{} chunks, {} to compute, {} parallel'.format(
        len(jobs), len(todo), args.jobs), flush=True)

    t0 = time.time()
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(run_job, args.exe, chunk, values, npt): chunk
                for _, chunk, values, npt in todo}
        for n, fut in enumerate(as_completed(futs), 1):
            try:
                fut.result()
            except Exception as err:
                failed.append(futs[fut])
                print('FAILED', err, file=sys.stderr, flush=True)
            if n % 20 == 0 or n == len(futs):
                print('[{}/{}] {:.0f} s'.format(n, len(futs),
                                                time.time() - t0), flush=True)
    if failed:
        sys.exit('{} chunks failed'.format(len(failed)))

    # assemble one file per curve
    curves = {}
    for key, chunk, _, _ in jobs:
        curves.setdefault(key, []).append(np.loadtxt(chunk, ndmin=2))
    for key, parts in curves.items():
        np.savetxt(os.path.join(outdir, key + '.dat'),
                   np.vstack(parts), fmt='%10.4f %14.6e %14.6e')

    rows = compare(cfg, outdir)
    plot(cfg, outdir)
    ok = all(r['max_abs_rel_dev'] <= cfg['tolerance'] for r in rows)
    lines = ['# Reproduction of arXiv:1607.00990 (signal)', '',
             'Config: `{}`, tolerance {:.1%}. Result: **{}**'.format(
                 args.config, cfg['tolerance'], 'PASS' if ok else 'FAIL'),
             '',
             '| curve | max abs. rel. dev. | at m_aa [GeV] | mean rel. dev. '
             '| dip / bump, ours [GeV] | dip / bump, paper [GeV] |',
             '|---|---|---|---|---|---|']
    for r in rows:
        db = lambda v: '{:.2f} / {:.2f}'.format(*v) if v else '-'
        lines.append('| {} | {:.3%} | {:.1f} | {:+.3%} | {} | {} |'.format(
            r['curve'], r['max_abs_rel_dev'], r['at'], r['mean_rel_dev'],
            db(r.get('dip_bump_ours')), db(r.get('dip_bump_paper'))))
    with open(os.path.join(outdir, 'summary.md'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump({'pass': ok, 'tolerance': cfg['tolerance'], 'curves': rows},
                  f, indent=1)
    print('\n'.join(lines))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
