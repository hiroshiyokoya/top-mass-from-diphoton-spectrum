#!/usr/bin/env python3
"""Extract the curves of the d sigma / d m_aa figures of arXiv:1607.00990.

The figures in the arXiv source are vector PDFs written by Grace.  They are
converted to SVG with pdftocairo (poppler-utils), and the stroked curves are
mapped to data coordinates with the major tick marks of each frame.  The
result is exact up to the resolution of the PDF (well below 0.1% of the
axis range), not a hand digitisation.

    python3 scripts/tools/digitize_1607_00990.py SRC_DIR OUT_DIR

SRC_DIR holds the figure PDFs of the arXiv e-print of 1607.00990
(LHC13TeV_LO_NLOPDF.pdf, LHC13TeV_NLO_NLOPDF.pdf, LHC13TeV_Temp_NLOPDF.pdf,
FCC100TeV_Temp_NLOPDF.pdf).  One CSV per curve is written to OUT_DIR.
"""
import os
import re
import subprocess
import sys

BLACK = 'rgb(0%, 0%, 0%)'
RED = 'rgb(100%, 0%, 0%)'
BLUE = 'rgb(0%, 0%, 100%)'
GREEN = 'rgb(0%, 54.492188%, 0%)'
MAGENTA = 'rgb(100%, 0%, 100%)'
DASH = {  # Grace line styles used in Fig. 1
    'one_loop': '160.344 68.544',
    'mu20': '22.644 68.544 114.444 68.544 22.644 68.544',
    'mu40': '22.644 68.544',
    'mu80': '114.444 68.544',
    'mu160': None,
}

# figure -> (x range, y major tick values bottom-to-top, curves)
#   curve: (csv name, stroke colour, dash pattern)
FIGURES = {
    'LHC13TeV_NLO_NLOPDF': (
        (330., 360.), [0.14 + 0.02 * i for i in range(8)],
        [('fig1R_one_loop', BLACK, DASH['one_loop'])] +
        [('fig1R_GNLO_' + k, RED, DASH[k])
         for k in ('mu20', 'mu40', 'mu80', 'mu160')]),
    'LHC13TeV_LO_NLOPDF': (
        (330., 360.), [0.14 + 0.02 * i for i in range(8)],
        [('fig1L_one_loop', BLACK, DASH['one_loop'])] +
        [('fig1L_GLO_' + k, BLUE, DASH[k])
         for k in ('mu20', 'mu40', 'mu80', 'mu160')]),
    'LHC13TeV_Temp_NLOPDF': (
        (300., 400.), [0.1 * i for i in range(6)],
        [('fig4L_mt{}'.format(m), c, None) for m, c in
         ((167, BLACK), (170, RED), (173, BLUE), (176, GREEN),
          (179, MAGENTA))]),
    'FCC100TeV_Temp_NLOPDF': (
        (300., 400.), [1.0 * i for i in range(7)],
        [('fig4R_mt{}'.format(m), c, None) for m, c in
         ((167, BLACK), (170, RED), (173, BLUE), (176, GREEN),
          (179, MAGENTA))]),
}


def svg_paths(fname):
    text = open(fname, encoding='utf-8').read()
    text = text[text.index('</defs>'):]  # skip glyph definitions
    out = []
    for m in re.finditer(r'<path ([^>]*)/>', text):
        attrs = ' ' + m.group(1)
        d = re.search(r' d="([^"]*)"', attrs).group(1)
        stroke = re.search(r'stroke="([^"]*)"', attrs)
        dash = re.search(r'stroke-dasharray="([^"]*)"', attrs)
        tr = re.search(r'transform="matrix\(([^)]*)\)"', attrs)
        toks = re.findall(r'[MLCZ]|-?[\d.]+', d)
        pts, cmd, i = [], None, 0
        while i < len(toks):
            if toks[i] in 'MLCZ':
                cmd = toks[i]
                i += 1
                continue
            n = 6 if cmd == 'C' else 2
            vals = [float(v) for v in toks[i:i + n]]
            i += n
            pts.append((vals[-2], vals[-1]))
        if tr:
            a, b, c, e, f, g = (float(v) for v in tr.group(1).split(','))
            pts = [(a * x + c * y + f, b * x + e * y + g) for x, y in pts]
        out.append({'stroke': stroke.group(1) if stroke else None,
                    'dash': dash.group(1) if dash else None, 'pts': pts})
    return out


def calibrate(paths, xrange_, yticks):
    frame = [p['pts'] for p in paths
             if p['stroke'] == BLACK and len(p['pts']) == 5][0]
    x0 = min(q[0] for q in frame)
    x1 = max(q[0] for q in frame)
    # major y ticks: horizontal segments of length ~12.2 on the left edge
    major = sorted({round(p['pts'][0][1], 3) for p in paths
                    if len(p['pts']) == 2
                    and abs(p['pts'][0][0] - x0) < 0.2
                    and abs(p['pts'][0][1] - p['pts'][1][1]) < 0.01
                    and 11. < abs(p['pts'][1][0] - p['pts'][0][0]) < 13.},
                   reverse=True)  # SVG y grows downwards
    if len(major) != len(yticks):
        raise RuntimeError('found {} major y ticks, expected {}'.format(
            len(major), len(yticks)))
    ya, yb = major[0], major[-1]
    va, vb = yticks[0], yticks[-1]

    def to_data(pt):
        x = xrange_[0] + (pt[0] - x0) / (x1 - x0) * (xrange_[1] - xrange_[0])
        y = va + (pt[1] - ya) / (yb - ya) * (vb - va)
        return x, y
    return to_data


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    for fig, (xrange_, yticks, curves) in FIGURES.items():
        svg = os.path.join(out, fig + '.svg')
        subprocess.run(['pdftocairo', '-svg', os.path.join(src, fig + '.pdf'),
                        svg], check=True)
        paths = svg_paths(svg)
        os.remove(svg)
        to_data = calibrate(paths, xrange_, yticks)
        for name, stroke, dash in curves:
            cand = [p for p in paths if p['stroke'] == stroke
                    and p['dash'] == dash and len(p['pts']) > 50]
            if len(cand) != 1:
                raise RuntimeError('{} {}: {} candidate paths'.format(
                    fig, name, len(cand)))
            pts = sorted(to_data(q) for q in cand[0]['pts'])
            with open(os.path.join(out, name + '.csv'), 'w') as f:
                f.write('# arXiv:1607.00990, figure {}.pdf, curve {}\n'
                        '# m_aa [GeV], dsigma/dm_aa [fb/GeV]\n'.format(
                            fig, name))
                for x, y in pts:
                    f.write('{:.4f},{:.6f}\n'.format(x, y))
            print('{}: {} points'.format(name, len(pts)))


if __name__ == '__main__':
    main()
