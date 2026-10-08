"""
Python module for "Top-quark Mass from DiPhoton mass spectrum"

Hiroshi Yokoya <hyokoya@gmail.com>
"""
import os
import re
import sys
import time
import yaml

import matplotlib.pyplot as plt
import numpy as np

import ROOT

def read_yaml(fname):
    print('input yml file :',fname)
    with open(fname,'r') as fopen:
        text = fopen.read()
    dict_ = yaml.safe_load(text)
    return dict_

def calc_simpson_integral(data):
    n = len(data[0])
    h = data[0][1] - data[0][0]
    ilist = list(range(0,n-1,2))
    S_part = [h/3.*(data[1][i] + 4*data[1][i+1] + data[1][i+2]) \
               for i in ilist]
    return sum(S_part)

TEMPLATE_NAME = re.compile(r'^Tab_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.dat$')

def get_mass_width_from_filename(fname):
    # Parse the file name only: the directory may contain '_'
    # (docs/REVIEW.md P6).
    match = TEMPLATE_NAME.match(os.path.basename(fname))
    if match is None:
        raise ValueError('not a template file name (Tab_<mt>_<Gt>.dat): '
                         + fname)
    return float(match.group(1)), float(match.group(2))

def parabola_minimum(values, chi2, window=4.0):
    """Best-fit value from chi^2 on the template grid (docs/REVIEW.md P2).

    A parabola is fitted to the chi^2 of the grid point with the smallest
    chi^2 and of its neighbours: always the two nearest ones, and further
    ones on each side while chi^2 - chi^2_min <= window.  Only the chi^2
    values are fitted; templates are never interpolated or extrapolated.

    Returns (best, sigma, status):
      'ok'    best = vertex of the parabola, sigma = its Delta chi^2 = 1
              half-width
      'edge'  the smallest chi^2 is at the first or last grid point:
              best = that grid point, sigma = nan
      'flat'  the parabola does not open upwards, or its vertex lies outside
              the fitted points: best = grid minimum, sigma = nan
    """
    values = np.asarray(values, dtype=float)
    chi2 = np.asarray(chi2, dtype=float)
    order = np.argsort(values)
    x, y = values[order], chi2[order]
    i0 = int(np.argmin(y))
    if i0 == 0 or i0 == len(x) - 1:
        return x[i0], np.nan, 'edge'
    lo, hi = i0 - 1, i0 + 1
    while lo > 0 and y[lo - 1] - y[i0] <= window:
        lo -= 1
    while hi < len(x) - 1 and y[hi + 1] - y[i0] <= window:
        hi += 1
    a, b, _ = np.polyfit(x[lo:hi + 1], y[lo:hi + 1], 2)
    if a <= 0:
        return x[i0], np.nan, 'flat'
    best = -b / (2 * a)
    if not x[lo] <= best <= x[hi]:
        return x[i0], np.nan, 'flat'
    return best, 1 / np.sqrt(a), 'ok'

def run_scan(tmdp, values, nloop, outfile):
    """Repeat pseudo-experiments and fit every template to each of them.

    `values` is the scanned quantity of each template (m_t or Gamma_t).
    Per pseudo-experiment the grid minimum and the parabola result of
    parabola_minimum are recorded, written to `outfile` and summarised.
    Returns a list of (grid_best, best, sigma, status).
    """
    rows = []
    for i in range(nloop):
        tmdp.genEvents()
        chi2 = tmdp.fit_templates()
        grid_best = values[int(np.argmin(chi2))]
        rows.append((grid_best,) + parabola_minimum(values, chi2))
        if (i + 1) % max(1, nloop // 10) == 0:
            ok = [r[1] for r in rows if r[3] == 'ok']
            print('{0}/{1}: current average {2}'.format(
                i + 1, nloop, np.mean(ok) if ok else np.nan))
    with open(outfile, 'w') as fout:
        fout.write('# grid_best  best  sigma  status   (parabola fit to chi^2;'
                   ' see TMDP.parabola_minimum)\n')
        for row in rows:
            fout.write('{0} {1:.4f} {2:.4f} {3}\n'.format(*row))
    summarize_scan(rows)
    return rows

def summarize_scan(rows):
    status = [r[3] for r in rows]
    ok = np.array([r[1:3] for r in rows if r[3] == 'ok'], dtype=float)
    grid = np.array([r[0] for r in rows], dtype=float)
    print('pseudo-experiments: {0} (ok {1}, edge {2}, flat {3})'.format(
        len(rows), status.count('ok'), status.count('edge'),
        status.count('flat')))
    if len(ok):
        print('parabola fit: mean {0:.4f}, spread {1:.4f}, '
              'mean Delta chi2 = 1 error {2:.4f}'.format(
                  ok[:, 0].mean(), ok[:, 0].std(), ok[:, 1].mean()))
    print('grid minimum: mean {0:.4f}, spread {1:.4f}'.format(
        grid.mean(), grid.std()))

def bin_fractions(hist, edges):
    """Fraction of the content of the ROOT histogram `hist` in each bin
    [edges[i], edges[i+1]], with the content spread uniformly inside each
    bin of `hist`.  Works for any binning of `hist`."""
    n = hist.GetNbinsX()
    hist_edges = np.array([hist.GetBinLowEdge(i) for i in range(1, n + 2)])
    content = np.array([hist.GetBinContent(i) for i in range(1, n + 1)])
    cumulative = np.concatenate([[0.], np.cumsum(content)])
    frac = np.diff(np.interp(edges, hist_edges, cumulative))
    return frac / frac.sum()

def read_template_from_file(fname):
    with open(fname,'r') as fopen:
        raw_data = fopen.read()
    list_data = list(map(float,raw_data.split()))
    np_data = np.array(list_data).reshape((-1,3)).T
    data = np_data.tolist()
    return data

def read_events_from_file(fname):
    with open(fname,'r') as fopen:
        raw_data = fopen.read()
    list_data = list(map(float,raw_data.split()))
    np_data = np.array(list_data)
    data = np_data.tolist()
    return data

class TMDP(object):
    """
    Python module for "Top-quark Mass from Diphoton Mass Spectrum"

    Hiroshi Yokoya <hyokoya@gmail.com>
    """

    bin_ = 200
    min_ = 300
    max_ = 400

    def __init__(self,yml):

        init = read_yaml(yml)

        # These parameters should not be modified, unless you calculate
        # all the template and background events by yourself.
        self.hmin = 300.
        self.hmax = 400.

        # Set the template directory
        self.dir = self.set_init('dir',init)

        self.rs = self.set_init('rs',init)
        self.lum = self.set_init('lum',init)
        self.corr = self.set_init('corr',init)
        self.kgg = self.set_init('kgg',init)

        self.sig_dir = self.set_init('sig_dir',init)
        self.sig_one = self.set_init('sig_one',init)
        self.sig_two = self.set_init('sig_two',init)
        self.sig_bg = self.sig_dir + self.sig_one + self.sig_two

        self.Nevnt = self.set_init('Nevnt',init)
        self.Nsig = self.Nevnt * self.kgg
        self.Nbg = self.Nevnt * (1 - self.kgg)

        self.hbin = self.set_init('hbin',init)
        self.fit_options = self.set_init('fitopt',init)
        self.fmin = self.set_init('fmin',init)
        self.fmax = self.set_init('fmax',init)

        self.files_dir = self.set_init('files_dir',init)
        self.files_one = self.set_init('files_one',init)
        self.files_two = self.set_init('files_two',init)
        self.files_sig = self.set_init('files_sig',init)
        self.files_temp = self.set_init('files_template',init)

        # Random numbers for the pseudo-data; optional 'seed' in the yml
        # makes a run reproducible.
        self.rng = np.random.default_rng(init.get('seed'))

        self.read_background()
        self.read_signal()

    def set_init(self,param,dict):
        if(param == 'Nevnt'):
            try:
                item = dict['Nevnt']
            except KeyError:
                item = self.sig_bg * self.lum * 10**3 \
                       * self.corr / (1-self.kgg)
            print('Nevnt :',item)
        else:
            try:
                item = dict[param]
                print(param,':',item)
            except KeyError:
                print(param,'is not set in the yml file.')
                item = None
        return item

    def read_background(self):
        self.dir_data = []
        for file in self.files_dir:
            self.fdir = self.dir + file
            self.data = read_events_from_file(self.fdir)
            self.dir_data.extend(self.data)
        self.h_dir = ROOT.TH1D('dir','dir',self.bin_,self.min_,self.max_)
        for i in self.dir_data:
            self.h_dir.Fill(i)
        self.h_dir.Scale(self.sig_dir / self.h_dir.Integral())

        self.one_data = []
        for file in self.files_one:
            self.fone = self.dir + file
            self.data = read_events_from_file(self.fone)
            self.one_data.extend(self.data)
        self.h_one = ROOT.TH1D('one','one',self.bin_,self.min_,self.max_)
        for i in self.one_data:
            self.h_one.Fill(i)
        self.h_one.Scale(self.sig_one / self.h_one.Integral())

        self.two_data = []
        for file in self.files_two:
            self.ftwo = self.dir + file
            self.data = read_events_from_file(self.ftwo)
            self.two_data.extend(self.data)
        self.h_two = ROOT.TH1D('two','two',self.bin_,self.min_,self.max_)
        for i in self.two_data:
            self.h_two.Fill(i)
        self.h_two.Scale(self.sig_two / self.h_two.Integral())

        self.hBG = self.h_dir + self.h_one + self.h_two
        print('read background data')

    def read_signal(self):
        self.hSig = ROOT.TH1D('sig','sig',1001,299.95,400.05)
        for file in self.files_sig:
            self.fsig = self.dir + file
            self.data = read_template_from_file(self.fsig)
            for i in range(len(self.data[0])):
                x, y = self.data[0][i], self.data[1][i]
                self.hSig.Fill(x,y)
        print('read signal data')

    def genEvents(self):
        """Fill hGenBG, hGenSig and hGen with one pseudo-dataset.

        The expected numbers of events per bin, mu_i, are the background
        (hBG) and signal (hSig) shapes integrated over the hbin bins of
        [hmin, hmax] and scaled to Nbg and Nsig; the observed numbers are
        Poisson(mu_i).  This replaces TH1::FillRandom, which with ROOT 6.34
        fills nothing when the binnings differ (docs/REVIEW.md P1), and lets
        the total number of events fluctuate (P4).
        """
        if not hasattr(self, 'mu_bg'):
            edges = np.linspace(self.hmin, self.hmax, self.hbin + 1)
            self.mu_bg = self.Nbg * bin_fractions(self.hBG, edges)
            self.mu_sig = self.Nsig * bin_fractions(self.hSig, edges)
        try:
            self.hGenBG.Reset()
        except AttributeError:
            self.hGenBG = ROOT.TH1D('Gen BG','Gen BG', \
                                    self.hbin,self.hmin,self.hmax)
            self.hGenBG.SetLineColor(25)
            self.hGenBG.SetFillColor(25)
        try:
            self.hGenSig.Reset()
        except AttributeError:
            self.hGenSig = ROOT.TH1D('Gen Sig','Gen Sig', \
                                     self.hbin,self.hmin,self.hmax)
            self.hGenSig.SetLineColor(38)
            self.hGenSig.SetFillColor(38)
        for hist, mu in ((self.hGenBG, self.mu_bg),
                         (self.hGenSig, self.mu_sig)):
            counts = self.rng.poisson(mu)
            for i, n in enumerate(counts, 1):
                hist.SetBinContent(i, n)
                hist.SetBinError(i, np.sqrt(n))
            hist.SetEntries(counts.sum())
        self.hGen = self.hGenSig + self.hGenBG
        self.hGen.SetTitle('Gen')
        self.hGen.SetName('Gen')

    def fit_templates(self):
        """Fit the current pseudo-data (hGen) with every template;
        returns the chi^2 of each fit, in the order of files_template."""
        chi2 = []
        for pfnc in self.list_pfnc:
            self.hGen.Fit(pfnc, self.fit_options)
            chi2.append(pfnc.GetChisquare())
        return np.array(chi2)

    def read_template(self):
        self.list_fname = [self.dir + temp for temp in self.files_temp]
        self.list_temp = [ GG2AA(fname, self.rs, self.hmin, self.hmax, \
                    self.hbin, self.Nevnt) for fname in self.list_fname ]

        self.list_mass = [ temp.mass for temp in self.list_temp ]
        self.list_width = [ temp.width for temp in self.list_temp ]

        self.list_pfnc = [ ROOT.TF1('', temp.PDF, self.fmin, self.fmax, 2) \
                           for temp in self.list_temp ]
        for pfnc in self.list_pfnc:
            pfnc.SetParNames('A0', 'Kgg')
            pfnc.SetParLimits(0, -100., 100.)
            pfnc.SetParLimits(1, 0., 1.)
            pfnc.SetParameters(50., .5)
        print('read template data')


class GG2AA(object):
    """
    Template for gg2aa distribution for given mass & width, 
    under certain theory setup (PDF,scales, etc.).
    """
    def __init__(self,fname,rs,hmin,hmax,hbin,Nevnt):
        self.fname = fname
        self.mass, self.width = get_mass_width_from_filename(fname)
        self.data = read_template_from_file(self.fname)
        self.norm = calc_simpson_integral(self.data)

        self.rs = rs
        self.hmin = hmin
        self.hmax = hmax
        self.hbin = hbin
        self.Nevnt = Nevnt

        self.rmax = self.hmax / self.rs
        self.rmin = self.hmin / self.rs
        self.rmax3 = pow(self.rmax,1/3)
        self.rmin3 = pow(self.rmin,1/3)

    def interpolate_func(self,x):
        try:
            i = int(10*(x-299.95))
            f = self.data[1][i]
        except IndexError:
            print('x: out of range')
            f = None
        return f

    def func(self,x):
        f = self.interpolate_func(x)
        return f

    def func_ATLAS(self,r,a):
        self.sATL = ( 3*pow(1-self.rmin3,1+a) * (2 + (1+a) \
            * (2 + (2+a)*self.rmin3) * self.rmin3) \
            - 3*pow(1-self.rmax3,1+a) * (2 + (1+a) * (2 + (2+a) \
            * self.rmax3) * self.rmax3) ) / ((1+a)*(2+a)*(3+a))
        ATL = pow((1-pow(r,1/3)),a) / self.sATL
        return ATL

    def PDF(self,x,par):
        r = x[0] / self.rs
        fgg = self.func(x[0]) / self.norm
        f = self.Nevnt/self.hbin*(self.hmax-self.hmin) \
            * ( (1-par[1])/self.rs*self.func_ATLAS(r,par[0]) \
                + par[1] * fgg )
        return f


if __name__ == '__main__':

    start_time = time.time()

    ## initialize TMDP class with an input yml file ##
    LHC = TMDP('fit/config/input_LHC13T.yml')

    #LHC.hBG.Draw()
    #LHC.hSig.Draw()

    ## generate events sample and put them in TH1D.hGen ##
    #LHC.genEvents()
    #LHC.hGen.Draw()

    ## initialize template functions for fitting ##
    LHC.read_template()

    middle_time = time.time()
    print("time for setup: {0}".format(middle_time-start_time) + "[sec]")

    Nloop = 1000

    rows = run_scan(LHC, LHC.list_mass, Nloop, 'output.dat')

    last_time = time.time()
    print("time for fitting: {0}".format(last_time-middle_time) + "[sec]")

    # Draw a histogram of the best-fit top-quark mass (parabola fits)
    plt.figure(1)
    plt.hist([r[1] for r in rows if r[3] == 'ok'], bins=40)
    plt.show()
