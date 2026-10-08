"""
Python module for "Top-quark Mass from DiPhoton mass spectrum"

Hiroshi Yokoya <hyokoya@gmail.com>
"""
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

def get_mass_width_from_filename(fname):
    separate = fname.split('_')
    mass = float(separate[1])
    width = float(separate[2][0:-4])
    return mass, width

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
        self.hGenBG.FillRandom(self.hBG ,int(self.Nbg) )
        self.hGenSig.FillRandom(self.hSig,int(self.Nsig))
        self.hGen = self.hGenSig + self.hGenBG
        self.hGen.SetTitle('Gen')
        self.hGen.SetName('Gen')

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
    LHC = TMDP('yaml/fit/input_LHC13T.yml')

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

    results = []
    list_best = []

    for i in range(Nloop):
        LHC.genEvents()
        result = []
        for mass, pfnc in zip(LHC.list_mass, LHC.list_pfnc):
            LHC.hGen.Fit(pfnc, LHC.fit_options)
            chi2 = pfnc.GetChisquare()
            result.append(chi2)
        best = LHC.list_mass[result.index(min(result))]
        list_best.append(best)
        results.append(result)
        if(i % (Nloop/10) == 0):
            print('{0}/{1}: current average {2}'.format( \
                   i, Nloop, np.mean(list_best)))
    # End of Nloop

    last_time = time.time()
    print("time for fitting: {0}".format(last_time-middle_time) + "[sec]")

    print(np.mean(list_best), np.std(list_best))

    with open('output.dat','w') as fout:
        fout.write(str(list_best))

    min_mass = min(LHC.list_mass)
    max_mass = max(LHC.list_mass)
    bin_widt = 1
    nbin = (max_mass-min_mass)/bin_widt

    # Draw a histogram of the best-fit top-quark mass
    plt.figure(1)
    plt.hist(list_best, bins=int(nbin), range=(min_mass,max_mass))
    plt.show()

