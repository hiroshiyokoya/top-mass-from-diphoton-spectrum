"""
> python3 -i ScanMass.py yaml/fit/input_LHC13T.yml Nloop
"""
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

import TMDP

start_time = time.time()

args = sys.argv

## initialize TMDP class with an input yml file ##
if len(args) < 2:
    yml = input('input yml file:')
else:
    yml = args[1]

if len(args) > 2:
    Nloop = int(args[2])
else:
    Nloop = 1

## initialize template functions for fitting ##
LHC = TMDP.TMDP(yml)
LHC.read_template()

middle_time = time.time()
print("time for setup: {0}".format(middle_time-start_time) + "[sec]")

results = []
list_best = []

# Nloop of pseudo-experiments
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

with open('outMass.dat','w') as fout:
    fout.write(str(list_best))

min_mass = min(LHC.list_mass)
max_mass = max(LHC.list_mass)
bin_widt = 1
nbin = (max_mass-min_mass)/bin_widt

# Draw a histogram of the best-fit top-quark mass
plt.figure(1)
plt.hist(list_best, bins=int(nbin), range=(min_mass,max_mass))
plt.show()
