"""
> python3 -i fit/Scan2D.py fit/config/input2D_FCC100.yml
"""
import sys
import time

import TMDP

start_time = time.time()

args = sys.argv

## initialize TMDP class with an input yml file ##
if len(args) < 2:
    yml = input('input yml file:')
else:
    yml = args[1]

## initialize template functions for fitting ##
LHC = TMDP.TMDP(yml)
LHC.read_template()

middle_time = time.time()
print("time for setup: {0}".format(middle_time-start_time) + "[sec]")

results = []
list_best = []

LHC.genEvents()
result = []
for mass, width, pfnc in zip(LHC.list_mass, LHC.list_width, LHC.list_pfnc):
    LHC.hGen.Fit(pfnc, LHC.fit_options)
    chi2 = pfnc.GetChisquare()
    print(mass, width, chi2)
    result.append([mass, width, chi2])

last_time = time.time()

print ("time for fitting:{0}".format(last_time-middle_time) + "[sec]")

with open('out2D.dat','w') as fout:
    for res in result:
        fout.write('{0} {1} {2}\n'.format(res[0], res[1], res[2]))

