"""
> python3 -i fit/ScanWidth.py fit/config/inputWidth_LHC13T.yml Nloop
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

## initialize TMDP class with an input yml file ##
LHC = TMDP.TMDP(yml)
LHC.read_template()

middle_time = time.time()
print("time for setup: {0}".format(middle_time-start_time) + "[sec]")

# Nloop of pseudo-experiments; the best value of each is the vertex of a
# parabola fitted to chi^2 near its minimum (TMDP.parabola_minimum)
rows = TMDP.run_scan(LHC, LHC.list_width, Nloop, 'outWidth.dat')

last_time = time.time()
print("time for fitting: {0}".format(last_time-middle_time) + "[sec]")

# Draw a histogram of the best-fit top-quark width (parabola fits)
plt.figure(1)
plt.hist([r[1] for r in rows if r[3] == 'ok'], bins=40)
plt.show()
