import os
import pathlib
import numpy as np

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import fea 
import graph

num_runs = 5 # number of runs

L = 10000 # total length of fiber

amplitudes = [L/10,L/20,L/40]

num_kappa_tilde = 4

runs = []
for amplitude in amplitudes:
    for run in range(num_runs):
        for kappa_tilde in range(num_kappa_tilde):
            runs.append(f'python3 simulate_discretized_sin.py {int(amplitude)} {run} {kappa_tilde}')

os.system(runs[int(sys.argv[1])-1])