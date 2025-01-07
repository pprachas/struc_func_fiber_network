import os
import pathlib
import numpy as np

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import fea 
import graph

num_runs = 20 # number of runs

num_links = [10,20,30,40,50,60, 70, 80] # number of links

L = 10000 # total length of fiber
widths = [L/10,L/20,L/40,L/50]

num_kappa_tilde = 4

runs = []
for num_link in num_links:
    for w in widths:
        for run in range(num_runs):
            for kappa_tilde in range(num_kappa_tilde):
                runs.append(f'python3 simulate_random_chain.py {num_link} {run} {w} {kappa_tilde}')


os.system(runs[int(sys.argv[1])-1])
                
