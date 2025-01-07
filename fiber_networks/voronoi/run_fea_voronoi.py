import os
import pathlib
import numpy as np
import networkx as nx

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import graph

num_runs = 20 # number of runs

num_seed = [100,200,300,400,500]# number of links

L = 10000 # total length of fiber
num_kappa_tilde = 4

combined_list = []

for kappa_tilde in range(num_kappa_tilde):
    for n in num_seed:
        for run in range(num_runs):
            combined_list.append([n, run, kappa_tilde])

combined_list=[combined_list[ii-1] for ii in [387]]
args = combined_list[int(sys.argv[1])-1]

n = args[0]

print(f'python3 fea_voronoi.py {args[0]} {args[1]} {args[2]}')

os.system(f'python3 fea_voronoi.py {args[0]} {args[1]} {args[2]}')

