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
args = combined_list[int(sys.argv[1])-1]

n = args[0]

init_list=np.empty([num_runs, num_kappa_tilde]) 
init_list.fill('NaN') 

# Make sure path exists
pathlib.Path(f'./phase_diagram/crit_strain/').mkdir(parents=True, exist_ok=True)

# file names to save results
f_crit_name = f'phase_diagram/crit_strain/n{int(n)}.txt'

if os.path.isfile(f_crit_name):
    pass
else:
    np.savetxt(f_crit_name, init_list)

print(f'python3 crit_strain_voronoi.py {args[0]} {args[1]} {args[2]}')

os.system(f'python3 crit_strain_voronoi.py {args[0]} {args[1]} {args[2]}')