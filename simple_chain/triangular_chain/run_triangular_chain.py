import os
import pathlib
import numpy as np
import networkx as nx

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import fea 
import graph

num_runs = 4 # number of runs

L = 10000 # total length of fiber

amplitude = [L/10,L/20,L/40][int(sys.argv[1])-1] # number of links

num_kappa_tilde = 4


init_list=np.empty([num_runs,num_kappa_tilde]) 
init_list.fill('NaN') 

pathlib.Path('./phase_diagram/crit_strain/').mkdir(parents=True, exist_ok=True)
pathlib.Path('./phase_diagram/init_cont/').mkdir(parents=True, exist_ok=True)

# file names to save results
f_crit_name = f'phase_diagram/crit_strain/a{int(amplitude)}.txt'
f_cont_name = f'phase_diagram/init_cont/a{int(amplitude)}.txt'

# file names to save results
f_crit_name = f'phase_diagram/crit_strain/a{int(amplitude)}.txt'
f_cont_name = f'phase_diagram/init_cont/a{int(amplitude)}.txt'

np.savetxt(f_crit_name, init_list)
np.savetxt(f_cont_name, init_list)

for run in range(num_runs):
    for kappa_tilde in range(num_kappa_tilde):
        os.system(f'python3 triangular_chain.py {amplitude} {run} {kappa_tilde}')