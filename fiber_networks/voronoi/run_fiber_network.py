import os
import pathlib
import numpy as np
import networkx as nx

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import fea 
import graph

num_runs = 30 # number of runs

num_seeds = [100,200,300,400][int(sys.argv[1])-1] # number of links

L = 10000 # total length of fiber
char_length = np.loadtxt('mesh_refinement/final_char.txt', dtype = int) # extra careful mesh refinement
num_kappa_tilde = 4
num_paths = 5

init_list=np.empty([num_runs,num_kappa_tilde])
init_list.fill('NaN') 

init_list_path = np.empty([num_runs,num_paths])
init_list_path.fill('NaN') 
pathlib.Path(f'./phase_diagram/crit_strain/').mkdir(parents=True, exist_ok=True)
pathlib.Path(f'./phase_diagram/init_cont/').mkdir(parents=True, exist_ok=True)

# file names to save results
f_crit_name = f'phase_diagram/crit_strain/n{int(num_seeds)}.txt'
f_cont_name = f'phase_diagram/init_cont/n{int(num_seeds)}.txt'

np.savetxt(f_crit_name, init_list)
np.savetxt(f_cont_name, init_list_path)

for run in range(num_runs):
    for kappa_tilde in range(num_kappa_tilde):
        os.system(f'python3 fiber_network.py {num_seeds} {run} {kappa_tilde}')