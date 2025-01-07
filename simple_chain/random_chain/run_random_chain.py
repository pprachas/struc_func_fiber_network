import os
import pathlib
import numpy as np
import networkx as nx

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import graph

num_runs = 20 # number of runs

num_links = [10,20,30,40,50,60, 70, 80][int(sys.argv[1])-1] # number of links

L = 10000 # total length of fiber
widths = [L/10,L/20,L/40,L/50]
char_length = np.loadtxt('mesh_refinement/final_char.txt', dtype = int) # extra careful mesh refinement
num_kappa_tilde = 4

init_list=np.empty([num_runs,num_kappa_tilde]) 
init_list.fill('NaN') 

for w in widths:
    # Make sure path exists
    pathlib.Path(f'./phase_diagram/crit_strain/w{int(w)}').mkdir(parents=True, exist_ok=True)
    pathlib.Path(f'./phase_diagram/init_cont/w{int(w)}').mkdir(parents=True, exist_ok=True)

    # file names to save results
    f_crit_name = f'phase_diagram/crit_strain/w{int(w)}/n{num_links}.txt'
    f_cont_name = f'phase_diagram/init_cont/w{int(w)}/n{num_links}.txt'

    np.savetxt(f_crit_name, init_list)
    np.savetxt(f_cont_name, init_list)

for w in widths:
    for run in range(num_runs):
        for kappa_tilde in range(num_kappa_tilde):
            os.system(f'python3 random_chain.py {num_links} {run} {w} {kappa_tilde}')