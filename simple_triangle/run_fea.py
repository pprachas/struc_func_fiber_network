import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys


#--------Import functions---------------------#
import sys
sys.path.append('../utils')

import fea_HR as fea

kappa_tilde = 1e-6
L = 10000

# Name of mesh
mesh_name= f'mesh/simple_triangle.xdmf'
f_name = f'plots/sensitivity/no_reg'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

#generate_network(W,H,n,seed,char_length,num_segments,mesh_threshold, root_dir = 'sensitivity')

# run mesh and find critical point
l_outer = (L/2)*np.sqrt(2)
l_inner = L

L_char = L
r = 2*L_char*np.sqrt(kappa_tilde)


step_size = L/50
fea.run_fea(mesh_name, f_name, r, L, step_size, ['p','p'], 0.5, compute_individual=False, compute_residual = True, init_c0 = None, abs_tol = 1e-10, rel_tol = 1e-12) # make runs stop after strain stiffening