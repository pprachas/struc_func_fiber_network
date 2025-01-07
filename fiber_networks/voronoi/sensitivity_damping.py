import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys


#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
sys.path.append('..')

import fea_HR as fea
import graph
from generate_voronoi import generate_network

L = 10000
W=H=L
char_length = L

n = 500
seed = 1
mesh_threshold = 0.0

kappa_tilde = 1e-6

init_c0 = 1e-7 #[1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10][int(sys.argv[1])-1]
max_res = [1e-2, 5e-3,1e-3, 5e-4, 1e-4, 5e-5,1e-5, 5e-6,1e-6][int(sys.argv[1])-1] #[1e-2, 5e-3,1e-3, 5e-4, 1e-4, 5e-5][int(sys.argv[1])-1]
num_segments = 15


# Name of mesh
mesh_name= f'sensitivity/mesh/n{n}/voronoi{seed}.xdmf'
f_name = f'plots/sensitivity/damping/max_res{max_res}'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

#generate_network(W,H,n,seed,char_length,num_segments,mesh_threshold, root_dir = 'sensitivity')

# run mesh and find critical point
# Construct graphs
nodes = np.loadtxt(f'./graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'./graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)
step_size = L/1000.

crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['r','r'], abs_tol = 1e-10, rel_tol = 1e-12, init_c0 = init_c0, max_damping_res = max_res)

np.savetxt(f'{f_name}/crit_strain.txt', np.array([crit_strain]))

print('Starting FEA run!')

step_size = crit_strain*L/100.

fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], 1.5*crit_strain, compute_individual=False, max_damping_res = max_res, compute_residual = True, init_c0 = init_c0, abs_tol = 1e-10, rel_tol = 1e-12) # make runs stop after strain stiffening