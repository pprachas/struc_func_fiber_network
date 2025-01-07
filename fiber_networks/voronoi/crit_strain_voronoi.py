import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys
import time

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
sys.path.append('..')

import fea_HR as fea
import graph
from generate_voronoi import generate_network

t0 = time.time()
n = int(sys.argv[1]) # number of seed
L = 10000 # total length of fiber
W=H=L
char_length = L
mesh_threshold= 0 
num_segments = 15

crit_strain_all = []
init_contour = []
seed = int(sys.argv[2])
num_kappa_tilde = int(sys.argv[3])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][num_kappa_tilde]

max_res = 1e-4

# Name of mesh
mesh_name= f'mesh/n{n}/voronoi{seed}.xdmf'

# file names to write results
f_crit_name = f'phase_diagram/crit_strain/n{int(n)}.txt'
# f_cont_name = f'phase_diagram/init_cont/n{n}.txt'

# find critical point

nodes = np.loadtxt(f'./graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'./graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)

generate_network(W,H,n,seed,num_segments,mesh_threshold)

step_size = L/1000.
crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['r','r'], abs_tol = 1e-10, rel_tol = 1e-12, init_c0 = 1e-7, max_damping_res = max_res)
print(crit_strain)
# Save critical strain values
crit_array = np.loadtxt(f_crit_name)
crit_array[seed][num_kappa_tilde] = crit_strain

print(f'time: {time.time()-t0}')