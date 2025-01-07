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

L = 10000
n = 300
seed = 13

crit_strain_all = []
init_contour = []
kappa_tilde = 1e-4
num_seg = 10 #int(np.loadtxt(f'mesh_refinement/mesh/n{n}/final_char{seed}.txt'))

# Name of mesh
mesh_name=f'mesh/n{n}/voronoi{seed}.xdmf'
f_name = f'visualize_network'
# f_name = f'plots/n{n}/kappa{kappa_tilde}/seed{seed}/data/quad3'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

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
step_size = L/2500.

fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], 1546./L, abs_tol = 1e-10, rel_tol = 1e-12, init_c0 = 1e-7) # make runs stop after strain stiffening