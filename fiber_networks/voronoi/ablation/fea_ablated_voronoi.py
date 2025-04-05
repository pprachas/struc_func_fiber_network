import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys


#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
sys.path.append('..')

import fea_HR as fea
import graph

L = 10000

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
path_num = int(sys.argv[4])
num_edge = int(sys.argv[5])
edge_num = int(sys.argv[6])

max_damping = 5.0e-5

# Name of mesh
mesh_name= f'mesh/n{n}/path{path_num}/num_edge{num_edge}/voronoi{seed}_edge{edge_num}.xdmf'

# file names to save results
f_name = f'fea_run/ablated/n{n}/path{path_num}/num_edge{num_edge}/edge{edge_num}/kappa_tilde{kappa_tilde}/seed{seed}'

print(f_name)
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

# run mesh and find critical point
nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)

crit_strains = np.loadtxt(f'../phase_diagram/crit_strain/n{int(n)}.txt')
step_size = crit_strains[seed][num_kappa_tilde]*L/100.

fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], crit_strains[seed][num_kappa_tilde]*1.5, compute_individual = True, abs_tol = 1e-10, rel_tol = 1e-12, init_c0 = 1e-5, max_damping_res = max_damping) # make runs stop after strain stiffening