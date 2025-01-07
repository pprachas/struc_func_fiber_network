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

n = 100
seed = 0
mesh_threshold = 0.0

kappa_tilde = 1e-6

init_c0 = 1e-7
max_damping = [1e-2, 5e-3,1e-3, 5e-4, 1e-4, 5e-5][int(sys.argv[1])-1]
num_segments = 20


# Name of mesh
mesh_name= f'residual/mesh/n{n}/voronoi{seed}.xdmf'
f_name = f'plots/residual/damping/damping{max_damping}'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)

generate_network(W,H,n,seed,char_length,num_segments,mesh_threshold, root_dir = 'residual')

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

crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['r','r'], abs_tol = 1e-10, rel_tol = 1e-12, init_c0 = init_c0, max_damping = max_damping)

np.savetxt(f'{f_name}/crit_strain.txt', np.array([crit_strain]))


crit_strains = np.loadtxt(f'phase_diagram/crit_strain/n{int(n)}.txt')
step_size = crit_strains[seed][1]*L/100.

print('Starting FEA run!')

fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], 1.5*crit_strains[seed][1], compute_individual=False, max_damping = max_damping, residual = True, init_c0 = init_c0, abs_tol = 1e-10, rel_tol = 1e-12) # make runs stop after strain stiffening 1.3*0.12541075747605646
print(f'max damping: {max_damping}')
step_size = crit_strain*L/100.
fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], 1.5*crit_strain, compute_individual=False, max_damping = max_damping, residual = True, init_c0 = init_c0, abs_tol = 1e-10, rel_tol = 1e-12) # make runs stop after strain stiffening