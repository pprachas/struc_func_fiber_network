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
import generate_simple_chain as sc

L = 10000

amplitude = float(sys.argv[1])
run = int(sys.argv[2])
lmbda =  (L/np.array([1,5,10,20,40]))[run]
num_kappa_tilde = int(sys.argv[3])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][num_kappa_tilde]
num_segments = int(np.loadtxt('mesh_refinement/final_char.txt'))*2

print(num_segments)
num_points = int(3*L/(lmbda/2))
# Name of mesh
mesh_name =  f'mesh/a{int(amplitude)}/simple_chain_lmbda{int(lmbda)}.xdmf'  

# file names to save results
f_crit_name = f'phase_diagram/crit_strain/a{int(amplitude)}.txt'
f_cont_name = f'phase_diagram/init_cont/a{int(amplitude)}.txt'

# run mesh and find critical point
sc.create_discretized_sin_simple_chain(amplitude, lmbda, num_points, L, num_segments)
link_points = np.loadtxt(f'link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')

G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)
step_size = L/100.

crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['p','p'], abs_tol = 1e-8, rel_tol = 1e-10)

crit_array = np.loadtxt(f_crit_name)
crit_array[run][num_kappa_tilde] = crit_strain

print(crit_array.shape)
# Save critical strain values
crit_array = np.loadtxt(f_crit_name)
crit_array[run][num_kappa_tilde] = crit_strain
np.savetxt(f_crit_name,crit_array)

# Save initial contour distance

cont_array = np.loadtxt(f_cont_name)
print(crit_array)
cont_array[run][num_kappa_tilde] = init_contour
np.savetxt(f_cont_name,cont_array)