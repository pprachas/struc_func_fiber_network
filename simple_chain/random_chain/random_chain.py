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
import generate_simple_chain as sc

t0 = time.time()
num_links = int(sys.argv[1]) # number of links
L = 10000 # total length of fiber
num_segments = np.loadtxt('mesh_refinement/final_char.txt', dtype = int) # extra careful with mesh refinement

crit_strain_all = []
init_contour = []
run = int(sys.argv[2])
w = float(sys.argv[3])
num_kappa_tilde = int(sys.argv[4])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][num_kappa_tilde]

# Name of mesh
mesh_name= f'mesh/w{int(w)}/n{num_links}/random_chain{run}.xdmf'

# file names to write results
f_crit_name = f'phase_diagram/crit_strain/w{int(w)}/n{num_links}.txt'
f_cont_name = f'phase_diagram/init_cont/w{int(w)}/n{num_links}.txt'

# run mesh and find critical point
sc.create_random_simple_chain(num_links, L, w, num_segments,run)
link_points = np.loadtxt(f'link_points/w{int(w)}/n{num_links}/link_points{run}.txt')

G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour = np.sum(edge_dist)
L_char = np.mean(edge_dist) 
r = 2*L_char*np.sqrt(kappa_tilde)
print(r/L)
print(init_contour/L)
step_size = L/1000.
crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size,['p','p'], abs_tol = 1e-8, rel_tol = 1e-10)
print(crit_strain, init_contour)
# Save critical strain values
crit_array = np.loadtxt(f_crit_name)
crit_array[run][num_kappa_tilde] = crit_strain
np.savetxt(f_crit_name,crit_array)

# Save initial contour distance

cont_array = np.loadtxt(f_cont_name)
print(crit_array)
cont_array[run][num_kappa_tilde] = init_contour
np.savetxt(f_cont_name,cont_array)

print(f'time: {time.time()-t0}')