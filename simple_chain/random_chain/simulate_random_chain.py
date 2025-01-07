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

num_links = int(sys.argv[1]) # number of links
run = int(sys.argv[2])
w = float(sys.argv[3])
num_kappa_tilde = int(sys.argv[4])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][num_kappa_tilde]

# Name of mesh
mesh_name = f'mesh/w{int(w)}/n{num_links}/random_chain{run}.xdmf'

# file names to save results
f_name = f'fea_run/w{int(w)}/n{num_links}/kappa_tilde{kappa_tilde}/seed{int(run)}'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)


# run mesh and find critical point
link_points = np.loadtxt(f'link_points/w{int(w)}/n{num_links}/link_points{run}.txt')

G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)

crit_strains = np.loadtxt(f'phase_diagram/crit_strain/w{int(w)}/n{num_links}.txt')
step_size = crit_strains[run][num_kappa_tilde]*L/50.

fea.run_fea(mesh_name, f_name, r, L, step_size, ['p','p'], crit_strains[run][num_kappa_tilde]*1.5, compute_individual = True, abs_tol = 1e-8, rel_tol = 1e-10) # make runs stop after strain stiffening