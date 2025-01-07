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
lmbda =  (L/np.array([5,10,20,40]))[run]

crit_strain_all = []
init_contour = []
num_kappa_tilde = int(sys.argv[3])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][num_kappa_tilde]


# Name of mesh
mesh_name= f'mesh/a{int(amplitude)}/simple_chain_lmbda{int(lmbda)}.xdmf'

# file names to save results
f_name = f'fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'

pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)


# run mesh and find critical point
link_points = np.loadtxt(f'link_points/a{int(amplitude)}/lmbda{int(lmbda)}.txt')

G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour= np.sum(edge_dist)
L_char = (lmbda/(2*L))*init_contour
r = 2*L_char*np.sqrt(kappa_tilde)


crit_strains = np.loadtxt(f'phase_diagram/crit_strain/a{int(amplitude)}.txt')
step_size = crit_strains[run][num_kappa_tilde]*L/50.
fea.run_fea(mesh_name, f_name, r, L, step_size, ['p','p'], crit_strains[run][num_kappa_tilde]*1.5, compute_individual=True, abs_tol = 1e-8, rel_tol = 1e-10) # make runs stop after strain stiffening