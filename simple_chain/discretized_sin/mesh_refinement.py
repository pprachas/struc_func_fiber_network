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

amplitude = L/10
lmbda =  L/40

kappa_tilde = 1e-3

num_points = int(3*L/(lmbda/2))
# Name of mesh
pathlib.Path('./mesh_refinement').mkdir(parents=True, exist_ok=True)
mesh_name =  f'mesh/a{int(amplitude)}/simple_chain_lmbda{int(lmbda)}.xdmf'  

crit_prev = 1e-16
diff = []
crit_strain_all = []


for ii in range(5,20):
    # run mesh and find critical point
    sc.create_discretized_sin_simple_chain(amplitude, lmbda, num_points, L, ii)
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

    crit_strain_all.append(crit_strain)

    diff.append(np.abs((crit_strain-crit_prev)/(crit_prev)))

    print('Percent difference of critical strain:', diff)
    print('Critical Strain:',crit_strain_all)

    crit_prev = crit_strain

    plt.plot(diff)
    plt.savefig('mesh_refinement')

    if diff[-1] < 5e-3:
        print('Final argument:', ii)
        np.savetxt('mesh_refinement/final_char.txt', [ii], fmt = '%i')
        break
