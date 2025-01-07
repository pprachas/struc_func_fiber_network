import numpy as np
import pathlib
import networkx as nx
import matplotlib.pyplot as plt

#--------Import functions---------------------#
import sys
sys.path.append('..')
sys.path.append('../../utils')
import generate_simple_chain as sc

import fea_HR as fea
import graph

crit_prev = 1e-6
diff = []
crit_strain_all = []

pathlib.Path('./mesh_refinement').mkdir(parents=True, exist_ok=True)

for ii in range(3,20):
    # parameters of chain
    num_links = 80 # number of links
    L = 10000 # total length of fiber
    w =  int(L/10) # maximum possible width
    seed = 0
    mesh_name= f'mesh/w{w}/n{num_links}/random_chain{seed}.xdmf'

    kappa_tilde = 1e-6

    # run mesh and find critical point
    sc.create_random_simple_chain(num_links, L, w, ii,seed)

    link_points = np.loadtxt(f'link_points/w{w}/n{num_links}/link_points{seed}.txt')
    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    init_contour = np.sum(edge_dist)
    L_char = np.mean(edge_dist) 
    r = 2*L_char*np.sqrt(kappa_tilde)
    step_size = L/100
    crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['p','p'], abs_tol = 1e-10, rel_tol = 1e-12)

    crit_strain_all.append(crit_strain)


    diff.append(np.abs((crit_strain-crit_prev)/(crit_prev)))

    print('Percent difference of critical strain:', diff)
    print('Critical Strain:',crit_strain_all)

    crit_prev = crit_strain

    plt.figure()
    plt.semilogy(diff)
    plt.savefig('mesh_refinement')

    if diff[-1] < 5e-3:
        print('Final argument:', ii)
        np.savetxt('mesh_refinement/final_char.txt', [ii], fmt = '%i')
        np.savetxt('mesh_refinement/crit_strain.txt', crit_strain_all)
        break
