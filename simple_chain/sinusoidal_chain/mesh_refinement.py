import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx

#--------Import functions---------------------#
import sys
sys.path.append('..')
sys.path.append('../../utils')
import generate_simple_chain as sc

import fea_HR as fea
import graph

crit_prev = 1e-16
contour_prev = 1e-16
diff = []
diff_contour = []

pathlib.Path('./mesh_refinement').mkdir(parents=True, exist_ok=True)

contour_all = []
crit_strain_all = []

for ii in range(2,20):
    # parameters of chain (large amplitude small wavelength for extreme case)
    L = 10000 # total length of fiber
    lmbda = L/40 # wavelength 
    amplitude = L/10 #amplitude
    num_points = int(2*ii*(L/(lmbda/2))) 

    sc.create_sinusoidal_simple_chain(amplitude, lmbda, num_points, L)

    mesh_name =  f'mesh/a{int(amplitude)}/simple_chain_lmbda{int(lmbda)}.xdmf'  

    kappa_tilde = 1e-6

    link_points = np.loadtxt(f'link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')

    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    init_contour = np.sum(edge_dist)
    L_char = (lmbda/(2*L))*init_contour
    r = 2*L_char*np.sqrt(kappa_tilde)
    step_size = L/100

    crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size,['p','p'], abs_tol = 1e-8, rel_tol = 1e-10)

    crit_strain_all.append(crit_strain)

    diff.append(np.abs((crit_strain-crit_prev)/crit_prev))
    diff_contour.append(np.abs((init_contour-contour_prev)/contour_prev))
    contour_all.append(init_contour)

    print('Critical Strain:', crit_strain_all)
    print('Percent difference of critical strain:', diff)

    print('Initial Contour:', contour_all)
    print('Percent Difference of contour length:', diff_contour)

    crit_prev = crit_strain
    contour_prev = init_contour

    plt.figure()
    plt.semilogy(diff)
    plt.savefig('mesh_refinement')

    if diff[-1] < 5e-3 and diff_contour[-1] < 5e-3:
        print('Final argument:', ii)
        np.savetxt('mesh_refinement/crit_strain.txt', crit_strain_all)
        np.savetxt('mesh_refinement/final_char.txt', [ii], fmt = '%i')
        break