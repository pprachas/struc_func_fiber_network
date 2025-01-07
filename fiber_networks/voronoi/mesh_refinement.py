import numpy as np
import pathlib
import networkx as nx
import matplotlib.pyplot as plt
from dolfin import *
#--------Import functions---------------------#
import sys
sys.path.append('..')
sys.path.append('../../utils')
from generate_voronoi import generate_network

import graph
import fea_HR as fea

crit_prev = 1e-16
diff = []
crit_strain_all = []
num_nodes = []

mesh_threshold= 0 # [0.0,0.1,0.5,1.0,1.5][arg] # meshing parameter (ignore shorter fibers)

# Chain parameters
L = 10000 # total length of fiber
W = L
H = L
n = 500 # voronoi seeds
char_length = L
seed = 0

for ii in range(3,20):
    num_segments = ii
    #pathlib.Path(f'./mesh_refinement/n{n}/seed{seed}/').mkdir(parents=True, exist_ok=True)
    mesh_name= f'mesh_refinement/damped/mesh/n{n}/voronoi{seed}.xdmf'

    kappa_tilde = 1e-6
    f_name = f'plots/n{n}/kappa{kappa_tilde}/seed{seed}/data'

    # run mesh and find critical point
    generate_network(W,H,n,seed,char_length,num_segments,mesh_threshold,f'mesh_refinement/damped')

    nodes = np.loadtxt(f'graph/nodes/n{n}/seed{seed}.txt')
    edges = np.loadtxt(f'graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
    G = graph.create_fiber_network(nodes,edges)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    L_char = np.mean(edge_dist) 
    r = 2*L_char*np.sqrt(kappa_tilde)

    print(L/L_char)
    step_size = L/1000.
    mesh,_,_ = fea.read_mesh(mesh_name)

    crit_strain, diff_energy = fea.run_critical_strain(mesh_name, r, L, step_size, ['r','r'], relax = 1.0, quad = 3, rel_tol = 1e-12, abs_tol = 1e-10, init_c0 = 1e-7, max_damping_res = 1e-4)

    diff.append(np.abs((crit_strain-crit_prev)/(crit_prev)))
    num_nodes.append(len(mesh.coordinates()))
    crit_strain_all.append(crit_strain)

    print('Percent difference of critical strain:', diff)
    print('Critical Strain:',crit_strain_all)
    print('Number of nodes',num_nodes)

    crit_prev = crit_strain

    plt.figure()
    plt.loglog(num_nodes,diff)

    if diff[-1] < 0.005:
        print('Final argument:', ii)
        np.savetxt(f'mesh_refinement/damped/mesh/n{n}/final_char{seed}.txt', [ii], fmt = '%i')
        np.savetxt(f'mesh_refinement/damped/mesh/n{n}/crit_strain{seed}.txt',np.array(crit_strain_all))
        print('Final argument:', ii)
        break
        #     fea.run_fea(mesh_name, f_name, r, L, step_size, ['r','r'], 1.01*crit_strain, compute_individual=False, relax = 0.99, quad = 3, c0=1e-8) # make runs stop after strain stiffening

        #     np.savetxt(f'mesh_refinement/damped/mesh/n{n}/num_nodes{seed}.txt', np.array(num_nodes))
        #     np.savetxt(f'mesh_refinement/damped/mesh/n{n}/diff{seed}.txt', np.array(diff))
        #     np.savetxt(f'mesh_refinement/damped/mesh/n{n}/final_char{seed}.txt', [ii], fmt = '%i')
        #     np.savetxt(f'mesh_refinement/damped/mesh/n{n}/crit_strain{seed}.txt',np.array(crit_strain_all))
        #     break
        # except RuntimeError:
        #     print(f'Force displacement curve did not converge at {ii} segments!')
        #     continue