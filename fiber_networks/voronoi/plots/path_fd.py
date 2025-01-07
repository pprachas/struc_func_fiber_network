import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from scipy.integrate import simpson
from matplotlib.ticker import PercentFormatter
from scipy.spatial.distance import jensenshannon
#----------Import functions---------------------#
import sys
import os
sys.path.append('../../../utils')
import graph

plt.style.use('jeff_style.mplstyle')

# fiber network parameters
W = 1e4
H = 1e4
n = 500 # n= 100,200,300,400,500
E = 1.0 # Stiffness modulus
seed = 1
L_0 = 10000
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]
num_kappa_tilde = 1

# Construct graphs
nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))

G_boundary = graph.add_boundary_nodes(G.copy(),0,H,W)

# Note: the last two nodes are the sources and sinks respectively

n_nodes = max(list(G.nodes))
paths = graph.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')

pos = nx.get_node_attributes(G,'pos')

path_edges_idx = []
for ii in range(len(paths)):
    path_edges_idx.append([])

nx.set_edge_attributes(G,0,'shortest_path')

for count,path in enumerate(paths):
    path_edges = list(zip(path,path[1:]))
    for edge in path_edges:
        G[edge[0]][edge[1]]['shortest_path'] = 1
        path_edges_idx[count].append(G[edge[0]][edge[1]]['index'])

shortest_path_vec = np.array(list(nx.get_edge_attributes(G,'shortest_path').values()))

theta_N = []
strain = []
L_C_all = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []

L_char = np.mean(edge_dist)

crit_strain = np.loadtxt(f'../phase_diagram/crit_strain/n{int(n)}.txt')[seed][:]

for kappa_tilde in kappa_tildes:
    f_root = f'../fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}'

    disp = np.loadtxt(f'{f_root}/disp.txt')
    cont = np.loadtxt(f'{f_root}/cont.txt')
    bend = np.loadtxt(f'{f_root}/bend_energy.txt')
    stretch = np.loadtxt(f'{f_root}/stretch_energy.txt')
    shear = np.loadtxt(f'{f_root}/shear_energy.txt')
    force = np.loadtxt(f'{f_root}/force.txt')

    orientation = np.loadtxt(f'{f_root}/rot.txt')

    

    r = 2*L_char*np.sqrt(kappa_tilde)
    S = np.pi*r**2
    I = np.pi*r**4/4

    L_C_all.append(cont)
    theta_N.append(orientation)
    strain.append(disp/L_0)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)
    normalized_force.append(force/(E*S))

plt.figure(figsize=(3,3))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'Normalized Applied Strain $\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'Normalized Force $\displaystyle{\sfrac{\sigma_{yy}}{E}}$')
plt.title('Dense Fiber Network')

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

for ii in range(len(kappa_tildes)):
    stretching = 0
    bending = 0
    for edges in path_edges_idx:
        init_contour = np.sum(edge_dist[edges])

        L_C0 = init_contour/L_0
        # stretching component
        L_C = np.sum(L_C_all[ii][:,edges], axis=1)
        stretching += (L_C/init_contour)-1
        L = L_0*(1+strain[ii])
        init_theta = np.arccos(L_C0**(-1))

        #bending component
        sin_theta = np.sqrt((L_C)**2-L**2)/(L_C)
        theta = np.arccos(L/L_C)
        bending += 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta
    plt.loglog(strain[ii]/crit_strain[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 1.75, c = color[ii])
    plt.loglog(strain[ii]/crit_strain[ii],stretching, lw = 0.5, c = color[ii], marker = 'o', fillstyle='none', markersize=5, markeredgewidth = 0.5)
    plt.loglog(strain[ii]/crit_strain[ii],stretching+bending, lw = 0.5, c = color[ii], marker = '*', fillstyle='none', markersize=5, markeredgewidth = 0.5)

plt.tight_layout()
plt.savefig(f'fiber_rom{n}.pdf')

plt.figure()
plt.plot(strain[ii]/crit_strain[ii],stretching, lw = 0.5, c = color[ii], marker = '*', fillstyle='none', markersize=5, markeredgewidth = 0.5)
plt.plot(strain[ii]/crit_strain[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 1.75, c = color[ii])

print(normalized_force[ii])
# print(normalized_force[ii][-1], stretching[ii])
# print((normalized_force[ii] - stretching[ii]))
print(np.mean(np.abs((normalized_force[ii] - stretching)/normalized_force[ii])))
# print(np.mean(np.abs((normalized_force[ii] - (stretching+bending))/normalized_force[ii])))

print(f_root)
# print(np.trapz(normalized_force[ii] - stretching, strain[ii]/crit_strain[ii]))
# print(np.trapz(normalized_force[ii] - (stretching+bending), strain[ii]/crit_strain[ii]))
#plt.legend(ncol=2, fontsize="12" )
plt.show()