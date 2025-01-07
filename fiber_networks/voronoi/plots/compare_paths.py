import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from matplotlib.ticker import PercentFormatter
import pandas as pd
#----------Import functions---------------------#
import sys
import os
sys.path.append('../../../utils')
import graph

plt.style.use('jeff_style.mplstyle')

# fiber network parameters
W = 1e4
H = 1e4
n = 300
E = 1.0
seed = 1
L_0 = 10000
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]
num_kappa_tilde = 4

# Construct graphs
nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

G_boundary = graph.add_boundary_nodes(G.copy(),0,H,W)

# Load shortest Paths
f_cont = f'shortest_paths/n{n}.csv'
df=pd.read_csv(f_cont, index_col=0, header=0)

paths = df[f'seed {seed}'].dropna()
num_paths=len(paths)

paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])


ii=0
pos = nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))

L_char = np.mean(edge_dist)

theta_N = []
strain = []
L_C_all = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []

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

path_edges_idx = []
for ii in range(len(paths)):
    path_edges_idx.append([])

for count,path in enumerate(paths):
    path_edges = list(zip(path,path[1:]))
    for edge in path_edges:
        path_edges_idx[count].append(G[edge[0]][edge[1]]['index'])

flattened_paths = sum(path_edges_idx,[])
support_network = np.arange(G.number_of_edges())
flattened_paths = sum(path_edges_idx,[])
support_network = np.delete(support_network, flattened_paths)

crit_strain = np.loadtxt(f'../phase_diagram/crit_strain/n{int(n)}.txt')[seed][:]

#----------Plot Energy Partition-----------#
# total energy
for ii in range(len(kappa_tildes)):
    energy_path_total_all = 0
    energy_path_stretch_all = 0
    energy_path_bend_all = 0

    # total energy plot
    fig1,ax1 = plt.subplots(figsize = (3,2.5))
    ax1.set_title('Total Energy')
    ax1.set_xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
    ax1.set_ylabel(r'Percent Energy')

    # stretching energy plot
    fig2,ax2 = plt.subplots(figsize = (3,2.5))
    ax2.set_title('Stretching Energy')
    ax2.set_xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
    ax2.set_ylabel(r'Percent Energy')

    # bending energy plot
    fig3,ax3 = plt.subplots(figsize = (3,2.5))
    ax3.set_title('Bending Energy')
    ax3.set_xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
    ax3.set_ylabel(r'Percent Energy')

    # force chains
    for count,edges in enumerate(path_edges_idx):
        energy_path_total = np.sum(total_energy[ii][:,edges],axis=1)/np.sum(total_energy[ii],axis=1) # total energy
        energy_path_stretch = np.sum(stretch_energy[ii][:,edges],axis=1)/np.sum(stretch_energy[ii],axis=1)
        energy_path_bend =  np.sum(bend_energy[ii][:,edges],axis=1)/np.sum(bend_energy[ii],axis=1)
        energy_path_total_all += energy_path_total
        energy_path_stretch_all += energy_path_stretch
        energy_path_bend_all += energy_path_bend
    ax1.plot(strain[ii]/crit_strain[ii],energy_path_total_all, lw = 3, color = 'k')
    ax2.plot(strain[ii]/crit_strain[ii],energy_path_stretch_all, lw = 3, color = 'k')
    ax3.plot(strain[ii]/crit_strain[ii],energy_path_bend_all, lw = 3, color = 'k')

    # support network
    energy_support_total = np.sum(total_energy[ii][:,support_network],axis=1)/np.sum(total_energy[ii],axis=1)
    energy_support_stretch = np.sum(stretch_energy[ii][:,support_network],axis=1)/np.sum(stretch_energy[ii],axis=1)
    energy_support_bend = np.sum(bend_energy[ii][:,support_network],axis=1)/np.sum(bend_energy[ii],axis=1)

    ax1.plot(strain[ii]/crit_strain[ii],energy_support_total, lw = 3, color = (0.9,0.5,0.5))
    ax1.set_xlim([0,1.5])
    ax1.set_ylim([0,1])
    fig1.tight_layout()
    fig1.savefig(f'total_energy_kappa{kappa_tildes[ii]}.pdf')

    ax2.plot(strain[ii]/crit_strain[ii],energy_support_stretch, lw = 3, color = (0.9,0.5,0.5))
    ax2.set_xlim([0,1.5])
    ax2.set_ylim([0,1])
    fig2.tight_layout()
    fig2.savefig(f'stretch_energy_kappa{kappa_tildes[ii]}.pdf')

    ax3.plot(strain[ii]/crit_strain[ii],energy_support_bend, lw = 3, color = (0.9,0.5,0.5))
    ax3.set_xlim([0,1.5])
    ax3.set_ylim([0,1])
    fig3.tight_layout()
    fig3.savefig(f'bend_energy_kappa{kappa_tildes[ii]}.pdf')
    

plt.show()