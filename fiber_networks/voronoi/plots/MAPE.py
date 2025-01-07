import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from matplotlib.ticker import PercentFormatter
from scipy.spatial.distance import jensenshannon
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
n_all = [100,200,300,400,500]
E = 1.0 # Stiffness modulus
num_seed = 20
L_0 = 10000
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

theta_N = []
strain = []
L_C_all = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []

error_stretch_bend_all = []
error_stretch_all = []

area_stretch_bend_all = []
area_stretch_all = []


for n in n_all:
    f_cont = f'shortest_paths/n{n}.csv'
    df=pd.read_csv(f_cont, index_col=0, header=0)
    for seed in range(num_seed):
        # Construct graphs
        nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
        edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
        G = graph.create_fiber_network(nodes,edges)

        edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
        edge_dist = np.array(list(edge_dist))
        L_char = np.mean(edge_dist)

        paths = df[f'seed {seed}'].dropna()
        num_paths=len(paths)
        total_path_dist=0
        paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])

        pos = nx.get_node_attributes(G,'pos')

        path_edges_idx = []
        for jj in range(len(paths)):
            path_edges_idx.append([])

        nx.set_edge_attributes(G,0,'shortest_path')

        for count,path in enumerate(paths):
            path_edges = list(zip(path,path[1:]))
            for edge in path_edges:
                G[edge[0]][edge[1]]['shortest_path'] = 1
                path_edges_idx[count].append(G[edge[0]][edge[1]]['index'])


        shortest_path_vec = np.array(list(nx.get_edge_attributes(G,'shortest_path').values()))
        for ii,kappa_tilde in enumerate(kappa_tildes):

            r = 2*L_char*np.sqrt(kappa_tilde)
            S = np.pi*r**2
            I = np.pi*r**4/4

            f_root = f'../fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}'
            strain = np.loadtxt(f'{f_root}/disp.txt')/L_0
            cont = np.loadtxt(f'{f_root}/cont.txt')
            bend = np.loadtxt(f'{f_root}/bend_energy.txt')
            stretch = np.loadtxt(f'{f_root}/stretch_energy.txt')
            shear = np.loadtxt(f'{f_root}/shear_energy.txt')
            force = np.loadtxt(f'{f_root}/force.txt')

            normalized_force = force/(E*S)

            stretching=0
            bending=0
            init_theta_all = []
            for edges in path_edges_idx:
                init_contour = np.sum(edge_dist[edges])

                L_C0 = init_contour/L_0
                # stretching component
                L_C = np.sum(cont[:,edges], axis=1)
                stretching += (L_C/init_contour)-1
                L = L_0*(1+strain)
                init_theta=np.arccos(L_C0**(-1))
                init_theta_all.append(init_theta)

                #bending component
                sin_theta = np.sqrt((L_C)**2-L**2)/(L_C)
                theta = np.arccos(L/L_C)
                bending += 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta


            error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))
            error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

            crit_strain = np.loadtxt(f'../phase_diagram/crit_strain/n{int(n)}.txt')[seed][ii]
            area_stretch = np.trapz(np.abs(normalized_force - stretching), strain/crit_strain)
            area_stretch_bend = np.trapz(np.abs(normalized_force - (stretching+bending)), strain/crit_strain)
            
            error_stretch_all.append([np.mean(np.array(init_theta_all)),L_0/L_char,error_stretch,ii])
            error_stretch_bend_all.append([np.mean(np.array(init_theta_all)),L_0/L_char,error_bend_stretch,ii])

            area_stretch_all.append([np.mean(np.array(init_theta_all)),L_0/L_char,area_stretch,ii])
            area_stretch_bend_all.append([np.mean(np.array(init_theta_all)),L_0/L_char,area_stretch_bend,ii])

print(len(area_stretch_all))

f_name = './MAPE_data/'
pathlib.Path(f_name).mkdir(parents=True, exist_ok=True)
np.savetxt(f'{f_name}error_stretch.txt',np.array(error_stretch_all))
np.savetxt(f'{f_name}error_stretch_bend.txt',np.array(error_stretch_bend_all))

np.savetxt(f'{f_name}area_stretch.txt',np.array(area_stretch_all))
np.savetxt(f'{f_name}area_stretch_bend.txt',np.array(area_stretch_bend_all))