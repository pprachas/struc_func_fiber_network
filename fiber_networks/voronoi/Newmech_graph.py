from generate_voronoi import generate_network
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from matplotlib.ticker import PercentFormatter
import seaborn as sns
#----------Import functions---------------------#
import sys
import os
sys.path.append('../../utils')
import graph

plt.style.use('jeff_style.mplstyle')

# fiber network parameters
W = 1e4
H = 1e4
n = 100
seed = 0
L_0 = 10000
kappa_tilde = 1e-6

# Construct graphs
nodes = np.loadtxt(f'graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

G_boundary = graph.add_boundary_nodes(G.copy(),0,H,W)

# Note: the last two nodes are the sources and sinks respectively

n_nodes = max(list(G.nodes))
paths = graph.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')

fig,ax = plt.subplots(figsize=(7,7))

ii=0
pos = nx.get_node_attributes(G,'pos')


#---------Plotting shortest paths-----------------------#
# function for plotting edge color
def get_cmap(n, name='copper'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

cmap = get_cmap(len(paths))

plt.figure(figsize=(7,7))
ii=0
for count,path in enumerate(paths):
    ii+=1
    path_edges = list(zip(path,path[1:]))
    nx.draw_networkx_edges(G_boundary,pos, edgelist=path_edges, width = 4, edge_color=cmap(count), label = f'{count+1}')
nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = '0.6', style =  ':')
plt.title('Shortest Paths Between Constrained Boundaries')
plt.legend(title = 'Shortest path number', ncol = 6, loc = 'lower center')
plt.tight_layout()
plt.savefig('shortest_path.pdf')

#-------------Plotting energy distribution at onset of stiffening----------#
#----------plot time evolution---------------------#
# stretching_energy = np.loadtxt(f'fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}/stretch_energy.txt')
# bending_energy = np.loadtxt(f'fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}/bend_energy.txt')
# shear_energy = np.loadtxt(f'fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}/shear_energy.txt')

# total_energy = stretching_energy+bending_energy+shear_energy


# pathlib.Path(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots').mkdir(parents=True, exist_ok=True)

# print(stretching_energy.shape)
# jj = 0
# for ii in range(len(stretching_energy)): 
#     plt.figure(figsize = (7,7))
#     if ii%5 == 0: 
#         plt.title('Evolution of Stretching Energy Distribution')
#         percent_stretch = stretching_energy[ii]/total_energy[ii]

#         nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
#         edges = nx.draw_networkx_edges(G,pos, width = 3.0, edge_color = percent_stretch, edge_cmap=plt.cm.Reds, edge_vmin = 0.0, edge_vmax=1.0)
#         nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = 'k', style = ':')
#         plt.colorbar(edges)
#         plt.savefig(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots/{jj:03}.png')
#         plt.close()
#         jj+=1

# os.system(f'convert -loop 0 plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots/*.png plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots/animated.gif')

# # #------------Plotting onset of transition------------------#
# fiber_strain = np.loadtxt(f'fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}/disp.txt')/L_0
# crit_strain = np.loadtxt(f'phase_diagram/crit_strain/n{int(n)}.txt')[seed][3]

# idx = np.argmin(np.abs(fiber_strain-crit_strain))

# print(idx)

# plt.figure(figsize = (7,7))
# plt.title('Force Chains')

# percent_stretch = stretching_energy[-1]/total_energy[-1]

# nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
# edges = nx.draw_networkx_edges(G,pos, width = 3.0, edge_color = percent_stretch, edge_cmap=plt.cm.Reds, edge_vmin = 0.0, edge_vmax=1.0)
# nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = 'k', style = ':')
# cbar=plt.colorbar(edges)
# cbar.set_label('Percent Stretching Energy in Fiber')
# plt.savefig(f'onset.pdf')
