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

# fiber network parameters
W = 1e4
H = 1e4
n = int(sys.argv[1])
seed = int(sys.argv[2])-1

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

nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = 'k', style = ':')
ii=0
for count,path in enumerate(paths):
    ii+=1
    path_edges = list(zip(path,path[1:]))
    nx.draw_networkx_edges(G_boundary,pos, edgelist=path_edges, width = 4, edge_color=cmap(count), label = f'{count+1}')
plt.title('Shortest Paths')
plt.legend()

#-------------Plotting energy distribution at onset of stiffening----------#
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][int(sys.argv[3])-1]
# stretching_energy = np.loadtxt(f'./plots/bisection/roller-roller/n{n}/kappa{kappa_tilde}/seed{seed}/data/stretch_energy.txt')
# bending_energy = np.loadtxt(f'./plots/bisection/roller-roller/n{n}/kappa{kappa_tilde}/seed{seed}/data/bend_energy.txt')
# shear_energy = np.loadtxt(f'./plots/bisection/roller-roller/n{n}/kappa{kappa_tilde}/seed{seed}/data/shear_energy.txt')

# total_energy = stretching_energy+bending_energy+shear_energy

# percent_stretch = stretching_energy/total_energy

# print(stretching_energy.shape, G.number_of_edges())

# fig,ax = plt.subplots(figsize=(9,7))
# nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
# edges = nx.draw_networkx_edges(G,pos, width = 3.0, edge_color = percent_stretch, edge_cmap=plt.cm.Reds)
# nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = 'k', style = ':')
# plt.colorbar(edges)
# plt.title('Percent stretching energy distribution')

# plt.figure()
# sns.histplot(data = percent_stretch, bins = 50, stat = 'probability')
# plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
# plt.title('Histogram of percent stretching energy')

#----------plot time evolution---------------------#
stretching_energy = np.loadtxt(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/data/stretch_energy.txt')
bending_energy = np.loadtxt(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/data/bend_energy.txt')
shear_energy = np.loadtxt(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/data/shear_energy.txt')

total_energy = stretching_energy+bending_energy+shear_energy


pathlib.Path(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots').mkdir(parents=True, exist_ok=True)

print(stretching_energy.shape)
for ii in range(len(stretching_energy)): 
    plt.figure()
    percent_stretch = stretching_energy[ii]/total_energy[ii]

    nx.draw_networkx_nodes(G,pos, node_size = 15, ax=ax, node_color = 'r')
    edges = nx.draw_networkx_edges(G,pos, width = 3.0, edge_color = percent_stretch, edge_cmap=plt.cm.Reds, edge_vmin = 0.0, edge_vmax=1.0)
    nx.draw_networkx_edges(G,pos, width = 1.5, edge_color = 'k', style = ':')
    plt.colorbar(edges)
    plt.savefig(f'./plots/n{n}/kappa{kappa_tilde}/seed{seed}/plots/{ii:03}.png')
    plt.close()

os.system(f'convert -delay 1 -loop 0 plots/n200/kappa{kappa_tilde}/seed0/plots/*.png plots/n200/kappa{kappa_tilde}/seed0/plots/animated.gif')