import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
#----------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')

# fiber network parameters
W = 1e4
H = 1e4
n = int(sys.argv[1])
seed = int(sys.argv[2])
char_length = W/100.
n_paths = 1

# Construct graphs
nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)

#--------------Plot Graph-----------------------------#
fig,ax = plt.subplots(figsize=(7,7))
plt.axis('off')
pos = nx.get_node_attributes(G,'pos')
nx.draw_networkx_edges(G,pos, width = 1., edge_color = 'k')
# Note: the last two nodes are the sources and sinks respectively
plt.savefig(f'fiber_{n}_seed{seed}.pdf')

#--------------Shortest path demo---------------------#
fig,ax = plt.subplots(figsize=(7,7))

G_boundary = graph.add_boundary_nodes(G.copy(),0,H,W)
n_nodes = max(list(G.nodes))
path = nx.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')
pos = nx.get_node_attributes(G_boundary,'pos')
nx.draw_networkx_nodes(G_boundary,pos, node_size = 20, ax=ax, node_color = (0.9,0.5,0.5), nodelist = [len(G),len(G)+1] )
nx.draw_networkx_edges(G_boundary,pos, width = 1., edge_color = (0.9,0.5,0.5), style = ':')
nx.draw_networkx_nodes(G,pos, node_size = 10, ax=ax, node_color = 'r')
nx.draw_networkx_edges(G,pos, width = 1.0, edge_color = '0.5')
ii=0

path_edges = list(zip(path,path[1:]))
nx.draw_networkx_edges(G_boundary,pos, edgelist=path_edges, width = 1.5, edge_color = 'k')

plt.axis('off')
nodes_elements = [Line2D([0], [0], marker='o', color = 'r', ls = 'None', label='Crosslinks',
                          markerfacecolor='r', markersize=5), 
                Line2D([0], [0], marker='o', color = (0.9,0.5,0.5), ls = 'None', label='Boundary',
                        markerfacecolor=(0.9,0.5,0.5), markersize=5)]

edges_elements = [Line2D([0], [0], color = '0.5', label='Fibers'), 
                Line2D([0], [0], color = (0.9,0.5,0.5), ls = ':',label='Boundary'),
                Line2D([0], [0], color = 'k', linewidth = 4.0, label='Shortest Path')]


nodes = ax.legend(handles=nodes_elements, loc='lower right', title = 'Nodes')
ax.legend(handles=edges_elements, loc='lower left', title = 'Edges', bbox_to_anchor=(0.0, -0.03))
plt.gca().add_artist(nodes)
plt.savefig('shortest_path_demo.pdf')


G.remove_nodes_from(path)
G_boundary.remove_nodes_from(path[1:-1])
pos = nx.get_node_attributes(G_boundary,'pos')

fig,ax = plt.subplots(figsize=(7,7))
plt.axis('off')
nx.draw_networkx_nodes(G_boundary,pos, node_size = 20, ax=ax, node_color = (0.9,0.5,0.5), nodelist = [n_nodes+1,n_nodes+2] )
nx.draw_networkx_edges(G_boundary,pos, width = 1., edge_color = (0.9,0.5,0.5), style = ':')
nx.draw_networkx_nodes(G,pos, node_size = 20, ax=ax, node_color = 'r')
nx.draw_networkx_edges(G,pos, width = 1.0, edge_color = '0.5')

plt.savefig('ablated_path_demo.pdf')

plt.show()