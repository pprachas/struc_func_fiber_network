import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
import pandas as pd

#----------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')

# shortest paths color scheme
def get_cmap(n, name='copper'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

L_0 = 1e4# fiber network parameters
n = 500
seed = 0
path_num = int(sys.argv[1])
num_edges = int(sys.argv[2])  # number of edges removed
edge_num = int(sys.argv[3]) # sample of edges

# Construct graphs
nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)
pos = nx.get_node_attributes(G,'pos')
n_nodes = max(list(G.nodes))

#----------Plot old and new paths-----------#

#-------Get shortest paths------------#
f_cont = f'../plots/shortest_paths/n{n}.csv'
df=pd.read_csv(f_cont, index_col=0, header=0)
paths = df[f'seed {seed}'].dropna()

paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])
paths = np.array(paths)

edges = np.loadtxt(f'sampled_edges/path{path_num}/num_edges{num_edges}.txt', dtype=int) 
edge_attr = np.array(list(nx.get_edge_attributes(G,'index').keys()))

# create ablated graph
G_ablated = G.copy()
edge = edge_attr[edges[edge_num]]
if num_edges == 1:
    G_ablated.remove_edge(*edge)
else:
    G_ablated.remove_edges_from(edge)
# compute ablated shortest path
G_boundary = graph.add_boundary_nodes(G_ablated.copy(),0,L_0,L_0)
paths_ablated = graph.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')
path_ablated_edges = list(zip(paths_ablated[path_num], paths_ablated[path_num][1:]))
path_edges = list(zip(paths[path_num],paths[path_num][1:]))

# plot fiber network
plt.figure(figsize=(3,3))
nx.draw_networkx_edges(G,pos, edge_color = 'k', style = ':', width = 0.75)
# plot og path
nx.draw_networkx_edges(G,pos, edge_color = 'k', edgelist = path_edges, width = 2.0, label='original path')
# plot new path
nx.draw_networkx_edges(G,pos, edge_color = (0.9,0.5,0.5), edgelist = path_ablated_edges, width = 2.0, label='new path')
# plot removed edge(s)
if num_edges == 1:
    nx.draw_networkx_edges(G,pos, width = 3.0, edgelist = [edge_attr[edges[edge_num]]], edge_color ='r', label = 'removed edge')
else:
    nx.draw_networkx_edges(G,pos, width = 3.0, edgelist = edge_attr[edges[edge_num]], edge_color ='r', label = 'removed edge')
plt.legend()

plt.savefig(f'n{n}_graph_edge{edge_num}.pdf')
plt.show()