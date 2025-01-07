import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation

from scipy.stats import entropy


#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
import graph


#---------------Random Chain-----------------------------------#
plt.style.use('jeff_style.mplstyle')

L_0 = 10000
kappa_tilde = 1e-6 #[1e-3,1e-4,1e-5,1e-6]

E = 1.0 # Stiffness moduli

num_link = 80 #amplitude 
w = L_0/10 #[L_0/10,L_0/20,L_0/40] #amplitude
num_run=1

print()

link_points = np.loadtxt(f'../random_chain/link_points/w{int(w)}/n{int(num_link)}/link_points{num_run}.txt')
G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
orientation = np.fromiter(nx.get_edge_attributes(G,'orientation').values(), dtype=float)

edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour = np.sum(edge_dist)

init_straightness = init_contour/L_0
init_theta = np.arccos(init_straightness**(-1))

print(init_straightness)

plt.figure(figsize=(2,5))
edges = nx.draw_networkx_edges(G,pos, width = 2)
plt.axis('equal')

plt.savefig('random_chain_high.pdf')

plt.show()

