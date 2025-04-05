import numpy as np 
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import gmsh
import meshio
import sys
import pathlib

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph
import meshing

n = 500 # density of network
seed = 0 # seed of initial Voronoi points

num_edges = [1,5,10,20]

rng = np.random.default_rng(seed=42)

nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')

# get graph
G = graph.create_fiber_network(nodes,edges)
edge_dist = np.array(list(nx.get_edge_attributes(G,'dist').values()))

#-------Get shortest paths------------#
f_cont = f'../plots/shortest_paths/n{n}.csv'
df=pd.read_csv(f_cont, index_col=0, header=0)
paths = df[f'seed {seed}'].dropna()

paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])

path_edges_idx = []

for ii in range(len(paths)):
    path_edges_idx.append([])

for count,path in enumerate(paths):
    path_edges = list(zip(path,path[1:]))
    for edge in path_edges:
        path_edges_idx[count].append(G[edge[0]][edge[1]]['index'])

# edges = [path_edges_idx[path_num][3]]

og_mesh = f'../mesh_msh/n{n}/voronoi{seed}.msh'

for path_num in range(2):
    pathlib.Path(f'sampled_edges/path{path_num}').mkdir(parents=True, exist_ok=True)
    for ii in num_edges:
        if ii==1:
            path = np.array(path_edges_idx[path_num])
            shuffled_path = rng.permutation(path)
            np.savetxt(f'sampled_edges/path{path_num}/num_edges{ii}.txt', shuffled_path[:20], fmt = '%i')

        else:
            path = np.array([path_edges_idx[path_num]]*20) #np.array([path_edges_idx[path_num]]*20)
            shuffled_path = rng.permuted(path, axis = 1)

            np.savetxt(f'sampled_edges/path{path_num}/num_edges{ii}.txt', shuffled_path[:,:ii], fmt = '%i')
        