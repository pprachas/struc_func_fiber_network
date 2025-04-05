import numpy as np 
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import gmsh
import meshio

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph
import meshing

n_all = [100,300,500]
seeds = [0,1,2]

num_edges = [1]

jj=0
#-------Get shortest paths------------#
for seed in seeds:
    for n in n_all:

        nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
        edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
        G = graph.create_fiber_network(nodes,edges)

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

        og_mesh = f'../mesh_msh/n{n}/voronoi{seed}.msh'

        for path_num in range(1): 
            print(path_edges_idx[path_num])
            for num_edge in num_edges:
                for edge_num, edge in enumerate(path_edges_idx[path_num]):
                    f_name = f'voronoi{seed}_edge{edge_num}'
                    #print(edge, edge_num)
                    jj+=1
                    meshing.create_ablated_mesh(f'n{n}/path{path_num}/num_edge{num_edge}',f_name,og_mesh,n,seed,[edge])

print(f'Number of simulations {jj}')