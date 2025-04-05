import numpy as np
import networkx as nx
import pandas as pd
import sys
import pathlib
import matplotlib.pyplot as plt
import pathlib

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

n_all = [100,300,500]
L_0 = 10000
seeds = [0,1,2]
num_edges = [1]

for n in n_all:
    for seed in seeds:
        nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
        edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
        G = graph.create_fiber_network(nodes,edges)

        #-------Get shortest paths------------#
        f_cont = f'../plots/shortest_paths/n{n}.csv'
        df=pd.read_csv(f_cont, index_col=0, header=0)
        paths = df[f'seed {seed}'].dropna()

        paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])
        paths = np.array(paths)
        edge_idx = nx.get_edge_attributes(G,'index')
        edge_idx = edge_idx.keys()

        # load sampled edges 
        for ii in num_edges:
            for path_num in range(1):
                ablated_path_dist=[]
                path = paths[path_num]
                path = list(zip(path,path[1:]))
                for edge in range(len(path)):
                    G_ablated = G.copy()
                    n_nodes = max(list(G.nodes))
                    edge=path[edge]
                    if ii == 1:
                        G_ablated.remove_edge(*edge)
                    else:
                        G_ablated.remove_edges_from(edge)

                    G_boundary = graph.add_boundary_nodes(G_ablated.copy(),0,L_0,L_0)
                    paths_ablated = graph.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')
                    # convert node sequence to edges
                    path_ablated_edges = list(zip(paths_ablated[path_num], paths_ablated[path_num][1:]))

                    path_ablated_idx = []
                    path_idx = []
                    test_ablated = 0
                    test = 0
                    #-------------------------------------------------#
                    for edge in path_ablated_edges:
                            path_ablated_idx.append(G_ablated[edge[0]][edge[1]]['index'])
                            test_ablated += G_ablated[edge[0]][edge[1]]['dist']
                            
                    next_path = list(zip(paths[path_num+1], paths[path_num+1][1:]))
                    for edge in next_path:
                            path_idx.append(G[edge[0]][edge[1]]['index'])
                            test += G[edge[0]][edge[1]]['dist']

                    edge_dist_ablated = np.array(list(nx.get_edge_attributes(G,'dist').values()))
                    edge_dist = np.array(list(nx.get_edge_attributes(G,'dist').values()))
                    ablated_path_dist.append(np.sum(edge_dist_ablated[path_ablated_idx])/np.sum(edge_dist[path_idx]))
                ablated_path_dist = np.array(ablated_path_dist)

                # Create folder
                pathlib.Path(f'path_length/n{n}/seed{seed}/path{path_num}/').mkdir(parents=True, exist_ok=True)


                np.savetxt(f'path_length/n{n}/seed{seed}/path{path_num}/num_edges{ii}.txt', ablated_path_dist)

plt.show()