import numpy as np
import networkx as nx
import pandas as pd
import sys
import pathlib


#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph


pathlib.Path(f'shortest_paths').mkdir(parents=True, exist_ok=True)


n_all = [100,200,300,400,500]
num_seed = 20
L_0 = 10000


for n in n_all: 
    paths_all = []
    index= []
    for seed in range(num_seed):
        nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
        edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')

        G = graph.create_fiber_network(nodes,edges)

        G_boundary = graph.add_boundary_nodes(G.copy(),0,L_0,L_0)

        # Note: the last two nodes are the sources and sinks respectively

        n_nodes = max(list(G.nodes))
        paths = graph.shortest_path(G_boundary, n_nodes+1, n_nodes+2, 'dist')
        paths_all.append(paths)


        index.append(f'seed {seed}')

    column = [f'path {ii}' for ii in range(max([len(jj) for jj in paths_all]))]

    df = pd.DataFrame(paths_all, index=index, columns=column).T
  
    print(df)


    df.to_csv(f'shortest_paths/n{n}.csv')
