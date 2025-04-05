import os
import pathlib
import numpy as np
import sys
import pandas as pd 

num_seed = [100,300,500] # number of links
num_edges = [1]
seeds = [0,1,2]

L = 10000 # total length of fiber
num_kappa_tilde = 4

combined_list = []

for seed in seeds:
    for n in num_seed:
        #-------Get shortest paths------------#
        f_cont = f'../plots/shortest_paths/n{n}.csv'
        df=pd.read_csv(f_cont, index_col=0, header=0)
        paths = df[f'seed {seed}'].dropna()
        paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])

        path = paths[0]
        path = list(zip(path,path[1:]))
        for path_num in range(1):
            for num_edge in num_edges:
                for edge_num in range(len(path)):
                    combined_list.append([n, seed, 3, path_num, num_edge, edge_num])

combined_list = [combined_list[ii-1] for ii in [170, 181]]
args = combined_list[int(sys.argv[1])-1]

print(f'python3 fea_ablated_voronoi.py {args[0]} {args[1]} {args[2]} {args[3]} {args[4]} {args[5]}', flush=True)

os.system(f'python3 fea_ablated_voronoi.py {args[0]} {args[1]} {args[2]} {args[3]} {args[4]} {args[5]} ')