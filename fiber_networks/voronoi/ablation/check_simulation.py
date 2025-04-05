import os
import pandas as pd

n_all = [100,300,500]
seeds = [0,1,2]
num_edges = [1]
kappa_tilde = 1e-6
ii=0
jj=0
idx = []
f_not_converged = []
for seed in seeds:
    for n in n_all:
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
                    f_path = f'fea_run/ablated/n{n}/path{path_num}/num_edge{num_edge}/edge{edge_num}/kappa_tilde{kappa_tilde}/seed{seed}/bend_energy.txt'
                    jj+=1

                    if os.path.isfile(f_path):
                        pass
                    else:
                        ii+=1
                        idx.append(jj)
                        f_not_converged.append(f_path)


print(f'Total number of simulations: {jj}')
print(f'simulations not converged: {ii}')
print(f'Index of simulations not converges {idx}')
print('Files not converged')
for count,kk in enumerate(f_not_converged):
    print(f'{count}: {kk}')
