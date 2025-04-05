import numpy as np
import matplotlib.pyplot as plt
import sys
import networkx as nx
import matplotlib as mpl
import pandas as pd

#--------Import functions---------------------#
sys.path.append('../../../utils')
import graph

# Load files for random chains
n_all = [100,300,500] #int(sys.argv[1])
seeds = [0,1,2]
L = 10000
E = 1.0
num_edges = [1]
kappa_tilde = 1e-6

# def get_cmap(n, name='Spectral'):
#             '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
#             RGB color; the keyword argument name must be a standard mpl colormap name.'''
#             return plt.cm.get_cmap(name, n)

cmap = np.array([(1.0,0.0,0.0),(0,0,0.6)])
# Plot settings
plt.style.use('jeff_style.mplstyle')

for n in n_all:
    fig,ax = plt.subplots(1,3, figsize=(7,2.5)) 
    #fig.suptitle(fr'n={n}')
    for count,seed in enumerate(seeds):
        # compute ES to normalized force
        nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
        edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
        G = graph.create_fiber_network(nodes,edges)

        pos=nx.get_node_attributes(G,'pos')
        edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
        edge_dist = np.array(list(edge_dist))
        init_contour= np.sum(edge_dist)
        L_char = np.mean(edge_dist)
        r = 2*L_char*np.sqrt(kappa_tilde)
        S = np.pi*r**2
        ES = E*S

        #-------------Load og fiber network----------------#
        f_root = f'../fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}'
        disp = np.loadtxt(f'{f_root}/disp.txt')/L
        force = np.loadtxt(f'{f_root}/force.txt')

        #-------Get shortest paths------------#
        f_cont = f'../plots/shortest_paths/n{n}.csv'
        df=pd.read_csv(f_cont, index_col=0, header=0)
        paths = df[f'seed {seed}'].dropna()

        paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])
        paths = np.array(paths)
        edge_idx = nx.get_edge_attributes(G,'index')
        edge_idx = edge_idx.keys()

        #-----------Load ablated fiber network-------------#
        for path_num in range(1):
            path = paths[path_num]
            path = list(zip(path,path[1:]))
            for num_edge in num_edges:
                relative_length = np.loadtxt(f'path_length/n{n}/seed{seed}/path{path_num}/num_edges{num_edge}.txt')
                c = relative_length ==1
                # norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
                # cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.bwr)
                # cmap.set_array([])
                for edge in range(len(path)):
                    f_root_ablated = f'fea_run/ablated/n{n}/path{path_num}/num_edge{num_edge}/edge{edge}/kappa_tilde{kappa_tilde}/seed{seed}'
                    disp_ablated = np.loadtxt(f'{f_root_ablated}/disp.txt')/L
                    force_ablated = np.loadtxt(f'{f_root_ablated}/force.txt')/ES

                    if (c==1).all(): 
                        ax.ravel()[count].plot(disp_ablated,force_ablated, lw = 1.25, c = cmap[1])
                    else:
                        if c[edge] == 0:
                            ax.ravel()[count].plot(disp_ablated,force_ablated, lw = 1.25, c = cmap[0], zorder=5)
                        else:
                            ax.ravel()[count].plot(disp_ablated,force_ablated, lw = 1, c = cmap[1])
        ax.ravel()[count].plot(disp,force, lw = 1.5, c = '0.6', label = 'original network', zorder=10)
        ax.ravel()[count].set_title(f'Seed {seed}')
        ax.ravel()[count].set_xlabel(r'$\varepsilon_{applied}$')
        ax.ravel()[count].set_ylabel(r'$\sfrac{\sigma_{yy}}{E}$')

        plt.tight_layout()
        plt.savefig(f'n{n}_abalted.pdf')
plt.show()