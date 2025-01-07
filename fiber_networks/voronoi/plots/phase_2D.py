import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')

L_0 = 10000

# Load files for random chains
n_all = [100,200,300,400,500]
num_paths = 5
num_kappa_tilde = 4

kappa_tildes=[1e-3,1e-4,1e-5,1e-6]
L = 10000

num_seed = 20

init_cont_mean = []
crit_strain=[]

dist = []

first_contour = []
first_contour_all = []

l_c = []

num_paths_all = []
crit_strain_all = [] 

test = 0
for num_kappa_tilde in range(4):
    crit_strain = []
    first_contour = []
    test = 0
    for n in n_all:
        f_cont = f'shortest_paths/n{n}.csv'
        df=pd.read_csv(f_cont, index_col=0, header=0)
        for seed in range(num_seed):
            nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
            edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')

            # get graph
            G = graph.create_fiber_network(nodes,edges)
            edge_dist = np.array(list(nx.get_edge_attributes(G,'dist').values()))

            l_c.append(np.mean(edge_dist))

            paths = df[f'seed {seed}'].dropna()

            num_paths=len(paths)
            total_path_dist=0
            paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])
            for count,path in enumerate(paths):
                if count > 0:
                    break
                path_edges = list(zip(path,path[1:]))
                for edge in path_edges:
                    total_path_dist+=G[edge[0]][edge[1]]['dist']
                #print(total_path_dist)
            
            first_contour.extend([total_path_dist/L_0])
            num_paths_all.extend([num_paths])

        f_crit = f'../phase_diagram/crit_strain/n{n}.txt'
        crit_strain.extend(np.loadtxt(f_crit)[:,num_kappa_tilde])

    crit_strain_all.append(crit_strain)
    first_contour_all.append(first_contour)

first_contour_all=np.array(first_contour_all)-1
crit_strain_all=np.array(crit_strain_all)   

print(first_contour_all.shape)
print(crit_strain_all.shape)
l_c = np.array(l_c)    

pearson_corr = []

c = [(0.0,0.0,0.0), (0.5,0.5,0.5), (0.9,0.5,0.5),(0.8,0.2,0.2), (0.4,0.0,0.0)]
for ii in range(4):
    fig,ax = plt.subplots(figsize = (4,3))
    kappa_tilde = kappa_tildes[ii]
    # split for colors
    x = np.split(first_contour_all[ii],5)
    y = np.split(crit_strain_all[ii],5)

    for jj,n in enumerate(n_all):
        plt.plot(x[jj], y[jj], marker = 'o', fillstyle = 'none', ls = 'none', label = fr'$n = {n}$',c=c[jj])
        plt.xlabel(r'Initial Tortuosity $\tau_0$')
        plt.ylabel(r'Critical Strain $\varepsilon_{crit}$')
        plt.title(rf'$\tilde{{\kappa}} = 10^{{{np.log10(kappa_tilde):.0f}}}$')
    pearson_corr.append(np.corrcoef(first_contour_all[ii],crit_strain_all[ii])[0,1])
    #plt.legend()
    plt.tight_layout()
    plt.savefig(f'fiber_phase_2D_kappa{ii}.pdf')

#-------Pearson Correlation-------#
plt.figure(figsize=(8,2.5))
plt.title(r'Correlation between $\varepsilon_{crit}$ and $\tau_0$')
plt.semilogx(kappa_tildes, pearson_corr,marker = 'o', ls = 'none', c = 'k')
plt.xlabel(r'Dimensionless Bending Length $\tilde{\kappa}$')
plt.ylabel(r'Pearson Correlation')
plt.tight_layout()
plt.savefig('fiber_corr_kappa.pdf')

# Look at first seeds of sparse medium dense networks:
seeds = [100,200,300,400,500]


print(crit_strain_all.shape)


split_x = []
split_y = []

for ii in range(4):
    y = np.split(crit_strain_all[ii],5)
    split_y.append(y)
    split_x.append([kappa_tildes[ii]]*20)

split_x = np.array(split_x)
split_y = np.array(split_y)

print(split_x.shape, split_y.shape)
for ii in range(5):
    plt.figure(figsize=(3,3))
    plt.title(fr'$n = {seeds[ii]}$')
    plt.xlabel(r'Dimensionless Bending Length $\tilde{\kappa}$')
    plt.ylabel(r'Critical Strain $\varepsilon_{crit}$')
    plt.semilogx(split_x[:,0].ravel(), split_y[:,ii,0].ravel(), marker = 'o', ls = 'none', c = 'k') # only look at first seed for all n
    plt.tight_layout()
    plt.savefig(f'eps_kappa_n{seeds[ii]}.pdf')
plt.show()