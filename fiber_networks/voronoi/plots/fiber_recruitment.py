import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from scipy.interpolate import interp1d

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')

L_0 = 10000
E = 1.0

# Load files for random chains
n = int(sys.argv[1])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][int(sys.argv[3])]

seed = int(sys.argv[2])

f_cont = f'shortest_paths/n{n}.csv'
df=pd.read_csv(f_cont, index_col=0, header=0)

nodes = np.loadtxt(f'../graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'../graph/edges/n{n}/seed{seed}.txt', dtype = 'int')

# get graph
G = graph.create_fiber_network(nodes,edges)
edge_dist = np.array(list(nx.get_edge_attributes(G,'dist').values()))

index = np.array(list(nx.get_edge_attributes(G,'index').values()))
l_c= np.mean(edge_dist)


r = 2*l_c*np.sqrt(kappa_tilde)
A=np.pi*r**2

paths = df[f'seed {seed}'].dropna()
num_paths=len(paths)

path_edges_idx = []
contour_length_all = []
shortest_path_all = []

for ii in range(len(paths)):
    path_edges_idx.append([])

total_path_dist=0
paths = paths.apply(lambda x: [int(el) for el in x.strip("[]").split(",")])

for count,path in enumerate(paths):
    path_edges = list(zip(path,path[1:]))
    contour_length  = 0
    for edge in path_edges:
        path_edges_idx[count].append(G[edge[0]][edge[1]]['index'])
        contour_length += G[edge[0]][edge[1]]['dist']
    contour_length_all.append(contour_length)

shortest_path_all.append(path_edges_idx)

mean_contour=total_path_dist/num_paths


f_root = f'../fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}'
disp = np.loadtxt(f'{f_root}/disp.txt')
stretch = np.loadtxt(f'{f_root}/stretch_energy.txt')
bend = np.loadtxt(f'{f_root}/bend_energy.txt')
shear = np.loadtxt(f'{f_root}/shear_energy.txt')
force = np.loadtxt(f'{f_root}/force.txt')
cont = np.loadtxt(f'{f_root}/cont.txt')

strain=disp/L_0


f_crit = f'../phase_diagram/crit_strain/n{n}.txt'
crit_strain=np.loadtxt(f_crit)[seed][int(sys.argv[3])]

print(crit_strain)
mean_contour = np.array(mean_contour)/L_0-1
crit_strain=np.array(crit_strain)   
l_c = np.array(l_c)    
A = np.array(A)
EA=E*np.array(A)
contour_length_all=np.array(contour_length_all)/L_0

# #-------Get results from single chains------#

# Load files for random chains
w = [L_0/5,L_0/10,L_0/20,L_0/40,L_0/50,L_0/100]
n = [30,40,50,60,70]
w_r = [L_0/50,L_0/40,L_0/20,L_0/10]
crit_strain_random_all = []
init_cont_random_all = []
num_kappa_tilde = 4

# Load files for sinusoial chains
a = [250,500,1000] #[250,500,1000]
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

crit_strain_sin = []
init_cont_sin = []

crit_strain_tri = []
init_cont_tri = []

crit_strain_rand = []
init_cont_rand = []

kappa_tilde = []
kappa_tilde_rand = []
# Load all sinusoidal and triangular data
for ii in a:
    # Sinusoidal
    crit_strain_sin.append(np.loadtxt(f'../../../simple_chain/sinusoidal_chain/phase_diagram/crit_strain/a{ii}.txt')[:,int(sys.argv[3])])
    init_cont_sin.append(np.loadtxt(f'../../../simple_chain/sinusoidal_chain/phase_diagram/init_cont/a{ii}.txt')[:,int(sys.argv[3])]/L_0)

    # triangular
    crit_strain_tri.append(np.loadtxt(f'../../../simple_chain/triangular_chain/phase_diagram/crit_strain/a{ii}.txt')[:,int(sys.argv[3])])
    init_cont_tri.append(np.loadtxt(f'../../../simple_chain/triangular_chain/phase_diagram/init_cont/a{ii}.txt')[:,int(sys.argv[3])]/L_0)

# Load all random data
for ii in w_r:
    for jj in n:
        ii = int(ii)
        jj = int(jj)
        crit_strain_rand.append(np.loadtxt(f'../../../simple_chain/random_chain/phase_diagram/crit_strain/w{ii}/n{jj}.txt')[:,int(sys.argv[3])])
        init_cont_rand.append(np.loadtxt(f'../../../simple_chain/random_chain/phase_diagram/init_cont/w{ii}/n{jj}.txt')[:,int(sys.argv[3])]/L_0)

crit_strain_sin = np.array(crit_strain_sin).ravel()
init_cont_sin = np.array(init_cont_sin).ravel()

crit_strain_tri = np.array(crit_strain_tri).ravel()
init_cont_tri = np.array(init_cont_tri).ravel()

crit_strain_rand = np.array(crit_strain_rand).ravel()
init_cont_rand = np.array(init_cont_rand).ravel()


init_cont_all = np.concatenate((init_cont_sin,init_cont_tri, init_cont_rand))
crit_strain_all = np.concatenate((crit_strain_sin,crit_strain_tri, crit_strain_rand))

#----------Plot graph---------------#
colors = plt.cm.Reds(np.linspace(1,0,len(paths)+2)) # +2 added for color purposes

G = graph.create_fiber_network(nodes,edges)
pos = nx.get_node_attributes(G,'pos') # get nodal postitions

plt.figure(figsize=(3,3))

#nx.draw_networkx_nodes(G,pos, node_size = 10, node_color = 'r')
nx.draw_networkx_edges(G,pos, width = 0.75, edge_color = 'k', style = ':')
for count,path in enumerate(paths):
    path_edges = list(zip(path,path[1:]))
    nx.draw_networkx_edges(G,pos, edgelist=path_edges, width = 2.0, edge_color = colors[count], label = f'{count+1}')
plt.legend()
plt.savefig(f'n{int(sys.argv[1])}.pdf')

f = interp1d(init_cont_all, crit_strain_all, fill_value='extrapolate')
e_crits = f(contour_length_all)

plt.figure()
plt.plot(init_cont_all, crit_strain_all, ls = 'none', marker = 'o', fillstyle='none')
plt.plot(contour_length_all, e_crits, ls = 'none', marker='s')
plt.xlabel('tortuosity')
plt.ylabel('critical strain')
#-----------------Tangent Stiffness-------------------#
crit_strain=np.loadtxt(f_crit)[seed][int(sys.argv[3])]

E_t = np.gradient(force, strain) 
fig, ax1 = plt.subplots(figsize=(3,3))

ax1.plot(strain/crit_strain, force/EA, lw = 2, label = r'Normalized Force $\sfrac{\sigma_{yy}}{EA}$', c = 'k')
ax1.set_ylabel(r'Normalized Force $\sfrac{\sigma_{yy}}{E}$')
ax1.set_xlabel(r'$\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}$')

ax2 = ax1.twinx()


ax2.plot([], [], lw = 3, label = r'Normalized Force $\sfrac{\sigma_{yy}}{E}$', c = 'k')

ax2_c = '0.6'
ax2.plot(strain/crit_strain, E_t, lw = 2, c = ax2_c, label = r'Tangential Stiffness $E_t$')

ax2.set_xlabel(r'$\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}$')
ax2.ticklabel_format(axis='y',style='sci', scilimits = (0,0), useMathText=True)
ax2.get_yaxis().get_offset_text().set_position((1.2,-1.0))
ax2.set_ylabel(r'Tangential Stiffness $\sfrac{E_t}{E}$', color = ax2_c)
ax2.tick_params(axis='y', labelcolor=ax2_c)
plt.xlim([0,np.max(strain/crit_strain)])

print(e_crits)
e_crits, inverse = np.unique(e_crits, return_inverse=True)

print(inverse)
for num,eps in enumerate(e_crits):
    print(num,eps/crit_strain)
    if eps/crit_strain < 1.5:
        plt.axvline(eps/crit_strain, ls='--', c=colors[num], lw = 1.0)
        # plt.text(eps/crit_strain,0.1,f'{num+1}',size = 9, bbox=dict(boxstyle="round",
        #             facecolor = '0.9',
        #             edgecolor = 'k',
        #             ), transform=ax1.get_xaxis_transform())
plt.plot([],[],ls=':', c='0.5', lw = 0.75, label = r'Single Chain Critical Strain $\varepsilon_{crit}$') # for legeng


n = int(sys.argv[1])
kappa_tilde = [1e-3,1e-4,1e-5,1e-6][int(sys.argv[3])]

#plt.legend()
plt.title(rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tilde):.0f}}}$')
plt.tight_layout()


plt.savefig(f'fiber_recruitment_n{n}_kappa_tilde{kappa_tilde}.pdf')


plt.show()