import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L = 10000 # total length of fiber
lmbda = L/10 # wavelength 
amplitude = L/40 #amplitude
E = 1e6 # Stiffness moduli

strain = []
straightness = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []
normalized_force = []

# Load data
for kappa_tilde in kappa_tildes:
    link_points = np.loadtxt(f'../../sinusoidal_chain/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')

    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    init_contour = np.sum(edge_dist)
    L_char = (lmbda/L)*init_contour
    r = 2*L_char*np.sqrt(kappa_tilde)
    S = np.pi*r**2
    
    f_prefix_name = f'sinusoidal_{kappa_tilde:.0e}'
    disp = np.loadtxt(f'data/{f_prefix_name}/disp.txt')
    cont = np.loadtxt(f'data/{f_prefix_name}/cont.txt')
    bend = np.loadtxt(f'data/{f_prefix_name}/bend.txt')
    stretch = np.loadtxt(f'data/{f_prefix_name}/stretch.txt')
    shear = np.loadtxt(f'data/{f_prefix_name}/shear.txt')
    force = np.loadtxt(f'data/{f_prefix_name}/force.txt')

    strain.append(disp/L)
    straightness.append(cont/L)#straightness.append((L+disp)/cont)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)
    normalized_force.append(force/(E*S))


# Plot coutour-displacement
fig, ax=plt.subplots(figsize=(9,7))
#ax.set_ylabel(r'$test$')
ax.set_xlabel(r'$\varepsilon_{applied}$')
ax.set_ylabel(r'$\frac{L_C}{L}$')

axins = zoomed_inset_axes(ax, zoom = 3.0, loc = 2, bbox_to_anchor=(0.12,0,1,1), bbox_transform=ax.transAxes)
for ii in range(len(kappa_tildes)):
    # main plot
    ax.plot(strain[ii],straightness[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None')
    
    # inset plot
    axins.plot(strain[ii],straightness[ii])
    axins.set_xlabel(r'$\varepsilon_{applied}$')
    axins.set_ylabel(r'$\frac{L_C}{L}$')
    
    
    axins.set_xlim(0.35,0.5)
    axins.set_ylim(1.455, 1.55)
    


mark_inset(ax,axins, loc1 = 3, loc2 = 4, edgecolor='k', ls = ':')
ax.legend(loc = 'lower right')

plt.savefig('sinusoidal_cont.pdf')
plt.show()