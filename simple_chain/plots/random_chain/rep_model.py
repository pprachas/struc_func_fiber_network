import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes
import matplotlib.gridspec as gridspec
from scipy.stats import skew

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
w = L_0/10 # width
num_link = 30 #amplitude
E = 1 # Stiffness moduli
num_run = 4

theta_N = []
strain = []
L_C = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []
normalized_force = []
crit_strains = []

# Create graph
link_points = np.loadtxt(f'../../random_chain/link_points/w{int(w)}/n{int(num_link)}/link_points{num_run}.txt')

G = graph.create_chain_graph(link_points)

edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour = np.sum(edge_dist)

init_theta_N = nx.get_edge_attributes(G,'orientation').values()
init_theta_N = np.pi/2-np.array(list(init_theta_N))
L_C0 = init_contour/L_0

L_char = np.mean(edge_dist)


# Load data
for count,kappa_tilde in enumerate(kappa_tildes):
    f_prefix_name = f'../../random_chain/fea_run/w{int(w)}/n{num_link}/kappa_tilde{kappa_tilde}/seed{int(num_run)}'
    disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
    cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
    bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
    stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
    shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
    force = np.loadtxt(f'{f_prefix_name}/force.txt')

    orientation = np.loadtxt(f'{f_prefix_name}/rot.txt')

    crit_strain = np.loadtxt(f'../../random_chain/phase_diagram/crit_strain/w{int(w)}/n{num_link}.txt')
    crit_strains.append(crit_strain[num_run][count])

    r = 2*L_char*np.sqrt(kappa_tilde)
    S = np.pi*r**2
    I = np.pi*r**4/4

    L_C.append(cont)
    theta_N.append(orientation)
    strain.append(disp/L_0)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)
    normalized_force.append(force/(E*S))

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

plt.figure(figsize=(3,3))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Random Chain (Log scale)')

for ii in range(len(kappa_tildes)):
    # stretching component
    stretching = (L_C[ii]/init_contour)-1
    L = L_0*(1+strain[ii])
    init_theta = np.arccos(L_C0**(-1))

    #bending component
    sin_theta = np.sqrt((L_C[ii])**2-L**2)/(L_C[ii])
    theta = np.arccos(L/(L_C[ii]))
    bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta

    
    plt.loglog(strain[ii]/crit_strains[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 1.75, c = color[ii])
    plt.loglog(strain[ii]/crit_strains[ii],stretching, lw = 0.5, c = color[ii], marker = 'o', fillstyle='none', markersize=5, markeredgewidth = 0.5)
    plt.loglog(strain[ii]/crit_strains[ii],stretching+bending, lw = 0.5, c = color[ii], marker = '*', fillstyle='none', markersize=5, markeredgewidth = 0.5)

plt.tight_layout()
plt.savefig('random_rom.pdf')

plt.figure(figsize=(7,7))
plt.loglog(strain[ii], bending)

plt.figure()
plt.plot((L_char/edge_dist))


plt.figure()
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Random Chain (Linear scale)')

for ii in range(len(kappa_tildes)):
    # stretching component
    stretching = (L_C[ii]/init_contour)-1
    L = L_0*(1+strain[ii])

    #bending component
    sin_theta = np.sqrt((L_C[ii])**2-L**2)/(L_C[ii])
    theta = np.arccos(L/(L_C[ii]))
    bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta

    plt.plot(strain[ii]/crit_strains[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 3)
    plt.plot(strain[ii]/crit_strains[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0, 0.5))
    plt.plot(strain[ii]/crit_strains[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1))

plt.plot(strain[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0 ,0.5),label = 'Representative model (stretching)')
plt.plot(strain[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1),label = 'Representative model (stretching + bending)')

plt.legend()
plt.savefig('random_rep_linear.pdf')

plt.figure()
pos = nx.get_node_attributes(G, 'pos')
nx.draw_networkx_edges(G,pos,edge_color='k')
plt.axis('equal')
plt.savefig('random_chain.pdf')

plt.show()
