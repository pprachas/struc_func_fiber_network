import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('../jeff_style.mplstyle')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


kappa_tildes = [1e-3,1e-4,1e-5,1e-6]


L_0 = 10000 # total length of fiber
lmbda = L_0/10 #L_0/np.array([1,5,10,20,40]) # wavelength
num_run = 1
amplitude = [L_0/10,L_0/20,L_0/40][num_run] #amplitude

E = 1.0 # Stiffness moduli

strain = []
straightness = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []
normalized_force = []
crit_strains = []

# Load data
for count,kappa_tilde in enumerate(kappa_tildes):
    link_points = np.loadtxt(f'../../discretized_sin/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')

    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    init_contour = np.sum(edge_dist)
    L_char = np.mean(edge_dist)
    r = 2*L_char*np.sqrt(kappa_tilde)
    S = np.pi*r**2
    I = np.pi*r**4/4
    f_prefix_name = f'../../discretized_sin/fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'
    disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
    cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
    bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
    stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
    shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
    force = np.loadtxt(f'{f_prefix_name}/force.txt')

    orientation = np.loadtxt(f'{f_prefix_name}/rot.txt')
    crit_strain = np.loadtxt(f'../../discretized_sin/phase_diagram/crit_strain/a{int(amplitude)}.txt')
    crit_strains.append(crit_strain[num_run][count])

    strain.append(disp/L_0)
    straightness.append(cont/L_0)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)
    normalized_force.append(force/(E*S))

init_straightness = init_contour/L_0
init_theta = np.arccos(init_straightness**(-1))

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

plt.figure(figsize=(4,4))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Discretized Sinusoid (Log scale)')
#plt.plot(straightness[0], normalized_force[0], c = 'r')
for ii in range(len(kappa_tildes)):
    # stretching component
    # stretching component
    stretching = (straightness[ii]/init_straightness)-1
    L = L_0*(1+strain[ii])

    #bending component
    sin_theta = np.sqrt((straightness[ii]*L_0)**2-L**2)/(straightness[ii]*L_0)
    theta = np.arccos(L/(straightness[ii]*L_0))
    bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta

    area_stretch = np.trapz(np.abs(normalized_force[ii] - stretching), strain[ii]/crit_strains[ii])
    area_stretch_bend = np.trapz(np.abs(normalized_force[ii] - (stretching+bending)), strain[ii]/crit_strains[ii])

    print(area_stretch)
    plt.loglog(strain[ii]/crit_strains[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 1.75, c = color[ii])
    plt.loglog(strain[ii]/crit_strains[ii],stretching, lw = 0.5, c = color[ii], marker = 'o', fillstyle='none', markersize=5, markeredgewidth = 0.5)
    plt.loglog(strain[ii]/crit_strains[ii],stretching+bending, lw = 0.5, c = color[ii], marker = '*', fillstyle='none', markersize=5, markeredgewidth = 0.5)

    #plt.loglog(strain[ii],bending, marker = 'None', ls = ':',c = 'g', lw = 2.0, dashes=(2, 0.5))
#plt.plot(strain[ii],straightness[ii]/straightness[ii][0]-1, marker = 'None', ls = ':',c = 'orange', lw = 2.0, dashes=(2, 0.5), label = 'Anal
# ytical') # for label
# plt.loglog(strain[ii]/crit_strains[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0 ,0.5),label = 'Reduced order model (axial)')
# plt.loglog(strain[ii]/crit_strains[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1),label = 'Reduced order model (axial + rotation)')
plt.tight_layout()
plt.savefig('discretize_sin_rom.pdf')
plt.figure()

# for kappa
line1 = Line2D([0], [0], label=r'$10^{-3}$', color=(0.2,0.2,0.2), lw = 3)
line2 = Line2D([0], [0], label=r'$10^{-4}$', color=(0.5,0.5,0.5), lw = 3)
line3 = Line2D([0], [0], label=r'$10^{-5}$', color=(0.9,0.5,0.5), lw = 3)
line4 = Line2D([0], [0], label=r'$10^{-6}$', color=(0.5,0.0,0.0), lw = 3)

handles =[line1, line2, line3, line4]

plt.legend(handles=handles, title = r'Dimensionless Bending Length $\tilde{\kappa}$', ncol=4)
plt.savefig('kappa_legend.pdf')

plt.figure()

line5 = Line2D([0], [0], label='FEA', color='k', lw = 3.0, fillstyle='none')
line6 = Line2D([0], [0], label='axial',marker = 'o', c = 'k', lw = 0.5, fillstyle='none')
line7 = Line2D([0], [0], label='axial + rotation',marker = '*', c = 'k', lw = 0.5, fillstyle='none')

handles =[line5, line6, line7]

plt.legend(handles=handles, title = r'Model Type', ncol=3)
plt.savefig('model_legend.pdf')

plt.figure()
pos = nx.get_node_attributes(G, 'pos')
nx.draw_networkx_edges(G,pos,edge_color='k')
plt.axis('equal')
plt.savefig('discrite_chain.pdf')
plt.show()