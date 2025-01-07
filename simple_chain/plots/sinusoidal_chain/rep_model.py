import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes
import matplotlib.gridspec as gridspec

#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
lmbda = L_0/5 # wavelength 
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
    link_points = np.loadtxt(f'../../sinusoidal_chain/link_points/a{int(int(amplitude))}/lmbda{int(lmbda)}.txt')

    G = graph.create_chain_graph(link_points)

    edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
    edge_dist = np.array(list(edge_dist))
    init_contour = np.sum(edge_dist)
    L_char = (lmbda/(2*L_0))*init_contour
    r = 2*L_char*np.sqrt(kappa_tilde)
    S = np.pi*r**2
    I = np.pi*r**4/4
    f_prefix_name = f'../../sinusoidal_chain/fea_run/a{int(amplitude)}/kappa_tilde{kappa_tilde}/lmbda{int(lmbda)}'
    disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
    cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
    bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
    stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
    shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
    force = np.loadtxt(f'{f_prefix_name}/force.txt')

    orientation = np.loadtxt(f'{f_prefix_name}/rot.txt')
    crit_strain = np.loadtxt(f'../../sinusoidal_chain/phase_diagram/crit_strain/a{int(amplitude)}.txt')
    crit_strains.append(crit_strain[num_run][count])

    strain.append(disp/L_0)
    straightness.append(cont/L_0) #straightness.append((L+disp)/cont)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)
    normalized_force.append(force/(E*S))

init_straightness = init_contour/L_0
init_theta = np.arccos(init_straightness**(-1))

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

plt.figure(figsize=(3,3))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Sinusoidal Chain (Log scale)')
#plt.plot(straightness[0], normalized_force[0], c = 'r')
for ii in range(len(kappa_tildes)):
    # stretching component
    stretching = (straightness[ii]/init_straightness)-1
    L = L_0*(1+strain[ii])

    #bending component
    sin_theta = np.sqrt((straightness[ii]*L_0)**2-L**2)/(straightness[ii]*L_0)
    theta = np.arccos(L/(straightness[ii]*L_0))
    bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta

    plt.loglog(strain[ii]/crit_strains[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 1.75, c = color[ii])
    plt.loglog(strain[ii]/crit_strains[ii],stretching, lw = 0.5, c = color[ii], marker = 'o', fillstyle='none', markersize=5, markeredgewidth = 0.5)
    plt.loglog(strain[ii]/crit_strains[ii],stretching+bending, lw = 0.5, c = color[ii], marker = '*', fillstyle='none', markersize=5, markeredgewidth = 0.5)
    #plt.loglog(strain[ii],bending, marker = 'None', ls = ':',c = 'g', lw = 2.0, dashes=(2, 0.5))
#plt.plot(strain[ii],straightness[ii]/straightness[ii][0]-1, marker = 'None', ls = ':',c = 'orange', lw = 2.0, dashes=(2, 0.5), label = 'Anal
# ytical') # for label
# plt.loglog(strain[ii]/crit_strains[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0 ,0.5),label = 'Reduced order model (stretching)')
# plt.loglog(strain[ii]/crit_strains[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1),label = 'Reduced order model (stretching + bending)')

plt.tight_layout()
plt.savefig('sinusoidal_rom.pdf')
plt.figure(figsize=(7,7))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Triangular Chain (Linear scale)')
#plt.plot(straightness[0], normalized_force[0], c = 'r')
for ii in range(len(kappa_tildes)):
    # stretching component
    stretching = (straightness[ii]/init_straightness)-1
    L = L_0*(1+strain[ii])

    #bending component
    sin_theta = np.sqrt((straightness[ii]*L_0)**2-L**2)/(straightness[ii]*L_0)
    theta = np.arccos(L/(straightness[ii]*L_0))
    bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta

    error_stretch = np.mean(np.abs((normalized_force[ii] - stretching)/normalized_force[ii]))

    error_bend_stretch = np.mean(np.abs((normalized_force[ii] - (stretching+bending))/normalized_force[ii]))

    print(rf'MAPE kappa = {kappa_tildes[ii]:2e} | stretch only: {error_stretch}| stretch+bend: {error_bend_stretch}|')

    plt.plot(strain[ii],normalized_force[ii], label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', marker = 'None', lw = 3)
    plt.plot(strain[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0, 0.5))
    plt.plot(strain[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1))
    #plt.loglog(strain[ii],bending, marker = 'None', ls = ':',c = 'g', lw = 2.0, dashes=(2, 0.5))
#plt.plot(strain[ii],straightness[ii]/straightness[ii][0]-1, marker = 'None', ls = ':',c = 'orange', lw = 2.0, dashes=(2, 0.5), label = 'Anal
# ytical') # for label
plt.plot(strain[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0 ,0.5),label = 'Representative model (stretching)')
plt.plot(strain[ii],stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1),label = 'Representative model (stretching + bending)')

plt.legend()
plt.savefig('Sinusoidal_rep_linear.pdf')


plt.figure()
pos = nx.get_node_attributes(G, 'pos')
nx.draw_networkx_edges(G,pos,edge_color='k')
plt.axis('equal')
plt.savefig('sin_chain.pdf')
plt.show()