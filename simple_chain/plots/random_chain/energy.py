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
E = 1.0 # Stiffness moduli
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
for kappa_tilde in kappa_tildes:
    f_prefix_name = f'../../random_chain/fea_run/w{int(w)}/n{num_link}/kappa_tilde{kappa_tilde}/seed{int(num_run)}'
    disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
    cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
    bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
    stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
    shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
    force = np.loadtxt(f'{f_prefix_name}/force.txt')
    crit_strain = np.loadtxt(f'../../random_chain/phase_diagram/crit_strain/w{int(w)}/n{num_link}.txt')
    crit_strains.append(crit_strain[num_run][-1])


    orientation = np.loadtxt(f'{f_prefix_name}/rot.txt')

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

stretch_energy = np.array(stretch_energy)
bend_energy = np.array(bend_energy)
shear_energy = np.array(shear_energy)

total_energy = np.array(total_energy)

percent_stretch = stretch_energy/total_energy
percent_bend = bend_energy/total_energy
percent_shear = shear_energy/total_energy

#-------------Plot Force-Displacement Curve-------------------#
# plt.figure(figsize=(5,3))
# plt.gca().set_prop_cycle(None)
# plt.xlabel(r'$\displaystyle{\varepsilon_{applied}}$')
# plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
# plt.title('Random Chain Force-Displacement Curve')
# #plt.plot(straightness[0], normalized_force[0], c = 'r')

#----------Energy----------------------#
# y_min = -0.005
# y_max = 0.20

# ii = 0
# # stretching component
# stretching = (L_C[0]/init_contour)-1
# L = L_0*(1+strain[0])
# init_theta = np.arccos(L_C0**(-1))

# #bending component
# sin_theta = np.sqrt((L_C[ii])**2-L**2)/(L_C[ii])
# theta = np.arccos(L/(L_C[ii]))
# bending = 3*kappa_tildes[ii]*(init_theta-theta)*sin_theta


# # bending = np.sum(3*(kappa_tildes[ii]*delta_theta*cos_theta)*(L_char/edge_dist), axis = 1)#np.sum(3*kappa_tildes[ii]*delta_theta*cos_theta, axis = 1) #*init_contour/L_char
# plt.plot(strain[ii],normalized_force[ii], label = rf'FE Model', marker = 'None', lw = 5)
# #plt.loglog(strain[ii],stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0, 0.5))
# plt.plot(strain[ii],stretching+bending, marker = 'None', label = 'Reduced Order Model',ls = ':', lw = 4, dashes=(3, 1), c = (0.9,0.5,0.5))
# plt.fill_between(strain[0], y_min, y_max, where = strain[0] <= crit_strains[0], color = '0.95')
# plt.fill_between(strain[0], y_min, y_max, where = strain[0] >= crit_strains[0], color = '0.70')
# plt.xlim([0,strain[0][-1]])
# plt.ylim([y_min,y_max])
# #plt.loglog(strain[ii],bending, marker = 'None', ls = ':',c = 'g', lw = 2.0, dashes=(2, 0.5))
# #plt.plot(strain[ii],straightness[ii]/straightness[ii][0]-1, marker = 'None', ls = ':',c = 'orange', lw = 2.0, dashes=(2, 0.5), label = 'Anal
# # ytical') # for label

# plt.legend(loc='center left')
# plt.tight_layout()
# plt.savefig('random_fd_newmech.pdf')

#------Plot Energy Distribution----------------------#
y_min = -0.02
y_max = 1.05

plt.figure(figsize = (4,3))
plt.gca().set_prop_cycle(None)
plt.xlabel(r'Applied Strain $\displaystyle{\varepsilon_{applied}}$')
plt.ylabel(r'Strain Energy Ratio')
plt.title('Random Chain Strain Energy Ratio')
#plt.plot(strain[0], percent_shear[0], label = 'Shear')
plt.plot(strain[-1], percent_bend[-1],  label = 'Bending', lw = 4, c = 'k')
plt.plot(strain[-1], percent_stretch[-1],  label = 'Stretching', lw = 4, c = (0.9,0.5,0.5))
plt.plot(strain[-1], percent_shear[-1],  label = 'Shear', lw = 4, c = (0.5,0,0))
plt.axvline(crit_strains[-1], ls='--', lw=3, color='steelblue', label=r'Critical Strain $\varepsilon_{crit}$')

# prepend for plotting colors
strain_plot = strain[-1]
strain_plot = np.insert(strain_plot,0,0)

plt.fill_between(strain_plot, y_min, y_max, where = strain_plot <= crit_strains[-1], color = '0.40')
plt.fill_between(strain_plot, y_min, y_max, where = strain_plot >= crit_strains[-1], color = '0.70')
plt.xlim([0,strain_plot[-1]])
plt.ylim([-0.02,1.01])
plt.legend(loc='center left')
plt.tight_layout()
plt.savefig('energy_dist.pdf')

#---------\tilde{\kappa} stretching ------------------#
color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

plt.figure(figsize=(4,3))
for ii in range(len(kappa_tildes)):
    plt.plot(strain[ii]/crit_strains[ii], percent_stretch[ii], color=color[ii], lw=4)


plt.xlim([0,1.51])
plt.ylim([-0.02,1.01])
plt.xlabel(r'Normalized Applied Strain $\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel('Stretching Energy Ratio')
plt.title('Random Chain Stretching Energy Ratio')
plt.tight_layout()
plt.savefig('energy_kappa.pdf')

plt.show()


