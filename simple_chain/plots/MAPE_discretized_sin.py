import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition, mark_inset, zoomed_inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pathlib
from scipy.stats import entropy


#--------Import functions---------------------#
import sys
sys.path.append('../../utils')
import graph

plt.style.use('jeff_style.mplstyle')

kappa_tilde = 1e-6

L_0 = 10000 # total length of fiber
lmbdas = L_0/np.array([1,5,10,20,40]) # wavelength 
amplitudes = [L_0/10,L_0/20,L_0/40] #amplitude
E = 1e6 # Stiffness moduli

# marker = ['s','D','o','*']

# plt.figure()
# plt.xlabel(r'$\theta_0$')
# plt.ylabel(r'MAPE')
# plt.title('Triangular Chains')
# Load data

link_points = np.loadtxt(f'../discretized_sin/link_points.txt')

G = graph.create_chain_graph(link_points)

pos=nx.get_node_attributes(G,'pos')
orientation = np.fromiter(nx.get_edge_attributes(G,'orientation').values(), dtype=float)

edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))
init_contour = np.sum(edge_dist)
L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)
S = np.pi*r**2
I = np.pi*r**4/4
f_prefix_name = '../discretized_sin/fea_run'
disp = np.loadtxt(f'{f_prefix_name}/disp.txt')
cont = np.loadtxt(f'{f_prefix_name}/cont_total.txt')
bend = np.loadtxt(f'{f_prefix_name}/bend_total.txt')
stretch = np.loadtxt(f'{f_prefix_name}/stretch_total.txt')
shear = np.loadtxt(f'{f_prefix_name}/shear_total.txt')
force = np.loadtxt(f'{f_prefix_name}/force.txt')
bend_ind = np.loadtxt(f'{f_prefix_name}/bend_energy.txt')

strain = disp/L_0
straightness = cont/L_0 #straightness.append((L+disp)/cont)
bend_energy = bend
stretch_energy = stretch
shear_energy = shear
total_energy = bend+stretch+shear
normalized_force = force/(E*S)


init_straightness = init_contour/L_0
init_theta = np.arccos(init_straightness**(-1))

# stretching component
stretching = (straightness/init_straightness)-1
L = L_0*(1+strain)

#bending component
sin_theta = np.sqrt((straightness*L_0)**2-L**2)/(straightness*L_0)
theta = np.arccos(L/(straightness*L_0))
bending = 3*kappa_tilde*(init_theta-theta)*sin_theta

error_stretch = np.mean(np.abs((normalized_force - stretching)/normalized_force))


error_bend_stretch = np.mean(np.abs((normalized_force - (stretching+bending))/normalized_force))

plt.figure()
plt.loglog(strain, normalized_force, label = 'FEA solution')
plt.loglog(strain,stretching, marker = 'None', ls = ':',c = 'chocolate', lw = 1.75, dashes=(1.0, 0.5), label = 'Representative model (stretching)')
plt.loglog(strain,stretching+bending, marker = 'None', ls = ':',c = 'royalblue', lw = 1.75, dashes=(3, 1), label = 'Representative model (stretching+bending)')

print(f'MAPE stretch: {error_stretch}')
print(f'MAPE stretch: {error_bend_stretch}')

print(bend_ind[0].shape)
plt.legend()


plt.figure(figsize=(5,5))
edges = nx.draw_networkx_edges(G,pos, width = 3.5, edge_color = bend_ind[0], edge_cmap = plt.cm.OrRd, edge_vmin = 0, edge_vmax = bend_ind[0].max())
plt.colorbar(edges, label = 'Bending Energy')

# orientation
bend_plot = bend_ind[0]
norm = mpl.colors.Normalize(vmin=0, vmax=bend_plot.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.OrRd)
cmap.set_array([])
plt.figure(figsize=(7,7))
ax = plt.subplot(111,polar=True)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_ylim([0,0.15])
for ii in range(len(orientation)):
    ax.plot([0,orientation[ii]],[0,edge_dist[ii]/L_0], c = cmap.to_rgba(bend_plot[ii]), lw = 3)
cbar = plt.colorbar(cmap)
cbar.set_label('Bending Energy')

plt.show()
