import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import sys
#----------Import functions---------------------#
sys.path.append('../../utils')
import graph

# plot settings
plt.style.use('jeff_style.mplstyle')

n = 100
seed = 0
E = 1e6 # Stiffness moduli


L_0 = 10000 # total length of fiber
kappa_tilde = 1e-6

# Construct graphs
nodes = np.loadtxt(f'graph/nodes/n{n}/seed{seed}.txt')
edges = np.loadtxt(f'graph/edges/n{n}/seed{seed}.txt', dtype = 'int')
G = graph.create_fiber_network(nodes,edges)
pos = nx.get_node_attributes(G,'pos')
edge_dist = nx.get_edge_attributes(G,'dist').values() # get edge weights
edge_dist = np.array(list(edge_dist))

L_char = np.mean(edge_dist)
r = 2*L_char*np.sqrt(kappa_tilde)
S = np.pi*r**2
I = np.pi*r**4/4

# path to file
f_path = f'plots/n{n}/kappa{kappa_tilde}/seed{seed}/quad3'


fiber_force = np.loadtxt(f'{f_path}/force.txt')/(E*r*L_0)
fiber_strain = np.loadtxt(f'{f_path}/disp.txt')/L_0
fiber_stretch_total = np.loadtxt(f'{f_path}/stretch_total.txt')
fiber_bend_total = np.loadtxt(f'{f_path}/bend_total.txt')
fiber_shear_total = np.loadtxt(f'{f_path}/shear_total.txt')
fiber_total_energy = fiber_stretch_total+fiber_bend_total+fiber_shear_total
fiber_stretch = np.loadtxt(f'{f_path}/stretch_energy.txt')
fiber_bend = np.loadtxt(f'{f_path}/bend_energy.txt')
fiber_shear = np.loadtxt(f'{f_path}/shear_energy.txt')
fiber_energy = fiber_stretch+fiber_bend+fiber_shear

crit_strain = 0.12629596522176872

#-----------------------Plot f-d curve-------------------#
plt.figure(figsize=(5,3))
y_min = -1e-6
y_max = fiber_force[-1]
plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\varepsilon_{applied}}$')
plt.ylabel(r'$\sfrac{\sigma_{yy}}{E}$')
plt.title('Fiber Network Force-Displacement Curve')
plt.plot(fiber_strain,fiber_force, lw =5)
plt.fill_between(fiber_strain, y_min, y_max, where = fiber_strain <= crit_strain, color = '0.95')
plt.fill_between(fiber_strain, y_min, y_max, where = fiber_strain >= crit_strain, color = '0.70')
plt.xlim([0,fiber_strain[-1]])
plt.ylim([y_min,y_max])
plt.tight_layout()

plt.savefig('f-d_curve.pdf')

#----------------------Plot Energy-----------------#
plt.figure(figsize=(5,3))
y_min = -0.02
y_max = 1.05

plt.gca().set_prop_cycle(None)
plt.xlabel(r'$\displaystyle{\varepsilon_{applied}}$')
plt.ylabel(r'Strain Energy Ratio')
plt.title('Fiber Network Energy Ratio')
plt.plot(fiber_strain,fiber_bend_total/fiber_total_energy,  lw = 5, c = 'k', label = 'Bending')
plt.plot(fiber_strain,fiber_stretch_total/fiber_total_energy, lw = 5, label = 'Stretching', c = (0.9,0.5,0.5))
plt.plot(fiber_strain,fiber_shear_total/fiber_total_energy,  lw = 5, c = (0.5,0,0), label = 'Shear')
plt.fill_between(fiber_strain, y_min, y_max, where = fiber_strain <= crit_strain, color = '0.95')
plt.fill_between(fiber_strain, y_min, y_max, where = fiber_strain >= crit_strain, color = '0.70')
plt.xlim([0,fiber_strain[-1]])
plt.ylim([y_min,y_max])
plt.legend(loc='center left')
plt.tight_layout()

plt.savefig('energy.pdf')

plt.figure(figsize=(7,7))
nx.draw_networkx_edges(G,pos, width = 2, edge_color = 'k')
plt.savefig('graph.pdf')

idx = np.argmin(np.abs(fiber_strain-crit_strain))

print(idx)