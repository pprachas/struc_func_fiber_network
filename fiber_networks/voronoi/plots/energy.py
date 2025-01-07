import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import pathlib
from matplotlib.ticker import PercentFormatter
from scipy.spatial.distance import jensenshannon
#----------Import functions---------------------#
import sys
import os
sys.path.append('../../../utils')
import graph

plt.style.use('jeff_style.mplstyle')

# fiber network parameters
W = 1e4
H = 1e4
n = 500
E = 1 # Stiffness modulus
seed = 1
L_0 = 10000
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]
num_kappa_tilde = 1

theta_N = []
strain = []
L_C_all = []
bend_energy = []
stretch_energy = []
shear_energy = []
total_energy = []
normalized_force = []


crit_strain = np.loadtxt(f'../phase_diagram/crit_strain/n{int(n)}.txt')[seed][:]

for kappa_tilde in kappa_tildes:
    f_root = f'../fea_run/n{n}/kappa_tilde{kappa_tilde}/seed{seed}'
    disp = np.loadtxt(f'{f_root}/disp.txt')
    cont = np.loadtxt(f'{f_root}/cont.txt')
    bend = np.loadtxt(f'{f_root}/bend_total.txt')
    stretch = np.loadtxt(f'{f_root}/stretch_total.txt')
    shear = np.loadtxt(f'{f_root}/shear_total.txt')
    force = np.loadtxt(f'{f_root}/force.txt')

    orientation = np.loadtxt(f'{f_root}/rot.txt')


    L_C_all.append(cont)
    strain.append(disp/L_0)
    bend_energy.append(bend)
    stretch_energy.append(stretch)
    shear_energy.append(shear)
    total_energy.append(bend+stretch+shear)

stretch_energy_last = np.array(stretch_energy[-1])
bend_energy_last = np.array(bend_energy[-1])
shear_energy_last = np.array(shear_energy[-1])

total_energy_last = np.array(total_energy[-1])

percent_stretch = stretch_energy_last/total_energy_last
percent_bend = bend_energy_last/total_energy_last
percent_shear = shear_energy_last/total_energy_last

#--------------Plot Energy Distrivution-------------------#
plt.figure(figsize=(4,3))
plt.gca().set_prop_cycle(None)
color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]


# prepend for plotting colors
strain_plot = strain[-1]
strain_plot = np.insert(strain_plot,0,0)

y_min = -0.02
y_max = 1.05

plt.fill_between(strain_plot, y_min, y_max, where = strain_plot <= crit_strain[-1], color = '0.40')
plt.fill_between(strain_plot, y_min, y_max, where = strain_plot >= crit_strain[-1], color = '0.70')


plt.xlim([0,strain_plot[-1]])
plt.ylim([-0.02,1.01])

plt.xlabel(r'Applied Strain $\displaystyle{\varepsilon_{applied}}$')
plt.ylabel(r'Strain Energy Ratio')
plt.title('Dense Fiber Network')

plt.plot(strain[-1], percent_bend,  label = 'Bending', lw = 4, c = 'k')
plt.plot(strain[-1], percent_stretch,  label = 'Stretching', lw = 4, c = (0.9,0.5,0.5))
plt.plot(strain[-1], percent_shear,  label = 'Shear', lw = 4, c = (0.5,0,0))
plt.axvline(crit_strain[-1], ls='--', lw=3, color='steelblue', label=r'Critical Strain $\varepsilon_{crit}$')
plt.tight_layout()
plt.savefig('energy_dist_fiber.pdf')
#----------\tilde{\kappa} stretching ------------------#
color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

print(stretch_energy)
plt.figure(figsize=(4,3))
for ii in range(len(kappa_tildes)):
    stretch_energy_plot = np.array(stretch_energy[ii])

    total_energy_plot = np.array(total_energy[ii])

    percent_stretch = stretch_energy_plot/total_energy_plot


    plt.plot(strain[ii]/crit_strain[ii], percent_stretch, color=color[ii], lw=4)

plt.xlim([0,1.51])
plt.ylim([-0.02,1.01])
plt.xlabel(r'Normalized Applied Strain $\displaystyle{\sfrac{\varepsilon_{applied}}{\varepsilon_{crit}}}$')
plt.ylabel('Stretching Energy Ratio')
plt.title('Dense Network Stretching Energy Ratio')
plt.tight_layout()
plt.savefig('energy_kappa_fiber.pdf')
plt.show()