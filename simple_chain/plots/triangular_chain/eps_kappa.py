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
plt.style.use('jeff_style.mplstyle')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')


kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

L_0 = 10000 # total length of fiber
lmbda = L_0/5 # wavelength 
num_run = 1
amplitude = [L_0/10,L_0/20,L_0/40][num_run] #amplitude

crit_strains = []

# Load data
for count,kappa_tilde in enumerate(kappa_tildes):

    crit_strain = np.loadtxt(f'../../triangular_chain/phase_diagram/crit_strain/a{int(amplitude)}.txt')
    crit_strains.append(crit_strain[num_run][count])

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

plt.figure(figsize=(3,3))
plt.title(r'Triangular Chains (Semi-log $x$)')
plt.xlabel(r'Dimensionless Bending length $\tilde{\kappa}$')
plt.ylabel(r'Applied Strain $\varepsilon_{applied}$')
plt.semilogx(kappa_tildes,crit_strains, marker = 'o', color='k', markersize=10, ls = '--')
plt.fill_between(kappa_tildes,crit_strains, crit_strains[0], facecolor='0.7', interpolate=True)
plt.fill_between(kappa_tildes,crit_strains, crit_strains[-1], facecolor='0.4', interpolate=True)

plt.xlim([kappa_tildes[-1], kappa_tildes[0]])
plt.ylim([crit_strains[-1], crit_strains[0]])

plt.tight_layout()
plt.savefig('triangular_kappa.pdf')

plt.figure()

line = Line2D([0], [0], color='k', lw = 3.0, ls='none', marker='o', label=r'Critical Strain $\varepsilon_{crit}$')

plt.legend(handles=[line])

plt.savefig('crit_strain_legend.pdf')
plt.show()