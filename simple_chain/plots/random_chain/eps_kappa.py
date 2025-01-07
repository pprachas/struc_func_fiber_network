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
w = L_0/10 # width
num_link = 30 #amplitude
E = 1 # Stiffness moduli
num_run = 4

crit_strains = []


# Load data
for count,kappa_tilde in enumerate(kappa_tildes):

    crit_strain = np.loadtxt(f'../../random_chain/phase_diagram/crit_strain/w{int(w)}/n{num_link}.txt')
    crit_strains.append(crit_strain[num_run][count])

plt.figure(figsize=(3,3))
plt.title(r'Random Chains (Semi-log $x$)')
plt.xlabel(r'Dimensionless Bending length $\tilde{\kappa}$')
plt.ylabel(r'Critical Strain $\varepsilon_{crit}$')

plt.semilogx(kappa_tildes,crit_strains, ls='--', marker = 'o')
plt.fill_between(kappa_tildes,crit_strains, crit_strains[0], facecolor='0.7', interpolate=True)
plt.fill_between(kappa_tildes,crit_strains, crit_strains[-1], facecolor='0.4', interpolate=True)

plt.xlim([kappa_tildes[-1], kappa_tildes[0]])
plt.ylim([crit_strains[-1], crit_strains[0]])
plt.tight_layout()
plt.savefig('random_kappa.pdf')
plt.show()