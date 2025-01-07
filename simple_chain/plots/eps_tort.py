import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import optimize
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# plot settings
plt.style.use('jeff_style.mplstyle')

L = 10000

# Load files for random chains
w = [L/5,L/10,L/20,L/40,L/50,L/100]
n = [30,40,50,60,70]
w_r = [L/50,L/40,L/20,L/10]
crit_strain_random_all = []
init_cont_random_all = []
num_kappa_tilde = 4

# Load files for sinusoial chains
a = [250,500,1000]
kappa_tildes = [1e-3,1e-4,1e-5,1e-6]

crit_strain_sin = []
init_cont_sin = []

crit_strain_tri = []
init_cont_tri = []

crit_strain_rand = []
init_cont_rand = []

# Load all sinusoidal and triangular data
for ii in a:


    # Sinusoidal
    crit_strain_sin.append(np.loadtxt(f'../sinusoidal_chain/phase_diagram/crit_strain/a{ii}.txt'))
    init_cont_sin.append(np.loadtxt(f'../sinusoidal_chain/phase_diagram/init_cont/a{ii}.txt')/L)

    # triangular
    crit_strain_tri.append(np.loadtxt(f'../triangular_chain/phase_diagram/crit_strain/a{ii}.txt'))
    init_cont_tri.append(np.loadtxt(f'../triangular_chain/phase_diagram/init_cont/a{ii}.txt')/L)

# Load all random data
for ii in w_r:
    for jj in n:
        ii = int(ii)
        jj = int(jj)
        crit_strain_rand.append(np.loadtxt(f'../random_chain/phase_diagram/crit_strain/w{ii}/n{jj}.txt'))
        init_cont_rand.append(np.loadtxt(f'../random_chain/phase_diagram/init_cont/w{ii}/n{jj}.txt')/L)

crit_strain_rand = np.array(crit_strain_rand).swapaxes(0,2)
init_cont_rand = np.array(init_cont_rand).swapaxes(0,2)
crit_strain_tri = np.array(crit_strain_tri).swapaxes(0,2)
init_cont_tri = np.array(init_cont_tri).swapaxes(0,2)
crit_strain_sin=np.array(crit_strain_sin).swapaxes(0,2)
init_cont_sin=np.array(init_cont_sin).swapaxes(0,2)

c=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

tri1 = 2.5
tri2 = 1.0
slope  = np.array([tri1,tri2]) # to draw slope
mul = 2.0
triangle = np.array([[tri1,mul*tri1],[tri2,mul*tri2], [tri2,mul*tri1], [tri1, mul*tri1]])

print(triangle[:,0], triangle[:,1])
fig,ax1 = plt.subplots(figsize=(4,4))
min_x = np.min(init_cont_tri[-1,:].reshape(-1)-1)
max_x = np.max(init_cont_tri[-1,:].reshape(-1)-1)

min_y = np.min(crit_strain_tri[-1,:].reshape(-1))
max_y = np.max(crit_strain_tri[-1,:].reshape(-1))

ax1.loglog([min_x,max_x], [min_y,max_y], lw= 2, ls = '--', fillstyle = 'none', c='k', zorder=-1)
ax1.loglog(init_cont_tri[-1,:].reshape(-1)-1, crit_strain_tri[-1,:].reshape(-1), ls = 'None', marker = 's', fillstyle = 'none', c=(0.9,0.5,0.5), markeredgewidth=2.0, markersize = 10, label='Triangular')
ax1.loglog(init_cont_sin[-1,:].reshape(-1)-1, crit_strain_sin[-1,:].reshape(-1), ls = 'None', marker = 'o', fillstyle = 'none', c=(0.75,0,0), markeredgewidth=2.0, markersize = 7, label = 'Sinusoidal')
ax1.loglog(init_cont_rand[-1,:].reshape(-1)-1, crit_strain_rand[-1,:].reshape(-1), ls = 'None', marker = '*', fillstyle = 'none', c=(0.4,0,0), markeredgewidth=1.0, markersize = 10, zorder=0, label = 'Random')

ax1.plot(triangle[:,0], triangle[:,1], c='k', ls = ':', linewidth = 1.5)
ax1.fill_between(init_cont_tri[-1,:].reshape(-1)-1, crit_strain_tri[-1,:].reshape(-1), np.max(crit_strain_tri), facecolor='0.7', interpolate=True, zorder=-2)
ax1.fill_between(init_cont_tri[-1,:].reshape(-1)-1, crit_strain_tri[-1,:].reshape(-1), np.min(crit_strain_tri), facecolor='0.4', interpolate=True, zorder=-2)


ax1.set_xlabel(r'(Shifted) Initial Tortuosity $\tau_0-1$')
ax1.set_ylabel(r'Applied Strain $\varepsilon_{applied}$')
plt.xlim([min_x,max_x])
plt.ylim([min_y,max_y])
plt.text(2.5,0.5,'Stretch Dominated', fontsize=10)
plt.text(0.2,3.8,'Bend Dominated', fontsize=10)

plt.legend(loc = 'lower right', title=r'Critical Strain $\varepsilon_{crit}$',fontsize=10, title_fontsize=10)
plt.tight_layout()
plt.savefig('eps_tort.pdf')

print(init_cont_rand[-1,:].shape)
plt.show()