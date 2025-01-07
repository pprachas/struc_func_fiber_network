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
a = [250,500,1000] #[250,500,1000]
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

tri1 = 1.0
tri2 = 2.0
slope  = np.array([tri1,tri2]) # to draw slope
mul = 0.5
triangle = np.array([[tri1,mul*tri1],[tri2,mul*tri2], [tri2,mul*tri1], [tri1, mul*tri1]])

print(triangle[:,0], triangle[:,1])
fig,ax1 = plt.subplots(figsize=(7,6))
ax2 = inset_axes(ax1, width="38%", height="38%", loc=2, bbox_to_anchor=(0.07,-0.02,1.0,1), bbox_transform=ax1.transAxes)

for ii in range(len(kappa_tildes)): # loop over \kappa
    ax1.loglog(init_cont_sin[ii,:].reshape(-1)-1, crit_strain_sin[ii,:].reshape(-1), label = rf'${{\tilde{{\kappa}}}} = 10^{{{np.log10(kappa_tildes[ii]):.0f}}}$', ls = 'None', marker = 's', fillstyle = 'none', c=c[ii], markeredgewidth=0.5, markersize = 7)
    ax1.loglog(init_cont_tri[ii,:].reshape(-1)-1, crit_strain_tri[ii,:].reshape(-1), ls = 'None', marker = '^', fillstyle = 'none', c = c[ii], markeredgewidth=0.5, markersize = 7)
    ax1.loglog(init_cont_rand[ii,:].reshape(-1)-1, crit_strain_rand[ii,:].reshape(-1), ls = 'None', marker = 'h', fillstyle = 'none', c = c[ii], markeredgewidth=0.5, markersize = 7)


    ax2.plot(init_cont_sin[ii,:].reshape(-1), crit_strain_sin[ii,:].reshape(-1), ls = 'None', marker = 's', fillstyle = 'none', c=c[ii], markeredgewidth=0.5, markersize = 7)
    ax2.plot(init_cont_tri[ii,:].reshape(-1), crit_strain_tri[ii,:].reshape(-1), ls = 'None', marker = '^', fillstyle = 'none', c = c[ii], markeredgewidth=0.5, markersize = 7)
    ax2.plot(init_cont_rand[ii,:].reshape(-1), crit_strain_rand[ii,:].reshape(-1), ls = 'None', marker = 'h', fillstyle = 'none', c = c[ii], markeredgewidth=0.5, markersize = 7)

ax1.plot(triangle[:,0], triangle[:,1], c='k', ls = '--', linewidth = 1)

ax1.axis('equal')
# ax1.text(tri1+0.1, tri2-0.1, '1', fontsize = 12)
# ax1.text(tri1+0.12, tri2-0.467, '1', fontsize = 12)
ax1.set_xlabel(r'$\sfrac{L_{C_0}}{L}-1$')
ax1.set_ylabel(r'$\varepsilon_{crit}$')

ax2.set_xlabel(r'$\sfrac{L_{C_0}}{L}$')
ax2.set_ylabel(r'$\varepsilon_{crit}$')
ax2.grid('on')

# Custom Legend
color_elements = [Patch(facecolor=c[0],
                         label=r'$10^{-3} $'),
                   Patch(facecolor=c[1],
                         label=r'$10^{-4} $'),
                   Patch(facecolor=c[2],
                         label=r'$10^{-5} $'),
                    Patch(facecolor=c[3],
                         label=r'$10^{-6} $')]

symbol_elements = [Line2D([0], [0], marker='s', color = 'None', markeredgecolor = 'k', label='Sinusoidal   ',
                          markerfacecolor='None', markersize=7),
                    Line2D([0], [0], marker='^', color = 'None', markeredgecolor = 'k', label='Triangular     ',
                          markerfacecolor='None', markersize=7),
                    Line2D([0], [0], marker='h', color = 'None', markeredgecolor = 'k', label='Random   ',
                          markerfacecolor='None', markersize=7)]

colors=ax1.legend(handles=color_elements,loc='lower right', title = r'Dimensionless Bending Length $\tilde{\kappa}$', ncol = 4, columnspacing=1.0)
ax1.legend(handles=symbol_elements,loc='lower right', title = 'Single Chain Type', bbox_to_anchor=(1.0, 0.125), ncol = 3, columnspacing=0.42)
ax1.add_artist(colors)
ax1.grid('on')
plt.tight_layout()
#plt.savefig('phase_2D.pdf')

print(crit_strain_rand[ii,:].reshape(-1)[0])
plt.savefig('phase_2D.png')

print(init_cont_rand.shape)
plt.show()