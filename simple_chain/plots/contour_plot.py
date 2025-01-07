import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy import optimize
import matplotlib.ticker as mticker
from scipy.interpolate import griddata
from matplotlib.pyplot import contourf
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, ticker



# plot settings
plt.style.use('jeff_style.mplstyle')
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rc('font', size=12)


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

kappa_tilde = []
kappa_tilde_rand = []
# Load all sinusoidal and triangular data
for ii in a:


    # Sinusoidal
    crit_strain_sin.append(np.loadtxt(f'../sinusoidal_chain/phase_diagram/crit_strain/a{ii}.txt'))
    init_cont_sin.append(np.loadtxt(f'../sinusoidal_chain/phase_diagram/init_cont/a{ii}.txt')/L)

    # triangular
    crit_strain_tri.append(np.loadtxt(f'../triangular_chain/phase_diagram/crit_strain/a{ii}.txt'))
    init_cont_tri.append(np.loadtxt(f'../triangular_chain/phase_diagram/init_cont/a{ii}.txt')/L)

    kappa_tilde.append(kappa_tildes*4)



# Load all random data
for ii in w_r:
    for jj in n:
        ii = int(ii)
        jj = int(jj)
        crit_strain_rand.append(np.loadtxt(f'../random_chain/phase_diagram/crit_strain/w{ii}/n{jj}.txt'))
        init_cont_rand.append(np.loadtxt(f'../random_chain/phase_diagram/init_cont/w{ii}/n{jj}.txt')/L)
        kappa_tilde_rand.append(kappa_tildes*20)

crit_strain_sin = np.array(crit_strain_sin).reshape(-1)
init_cont_sin = np.array(init_cont_sin).reshape(-1)
kappa_tilde = np.array(kappa_tilde).reshape(-1)

crit_strain_tri = np.array(crit_strain_tri).reshape(-1)
init_cont_tri = np.array(init_cont_tri).reshape(-1)

crit_strain_rand = np.array(crit_strain_rand).reshape(-1)
init_cont_rand = np.array(init_cont_rand).reshape(-1)
kappa_tilde_rand = np.array(kappa_tilde_rand).reshape(-1)


init_cont_all = np.concatenate((init_cont_sin,init_cont_tri, init_cont_rand))
crit_strain_all = np.concatenate((crit_strain_sin,crit_strain_tri, crit_strain_rand))
kappa_tilde_all = np.concatenate((kappa_tilde,kappa_tilde,kappa_tilde_rand))


def log_tick_formatter(val, pos=None):
    return rf"$10^{{{int(val)}}}$"  # remove int() if you don't use MaxNLocator
    # return f"{10**val:.2e}"      # e-Notation

# ax.view_init(0,0,0)

x_bound = np.linspace(np.min(init_cont_all), np.max(init_cont_all), 100)
y_bound = np.logspace(-3, -6, 100)

x_plot,y_plot = np.meshgrid(x_bound,y_bound)

z_plot = griddata((init_cont_all.ravel(),kappa_tilde_all.ravel()), crit_strain_all, (x_plot,y_plot), method = 'linear')
# use nrearest data points for Nan (extrpolation) values
z_nn = griddata((init_cont_all.ravel(),kappa_tilde_all.ravel()), crit_strain_all, (x_plot,y_plot), method = 'nearest')

# idx_test = np.isnan(z_plot)
z_plot[np.isnan(z_plot)] = z_nn[np.isnan(z_plot)]

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111, projection='3d')

# ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
# ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
# ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

# ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
# ax.zaxis.set_major_locator(mticker.MaxNLocator(integer=True))


# plt.title('Sinusoidal Chains')
# ax.set_xlabel('Initial Waviness')
# ax.set_ylabel(r'$\tilde{\kappa}$')
# ax.set_zlabel('Critical Strain')
# ax.plot_surface(np.log10(x_plot-1),np.log10(y_plot), np.log10(z_plot))
# ax.scatter3D(np.log10(init_cont_sin-1),np.log10(kappa_tilde),np.log10(crit_strain_sin), marker = 's', s = 30, edgecolor = 'k', linewidth = 1,facecolor='none')

fig,ax =plt.subplots(figsize=(4,4))
ax.semilogy()
cs = ax.contourf(x_plot,y_plot, z_plot, cmap=cm.Reds)
cs2 = ax.contour(cs, levels=cs.levels[::2], colors='k')

ax.set_xlabel(r'Initial Tortuosity $\tau_0$')
ax.set_ylabel(r'Dimensionless Bending Length $\tilde{\kappa}$')


manual_loc = [(6.08,4.4e-5), (10.71,4.4e-5),(15.56,4e-5)]
cbar=fig.colorbar(cs, label=r'Critical Strain $\varepsilon_{crit}$')
ax.clabel(cs2, inline=True, manual=manual_loc)
cbar.add_lines(cs2)


plt.tight_layout()
plt.savefig('contour.pdf')

# plt.figure()
# plt.semilogy(init_cont_all.ravel(),kappa_tilde_all.ravel(), ls = 'none', marker = '.')
# plt.semilogy(x_plot[idx_test],y_plot[idx_test], ls = 'none', marker = 's', fillstyle = 'none')
plt.show()