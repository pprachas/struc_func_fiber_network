import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

#---------Triangular Chain--------------#
f_name = './MAPE_data/triangular_chain/'

error_stretch = np.loadtxt(f'{f_name}area_stretch.txt')
error_stretch_bend = np.loadtxt(f'{f_name}area_stretch_bend.txt')

error_stretch = np.split(error_stretch,4)
error_stretch_bend = np.split(error_stretch_bend,4)
plt.style.use('jeff_style.mplstyle')

fig,ax = plt.subplots(1,1, figsize=(3,3))
plt.xlabel(r'Approximation angle $\theta_{C_0}$')
plt.ylabel(r'Area Between Curves')
plt.title('Triangular Chains')

xlim = (0.46,0.467)
ylim = (0.0,6e-5)

axins = ax.inset_axes([0.1,0.6,0.3,0.3], xlim = xlim,ylim=ylim)

for ii in range(4):
    ax.plot(error_stretch[ii][:,0],error_stretch[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    ax.plot(error_stretch_bend[ii][:,0],error_stretch_bend[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)

    idx = np.argmin(error_stretch[ii][:,0])

    inset_plot = axins.plot(error_stretch[ii][idx,0],error_stretch[ii][idx,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    axins.plot(error_stretch_bend[ii][idx,0],error_stretch_bend[ii][idx,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)

axins.tick_params(labelsize=8)

mark_inset(ax, axins, loc1=3, loc2=4, ls = ':', color = '0.7')

plt.tight_layout()
plt.savefig('triangular_MAPE.pdf')

#---------Sinusoidal Chain--------------#
f_name = './MAPE_data/sinusoidal_chain/'

error_stretch = np.loadtxt(f'{f_name}area_stretch.txt')
error_stretch_bend = np.loadtxt(f'{f_name}area_stretch_bend.txt')

error_stretch = np.split(error_stretch,4)
error_stretch_bend = np.split(error_stretch_bend,4)

fig,ax = plt.subplots(1,1, figsize=(3,3))
plt.xlabel(r'Approximation angle $\theta_{C_0}$')
plt.ylabel(r'Area Between Curves')
plt.title('Sinusoidal Chains')

xlim = (0.498,0.5007)
ylim = (-3e-6,4.5e-5)

axins = ax.inset_axes([0.1,0.6,0.3,0.3], xlim = xlim,ylim=ylim, yticks = [0,3e-5])

for ii in range(4):
    ax.plot(error_stretch[ii][:,0],error_stretch[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    ax.plot(error_stretch_bend[ii][:,0],error_stretch_bend[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)

    idx = np.argmin(error_stretch[ii][:,0])

    inset_plot = axins.plot(error_stretch[ii][idx,0],error_stretch[ii][idx,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    axins.plot(error_stretch_bend[ii][idx,0],error_stretch_bend[ii][idx,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)
axins.tick_params(labelsize=8)
mark_inset(ax, axins, loc1=3, loc2=4, ls = ':', color = '0.7')

plt.tight_layout()
plt.savefig('sinusoidal_MAPE.pdf')


#---------Random Chain--------------#
f_name = './MAPE_data/random_chain/'

error_stretch = np.loadtxt(f'{f_name}area_stretch.txt')
error_stretch_bend = np.loadtxt(f'{f_name}area_stretch_bend.txt')

error_stretch_split = np.split(error_stretch,4)
error_stretch_bend_split = np.split(error_stretch_bend,4)


fig,ax = plt.subplots(1,1, figsize=(3,3))
plt.xlabel(r'Approximation angle $\theta_{C_0}$')
plt.ylabel(r'Area Between Curves')
plt.title('Random Chains')

xlim = (0.111,0.14)
ylim = (-1e-6,1.8e-5)

axins = ax.inset_axes([0.1,0.6,0.3,0.3], xlim = xlim,ylim=ylim)

for ii in range(4):
    plt.plot(error_stretch_split[ii][:,0],error_stretch_split[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=6, markeredgewidth = 0.5)
    plt.plot(error_stretch_bend_split[ii][:,0],error_stretch_bend_split[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=6, markeredgewidth = 0.5)

    idx = np.argmin(error_stretch_split[ii][:,0])

    inset_plot = axins.plot(error_stretch_split[ii][:,0],error_stretch_split[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    axins.plot(error_stretch_bend_split[ii][:,0],error_stretch_bend_split[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)
axins.tick_params(labelsize=8)
mark_inset(ax, axins, loc1=3, loc2=4, ls = ':', color = '0.7')

plt.tight_layout()
plt.savefig('random_MAPE.pdf')
#---------Discrete Sin--------------#
f_name = './MAPE_data/discretized_sin/'

error_stretch = np.loadtxt(f'{f_name}area_stretch.txt')
error_stretch_bend = np.loadtxt(f'{f_name}area_stretch_bend.txt')

error_stretch = np.split(error_stretch,4)
error_stretch_bend = np.split(error_stretch_bend,4)



fig,ax = plt.subplots(1,1, figsize=(4,4))
plt.xlabel(r'Approximation angle $\theta_{C_0}$')
plt.ylabel(r'Area Between Curves')
plt.title('Discretized Sinusoids')

xlim = (0.4810,0.482)
ylim = (-1e-6,3.8e-5)

axins = ax.inset_axes([0.1,0.6,0.3,0.3], xlim = xlim,ylim=ylim)

for ii in range(4):
    plt.plot(error_stretch[ii][:,0],error_stretch[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=6, markeredgewidth = 0.5)
    plt.plot(error_stretch_bend[ii][:,0],error_stretch_bend[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=6, markeredgewidth = 0.5)

    inset_plot = axins.plot(error_stretch[ii][:,0],error_stretch[ii][:,1], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
    axins.plot(error_stretch_bend[ii][:,0],error_stretch_bend[ii][:,1], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)
axins.tick_params(labelsize=8)
mark_inset(ax, axins, loc1=3, loc2=4, ls = ':', color = '0.7')

plt.tight_layout()

plt.savefig('discretized_MAPE.pdf')

#-------Legend----------#
plt.figure()


line5 = Line2D([0], [0], label='FEA', color='k', lw = 3.0, fillstyle='none', ls='none')
line6 = Line2D([0], [0], label='axial',marker = 'o', c = 'k', lw = 0.5, fillstyle='none', ls='none')
line7 = Line2D([0], [0], label='axial + rotation',marker = '*', c = 'k', lw = 0.5, fillstyle='none', ls='none')

handles =[line5, line6, line7]

plt.legend(handles=handles, title = r'Model Type', ncol=3)
plt.savefig('model_legend2.pdf')

plt.show()
