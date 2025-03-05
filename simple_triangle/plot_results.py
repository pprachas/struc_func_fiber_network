import numpy as np 
import matplotlib.pyplot as plt

L = 10000
n = 500
seed = 1

crit_strain_all = []
init_contour = []
kappa_tilde = 1e-6

max_res_all =  ['0.1', '0.01']
plot_coef = [1.0,1.0,]
magnitude = [-2,-3]

plt.style.use('jeff_style.mplstyle')

colors = plt.cm.Reds(np.linspace(1,0,10))
colors = [colors[6], colors[3]]

plt.figure(figsize=(3,3))
plt.title('Force-Displacement Curve')
for count,max_res in enumerate(max_res_all):
    f_name = f'plots/sensitivity/{max_res}'
    force = np.loadtxt(f'{f_name}/force.txt')
    strain = np.loadtxt(f'{f_name}/disp.txt')/L

    plt.plot(strain,force, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c = colors[count], lw=1.5)
#------Plot solution no regularization-------#
f_name = f'plots/sensitivity/no_reg'
force = np.loadtxt(f'{f_name}/force.txt')
strain = np.loadtxt(f'{f_name}/disp.txt')/L
plt.plot(strain,force, label = r'no regularization', c = 'k', ls = ':', lw=2.0)

plt.xlabel('Applied Strain')
plt.ylabel('Reaction Force')
plt.legend(title = r'$||\bm{R}_{eq}||/{F_{react}}$')
plt.tight_layout()
plt.savefig('tri_force.pdf')

plt.figure(figsize=(3,3))
plt.title('Equilibrium Residual')
for count,max_res in enumerate(max_res_all):
    f_name = f'plots/sensitivity/{max_res}'
    residual = np.loadtxt(f'{f_name}/residual.txt')
    disp = np.loadtxt(f'{f_name}/disp.txt')/L
    plt.semilogy(disp,residual, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c = colors[count], lw=1.5)

#------Plot solution no regularization-------#
f_name = f'plots/sensitivity/no_reg'
residual = np.loadtxt(f'{f_name}/residual.txt')
strain = np.loadtxt(f'{f_name}/disp.txt')/L
plt.plot(strain,residual, label = r'no regularization', c = 'k', ls = ':', lw=2.0)
plt.xlabel('Applied Strain')
plt.ylabel(r'$||\bm{R}_{eq}||$')
plt.tight_layout()
plt.savefig('tri_sensitivity_residual.pdf')

plt.figure(figsize=(3,3))
plt.title('Stretch Energy')
for count,max_res in enumerate(max_res_all):
    f_name = f'plots/sensitivity/{max_res}'
    stretch_energy = np.loadtxt(f'{f_name}/stretch_total.txt')
    bend_energy = np.loadtxt(f'{f_name}/bend_total.txt')
    shear_energy = np.loadtxt(f'{f_name}/shear_total.txt')
    disp = np.loadtxt(f'{f_name}/disp.txt')/L
    plt.plot(disp,stretch_energy+bend_energy+shear_energy, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c=colors[count], lw=1.5)

#------Plot solution no regularization-------#
f_name = f'plots/sensitivity/no_reg'
force = np.loadtxt(f'{f_name}/force.txt')
stretch_energy = np.loadtxt(f'{f_name}/stretch_total.txt')
bend_energy = np.loadtxt(f'{f_name}/bend_total.txt')
shear_energy = np.loadtxt(f'{f_name}/shear_total.txt')
strain = np.loadtxt(f'{f_name}/disp.txt')/L
plt.plot(strain,stretch_energy+bend_energy+shear_energy, label = r'no regularization', c = 'k', ls = ':', lw=2.0)
plt.xlabel('Applied Strain')
plt.ylabel('Strain Energy')
plt.tight_layout()
plt.savefig('tri_strain_energy.pdf')

plt.show()