import numpy as np 
import matplotlib.pyplot as plt

L = 10000
n = 500
seed = 1

crit_strain_all = []
init_contour = []
kappa_tilde = 1e-6

max_res_all = [5e-3,1e-3, 5e-4, 1e-4, 5e-5,1e-5]
plot_coef = [0.5,1.0,0.5,1.0,0.5,1.0]
magnitude = [-2,-3,-3,-4,-4,-5]

plt.style.use('jeff_style.mplstyle')

colors = plt.cm.Reds(np.linspace(1,0,10))


plt.figure(figsize=(3,3))
plt.title('Force-Displacement Curve')
for count,max_res in enumerate(max_res_all):
    f_name = f'sensitivity/damping/max_res{max_res}'
    force = np.loadtxt(f'{f_name}/force.txt')
    strain = np.loadtxt(f'{f_name}/disp.txt')/L

    plt.plot(strain,force, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c = colors[count], lw=1.5)
plt.xlabel('Applied Displacement')
plt.ylabel('Reaction Force')
plt.legend(title = r'$||\bm{R}_{eq}||/{F_{react}}$')
plt.tight_layout()
plt.savefig('sensitivity/force.pdf')

plt.figure(figsize=(3,3))
plt.title('Equilibrium Residual')
for count,max_res in enumerate(max_res_all):
    f_name = f'sensitivity/damping/max_res{max_res}'
    residual = np.loadtxt(f'{f_name}/residual.txt')
    disp = np.loadtxt(f'{f_name}/disp.txt')/L
    plt.semilogy(disp,residual, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c = colors[count], lw=1.5)
plt.xlabel('Displacement')
plt.ylabel(r'$||\bm{R}_{eq}||$')
plt.legend(title = r'$({||\bm{R}_{eq}||}/{F_{react}})$')
plt.tight_layout()
plt.savefig('sensitivity_residual.pdf')

plt.figure(figsize=(3,3))
plt.title('Stretch Energy')
for count,max_res in enumerate(max_res_all):
    f_name = f'sensitivity/damping/max_res{max_res}'
    stretch_energy = np.loadtxt(f'{f_name}/stretch_total.txt')
    bend_energy = np.loadtxt(f'{f_name}/bend_total.txt')
    shear_energy = np.loadtxt(f'{f_name}/shear_total.txt')
    disp = np.loadtxt(f'{f_name}/disp.txt')/L
    plt.plot(disp,stretch_energy+bend_energy+shear_energy, label = rf'${plot_coef[count]} \times 10^{{{magnitude[count]}}}$', c=colors[count], lw=1.5)
plt.xlabel('Displacement')
plt.ylabel('Strain Energy')
plt.legend(title = r'${||\bm{R}_{eq}||}/{F_{react}}$')
plt.tight_layout()
plt.savefig('sensitivity/strain_energy.pdf')

plt.figure(figsize=(3,3))
for max_res in max_res_all:
    f_name = f'sensitivity/damping/max_res{max_res}'
    crit_strain = np.loadtxt(f'{f_name}/crit_strain.txt')
    plt.semilogx(max_res**(-1), crit_strain, marker = 's', c = 'k', markersize = 5.0)
plt.xlabel(r'$({||\bm{R}_{eq}||}/{F_{react}})^{-1}$')
plt.ylabel(r'$\varepsilon_{crit}$')
plt.title(r'$\varepsilon_{crit}$ Convergence Study')
plt.tight_layout()
plt.savefig('sensitivity/critical_strain.pdf')

plt.show()