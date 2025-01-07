import numpy as np 
import matplotlib.pyplot as plt 

plt.style.use('jeff_style.mplstyle')
seed = 0
n = 500


#-------Plot strain value---------------#
plt.figure(figsize=(4,3))

crit_strain = np.loadtxt(f'../mesh_refinement/damped/mesh/n{n}/crit_strain{seed}.txt')

nodes = 3+np.arange(0,len(crit_strain))
plt.plot(nodes,crit_strain, c = 'k', marker='s')

plt.xticks(range(3,9))
plt.title(r'$\varepsilon_{crit}$ vs. Mesh Density')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'$\varepsilon_{crit}$ ')
plt.tight_layout()

plt.savefig('mesh_refinement_eps.pdf')


#------Plot difference in strain--------#
plt.figure(figsize=(4,3))

crit_strain = np.loadtxt(f'../mesh_refinement/damped/mesh/n{n}/crit_strain{seed}.txt')

nodes = 4+np.arange(0,len(crit_strain)-1)
plt.plot(nodes,np.diff(crit_strain)/crit_strain[:-1]*100, c = 'k', marker = 's')

plt.xticks(range(4,9))
plt.title(r'Change of $\varepsilon_{crit}$ vs. Mesh Density')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'Percent Change in $\varepsilon_{crit}$ ')

plt.tight_layout()
plt.savefig('mesh_refinement_diff.pdf')

plt.show()

