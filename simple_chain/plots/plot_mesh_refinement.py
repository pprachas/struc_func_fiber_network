import numpy as np
import matplotlib.pyplot as plt


plt.style.use('jeff_style.mplstyle')
#---------Triangular chain-----------#

crit_strain = np.loadtxt('../triangular_chain/mesh_refinement/crit_strain.txt')
num_ele = 3+np.arange(0,len(crit_strain))

# critical strain
plt.figure(figsize=(3,3))
plt.plot(num_ele,crit_strain, marker = 's') 
plt.xticks(range(3,len(crit_strain)+3))
plt.title('Triangular Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'$\varepsilon_{crit}$ ')
plt.tight_layout()
plt.savefig('triangular_mesh_eps.pdf')

# difference in critical strain
num_ele = 4+np.arange(0,len(crit_strain)-1)

plt.figure(figsize=(3,3))
plt.plot(num_ele,np.diff(crit_strain)/crit_strain[:-1]*100, marker = 's') 
plt.xticks(range(4,len(crit_strain)+3))
plt.title('Triangular Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'Percent Change in $\varepsilon_{crit}$ ')

plt.tight_layout()
plt.savefig('triangular_mesh_diff.pdf')

#---------Sinusoidal chain-----------#

crit_strain = np.loadtxt('../sinusoidal_chain/mesh_refinement/crit_strain.txt')
num_ele = 2*(2+np.arange(0,len(crit_strain)))-1

# critical strain
plt.figure(figsize=(3,3))
plt.plot(num_ele,crit_strain, marker = 's') 
plt.xticks(range(num_ele[0], num_ele[-1]+1))
plt.title('Sinusoidal Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'$\varepsilon_{crit}$ ')
plt.tight_layout()
plt.savefig('sinusoidal_mesh_eps.pdf')

# difference in critical strain
num_ele = 2*(3+np.arange(0,len(crit_strain)-1))-1

plt.figure(figsize=(3,3))
plt.plot(num_ele,np.diff(crit_strain)/crit_strain[:-1]*100, marker = 's') 
plt.xticks(range(num_ele[0], num_ele[-1]+1))
plt.title('Sinusoidal Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'Percent Change in $\varepsilon_{crit}$ ')

plt.tight_layout()
plt.savefig('sinusoidal_mesh_diff.pdf')


#---------Random chain-----------#

crit_strain = np.loadtxt('../random_chain/mesh_refinement/crit_strain.txt')
num_ele = 3+np.arange(0,len(crit_strain))

# critical strain
plt.figure(figsize=(3,3))
plt.plot(num_ele,crit_strain, marker = 's') 
plt.xticks(range(num_ele[0], num_ele[-1]+1))
plt.title('Random Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'$\varepsilon_{crit}$ ')
plt.tight_layout()
plt.savefig('random_mesh_eps.pdf')

# difference in critical strain
num_ele = 4+np.arange(0,len(crit_strain)-1)
plt.figure(figsize=(3,3))
plt.plot(num_ele,np.diff(crit_strain)/crit_strain[:-1]*100, marker = 's') 
plt.xticks(range(num_ele[0], num_ele[-1]+1))
plt.title('Random Chain')
plt.xlabel('Number of Elements per Fiber')
plt.ylabel(r'Percent Change in $\varepsilon_{crit}$ ')

plt.tight_layout()
plt.savefig('random_mesh_diff.pdf')

plt.show()