import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

color=[(0.2,0.2,0.2), (0.5,0.5,0.5), (0.9,0.5,0.5), (0.5,0.0,0.0)]

f_name = './MAPE_data/'

error_stretch = np.loadtxt(f'{f_name}area_stretch.txt')
error_stretch_bend = np.loadtxt(f'{f_name}area_stretch_bend.txt')

# error_stretch = np.split(error_stretch,4)
# error_stretch_bend = np.split(error_stretch_bend,4)
plt.style.use('jeff_style.mplstyle')

fig,ax = plt.subplots(1,1, figsize=(7,4))
plt.title('Size vs Area Between Curves')
plt.xlabel(r'Domain Size $\sfrac{L_0}{\ell_c}$')
plt.ylabel(r'Area Between Curves')
for ii in range(4):
      idx = np.where(error_stretch[:,-1] == ii)
      plt.plot(error_stretch[:,1][idx],error_stretch[:,2][idx], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
      plt.plot(error_stretch_bend[:,1][idx],error_stretch_bend[:,2][idx], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)
      print(np.mean(error_stretch[:,2][idx]))
plt.tight_layout()
plt.savefig('ABC_size.pdf')

fig,ax = plt.subplots(1,1, figsize=(7,4))
plt.title('Approximation Angle vs Area Between Curves')
plt.xlabel(r'Mean Approximation Angle $\langle \theta_{C_0} \rangle$')
plt.ylabel(r'Area Between Curves')
for ii in range(4):
      idx = np.where(error_stretch[:,-1] == ii)
      plt.plot(error_stretch[:,0][idx],error_stretch[:,2][idx], ls = 'None', fillstyle = 'none', marker = 'o', c = color[ii], markersize=7, markeredgewidth = 0.5)
      plt.plot(error_stretch_bend[:,0][idx],error_stretch_bend[:,2][idx], ls = 'None', fillstyle = 'none', marker='*', c=color[ii], markersize=7, markeredgewidth = 0.5)

plt.tight_layout()
plt.savefig('ABC_theta.pdf')

plt.show()