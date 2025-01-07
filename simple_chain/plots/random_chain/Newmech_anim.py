import numpy as np
import matplotlib.pyplot as plt
from dolfin import *
#--------Import functions---------------------#
import sys
sys.path.append('../../../utils')

import fea 

parameters['reorder_dofs_serial'] = False

plt.style.use('jeff_style.mplstyle')

L_0 = 10000 # total length of fiber
w = L_0/10 # width
num_link = 30 #amplitude
E = 1e6 # Stiffness moduli
num_run = 4
kappa_tilde = 1e-6

mesh_name=f'../../random_chain/mesh/w{int(w)}/n{num_link}/random_chain{num_run}.xdmf'
f_prefix_name = f'../../random_chain/fea_run/w{int(w)}/n{num_link}/kappa_tilde{kappa_tilde}/seed{int(num_run)}'

mesh, mf, dx = fea.read_mesh(mesh_name)
V,v,v_,dv,u,u_,theta,theta_ = fea.beam_function_space(mesh)

strain = np.loadtxt(f'{f_prefix_name}/disp.txt')/L_0


# Get dof coordinates:
x_dofs = V.sub(0).sub(0).dofmap().dofs()
y_dofs = V.sub(0).sub(1).dofmap().dofs()
theta_dofs = V.sub(1).dofmap().dofs()
dofs = V.tabulate_dof_coordinates()
dof_coords = dofs.reshape((-1, 2))

x_nodal_coord = dof_coords[x_dofs][:,0]
y_nodal_coord = dof_coords[y_dofs][:,1]

v0 = np.loadtxt(f'../../random_chain/fea_run/w{int(w)}/n{int(num_link)}/kappa_tilde{kappa_tilde}/seed{num_run}/v.txt')
# u1 = np.loadtxt('data/sinusoidal_1e-04/u.txt')
# u2 = np.loadtxt('data/sinusoidal_1e-05/u.txt')
# u3 = np.loadtxt('data/sinusoidal_1e-06/u.txt')

idx = np.argsort(y_nodal_coord)

# Initial shape
plt.figure(figsize=(2,7))
plt.plot(x_nodal_coord[idx], y_nodal_coord[idx], c = (0.0,0.0,0.0), lw = 6)

# plt.title(rf'$\varepsilon_{{applied}} = {strain[0]:.3f}$', fontsize=25)
plt.axis('equal')
plt.axis('off')
plt.ylim([0,np.max(strain)*L_0+L_0])

plt.savefig('deform/deform0000.png')
plt.savefig('init.pdf')

jj = 0
for ii in range(len(v0)):
    if ii%10 == 0:
        jj+=1
        # Plot displacement field
        disp_x0 = x_nodal_coord + v0[ii][x_dofs]
        disp_y0 = y_nodal_coord + v0[ii][y_dofs]

        idx = np.argsort(disp_y0)

        plt.figure(figsize=(2,7))
        plt.plot(disp_x0[idx], disp_y0[idx], c = (0.0,0.0,0.0), lw = 6)

        plt.axis('equal')
        plt.axis('off')
        plt.ylim([0,np.max(strain)*L_0+L_0])
        #plt.title(rf'$\varepsilon_{{applied}} = {strain[ii]:.3f}$', fontsize=25)
        plt.savefig(f'deform/deform{jj:04d}.png')
        plt.close()

crit_strain = 0.12074268431548262

idx = np.min(np.abs(strain-crit_strain))

print(idx)