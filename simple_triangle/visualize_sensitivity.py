import numpy as np
import pathlib
import matplotlib.pyplot as plt
import networkx as nx
import sys
from dolfin import *

#--------Import functions---------------------#
import sys
sys.path.append('../utils')

import fea_HR as fea

plt.style.use('jeff_style.mplstyle')

pathlib.Path(f'./sensitivity/deformed_network').mkdir(parents=True, exist_ok=True)

L = 10000

max_res_all = ['no_reg', '0.1', '0.01']

parameters["form_compiler"]["cpp_optimize"] = True
parameters['reorder_dofs_serial'] = False

colors = plt.cm.Reds(np.linspace(1,0,10)) # same color theme as force displacement curve
colors = ['k', colors[6], colors[3]]
markersize = [2.0,0.75,0.75]

plt.figure(figsize=(4,3))
for ii,max_res in enumerate(max_res_all):
    # Name of mesh
    mesh_name=f'mesh/simple_triangle.xdmf'
    f_name = f'plots/sensitivity/{max_res}'

    v_all = np.loadtxt(f'{f_name}/v.txt')

    mesh, mf, dx = fea.read_mesh(mesh_name)


    V, v, v_, dv, u, u_,du,dtheta,dn,dm, theta, theta_, n, n_, m, m_ = fea.beam_function_space(mesh)

    v_soln = v_all[-1]

    x_dofs = V.sub(0).sub(0).dofmap().dofs()
    y_dofs = V.sub(0).sub(1).dofmap().dofs()
    theta_dofs = V.sub(1).dofmap().dofs()
    dofs = V.tabulate_dof_coordinates()
    dof_coords = dofs.reshape((-1, 2))

    # notal values
    x_nodal_coord = dof_coords[x_dofs][:,0]
    y_nodal_coord = dof_coords[y_dofs][:,1]

    # Plot displacement field
    disp_x = x_nodal_coord + v_soln[x_dofs]
    disp_y = y_nodal_coord + v_soln[y_dofs]


    plt.plot(disp_x,disp_y, markersize = markersize[ii], marker = 's', ls = 'none', c = colors[ii])


plt.axis('off')
plt.axis('equal')
plt.tight_layout()

plt.savefig(f'tri_compare_sensitivity.pdf')


plt.show()
    