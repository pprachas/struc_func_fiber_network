from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import networkx as nx
from ufl import diag, Jacobian
import pathlib

#--------Import functions---------------------#
import sys
sys.path.append('../../utils')

import fea 
import graph

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters['reorder_dofs_serial'] = False



L = 10000 # total length of fiber
H = L
W = L
mesh_name=f'mesh/mikado.xdmf'
n_paths = 1 # number of shortest paths 

kappa_tilde = [1e-3,1e-4,1e-5,1e-6][int(sys.argv[1])-1]
r = 2*201.72070000308315*np.sqrt(kappa_tilde)

mesh, mf, dx = fea.read_mesh(mesh_name)

V,v,v_,dv,u,u_,theta,theta_ = fea.beam_function_space(mesh)

# define tangent vectors
R0, g01 = fea.compute_tangents(mesh)

Rot = fea.rotation_matrix(theta)

# Boundary Conditions
def top(x, on_boundary):
    return near(x[1],H,1e-6)

def bot(x, on_boundary):
    return near(x[1],0,1e-6)

def left(x,on_boundary):
    return near(x[0],0,1e-6)

def right(x,on_boundary):
    return near(x[0], W, 1e-6)


# finding the node number of one of the bottom to fix
x_dofs = V.sub(0).sub(0).dofmap().dofs()
dof_coords = V.tabulate_dof_coordinates().reshape((-1, 2))

bot_node_coords = []
for ii in x_dofs:
    if abs(dof_coords[ii,1]) <= 1e-6:
            bot_node_coords.append(dof_coords[ii,0])
        
bot_node_coords = np.array(bot_node_coords)

def bot_x(x, on_boundary):
    return near(x[0],np.min(bot_node_coords)) and near(x[1],0,1e-6)

# Mark subdomains for Dirichlet BCs
BC_bot_y = DirichletBC(V.sub(0).sub(1), Constant(0.0), bot) # fixed displacement
#BC_bot_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), bot) # fix x displacement
BC_left_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), left)
BC_right_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), right)
#BC_top_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # fixed x displacement on top

disp = Expression("t",t=0,degree=1)
BC_top_y = DirichletBC(V.sub(0).sub(1), disp, top)

bc = [BC_bot_y, BC_top_y, BC_left_x, BC_right_x]

# get strain measures:
defo, curv = fea.strains(R0,Rot,g01,u,theta)

# Geometrical properties (circle cross section; 0 Poisson's ratio)
r = Constant(r)
S = Constant(pi*r**2)
I = Constant(pi*r**4/4)
kappa = Constant(6/7)

# Stiffness moduli
E = Constant(1e6)
G = Constant(E/2)

# Beam stiffness parameters
ES = Constant(E*S)
GS = Constant(kappa*G*S)
EI = Constant(E*I)

# Constitutive Equations
C_N = diag(as_vector([ES, GS]))

# Applied Load:
F_max = Constant((0.0,0.0))
M_max = Constant(0.0)

elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

F_int = derivative(elastic_energy, v, v_)
F_ext = (M_max*theta_ + dot(F_max, u_)) * ds
residual = F_int - F_ext
tangent_form = derivative(residual, v, dv)

#-------------------Define Newton Solver----------------------------------#
problem = NonlinearVariationalProblem(residual, v, bc, tangent_form)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters["newton_solver"]
prm["linear_solver"] = "default"
prm["absolute_tolerance"] = 1e-9*E
prm["relative_tolerance"] = 1e-9#1e-8
prm['relaxation_parameter'] = 1.00
bend_total_all = []
stretch_total_all = []
shear_total_all = []

bend_total_frac = []
stretch_total_frac = []
shear_total_frac = []
num_fibers = int(mf.array().max()+1) # number of fibers in the network

# Energy functions to integrate
stretch = 0.5 * defo[0]*dot(C_N, defo)[0]
shear = 0.5 * defo[1]*dot(C_N, defo)[1]
bend = 0.5*EI*curv**2

pathlib.Path(f'./paraview/{kappa_tilde}').mkdir(parents=True, exist_ok=True)

out_file = XDMFFile(f'./paraview/{kappa_tilde}/mikado.xdmf')
out_file.parameters["functions_share_mesh"] = True
out_file.parameters["flush_output"] = True

reaction_force = []
bcRy = DirichletBC(V.sub(0).sub(1), Constant(1.), top)

dofs = V.tabulate_dof_coordinates()
dof_coords = dofs.reshape((-1, 2))

link_displacements_x = []
link_displacements_y = []
cont_length = []
shortest_path_energy_all = []
for ii in range(n_paths):
    cont_length.append([])
    shortest_path_energy_all.append([])

inc = []
step_size = L/10000.
step_size_init = step_size

# Function spaces for strains
T_translational = FunctionSpace(mesh, 'DG', 1)
axial = Function(T_translational, name = 'axial_strain')
shear = Function(T_translational, name = 'shear strain')

T_rotational = FunctionSpace(mesh, 'DG', 0)
rot = Function(T_rotational, name = 'rotational_strain')

T_energy = FunctionSpace(mesh, 'DG', 0)
energy = Function(T_energy, name = 'rotational_strain')

ii = 0

cond = []
# Iterate over load steps with Newton Solver -- use smaller steps if doesn't converge
while True:
    try:
        step_size = min(2*step_size,step_size_init)
        disp.t += step_size
        iteration,converged = solver.solve()
        v_converged = v.vector()[:] # save converged solution
        converged_prev = True
    except:
        print('solver fails!')
        disp.t -= step_size # go back to previously converged step
        step_size = min(step_size/2,step_size/1000) # use smaller step  
        v.vector()[:] = v_converged # use previously converged solution as initial guess
        converged_prev = False
        continue

    if converged:
        print(disp.t/L)
        #---------------------Compute energy------------------------------------------#
        stretch_total, bend_total, shear_total = fea.compute_energy(num_fibers, stretch, bend, shear, dx)
        total_energy = stretch_total + bend_total + shear_total
        bend_total_all.append(bend_total)
        stretch_total_all.append(stretch_total)
        shear_total_all.append(shear_total)

        bend_total_frac.append(bend_total/total_energy)
        stretch_total_frac.append(stretch_total/total_energy)
        shear_total_frac.append(shear_total/total_energy)

        reaction_force.append(fea.compute_rxn_force(V,bcRy, residual))
        
        inc.append(disp.t)

        epsilon = defo

        epsilon_axial = project(epsilon[0],T_translational)
        axial.assign(epsilon_axial)

        epsilon_shear = project(epsilon[1],T_translational)
        shear.assign(epsilon_shear)

        chi_z = project(curv,T_rotational)
        rot.assign(chi_z)

        out_file.write(v.sub(0,True),ii)
        out_file.write(axial,ii)
        out_file.write(shear,ii)
        out_file.write(rot,ii)
        
        ii+=1

        # # Assemble stiffness matrix
        # K = PETScMatrix()
        # f = PETScVector()

        # K=assemble(tangent_form, tensor = K)

        # for bci in bc:
        #     bci.apply(K)

        # eigenSolver = SLEPcEigenSolver(K)
        
        # eigenSolver.parameters['spectral_transform'] = 'shift-and-invert'
        # eigenSolver.parameters['spectral_shift'] = 0.0

        # eigenSolver.solve(1)
        # eigen_min=eigenSolver.get_eigenvalue(0)[0]
        # print(eigenSolver.get_eigenvalue(0)[0])

        # eigenSolver = SLEPcEigenSolver(K)
        # eigenSolver.parameters["spectrum"]="largest magnitude"
        # eigenSolver.solve(1)
        # eigen_max=eigenSolver.get_eigenvalue(0)[0]
        # print(eigen_max/eigen_min)

        # cond.append(eigen_max/eigen_min)

        # plt.figure()
        # plt.semilogy(np.array(inc)/L, cond)
        # plt.xlabel('Applied Strain')
        # plt.ylabel('Condition Number')
        # plt.tight_layout()
        # plt.savefig('cond.png')
        # plt.close()

        # for count,path in enumerate(shortest_path):
        #     shortest_path_length = 0
        #     shortest_path_energy = 0
        #     shortest_path_stretch_energy = 0

        #     for kk in path:
        #         shortest_path_length += assemble(sqrt(inner(g01+fea.tgrad(u,g01),g01+fea.tgrad(u,g01)))*dx(kk))
        #         shortest_path_energy += assemble(0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx(kk))
        #         shortest_path_stretch_energy += assemble(0.5 * (dot(defo, dot(C_N, defo))*dx(kk)))
        #     cont_length[count].append(shortest_path_length)
        #     shortest_path_energy_all[count].append(shortest_path_stretch_energy/shortest_path_energy)
        
        if disp.t/L >= 0.02:
            break
        # if stretch_total > bend_total :
        #     break

out_file.close()

# print('Stretch Energy:', stretch_total)
# print('Bend Energy:', bend_total)
# print('Energy Difference:', stretch_total-bend_total)

# eigen_min = np.array(eigen_min)
# eigen_max = np.array(eigen_max)
# cond = eigen_max/eigen_min


# plt.figure()
# plt.title('Energy vs apply strain')
# plt.xlabel('Applied Strain')
# plt.ylabel('Energy')
# plt.plot(np.array(inc)/L,bend_total_all, label = 'Bending')
# plt.plot(np.array(inc)/L,stretch_total_all, label = 'Stretching')
# plt.plot(np.array(inc)/L,shear_total_all, label = 'Shear')
# plt.legend()

# plt.figure()
# plt.title('Percent energy vs apply strain')
# plt.xlabel('Applied Strain')
# plt.ylabel('Percent Energy')
# plt.plot(np.array(inc)/L,bend_total_frac, label = 'Bending', marker = '.')
# plt.plot(np.array(inc)/L,stretch_total_frac, label = 'Stretching', marker = '.')
# plt.plot(np.array(inc)/L,shear_total_frac, label = 'Shear', marker = '.')
# #plt.axvline(x=0.2882455862326643,color="black", linestyle=':', label = 'Bending-stretching transition')
# plt.legend()

pathlib.Path(f'./plots').mkdir(parents=True, exist_ok=True)

plt.figure()
plt.title('Stress vs Strain')
plt.xlabel('Applied Strain')
plt.ylabel('Stress')
plt.plot(np.array(inc)/L,np.array(reaction_force)/float(ES), lw = 5, c = 'r')
plt.savefig(f'plots/fd{kappa_tilde}.png')

# plt.figure()
# plt.title('energy difference')
# plt.xlabel('Applied Strain')
# plt.ylabel('Stress')
# plt.plot(np.array(inc),np.array(bend_total_all)-np.array(stretch_total_all))

# plt.savefig(f'plots/n{n}/seed{seed}/fd.png')

# plt.figure()
# plt.title('straightness vs applied strain')
# plt.ylabel('Straightness')
# plt.xlabel('Applied Strain')
# for ii in range(n_paths):
#     plt.plot((np.array(inc))/L,((np.array(inc)+L)/np.array(cont_length[ii])), label = str(ii+1))

# plt.legend()
# plt.savefig(f'plots/n{n}/seed{seed}/straightness.png')

# plt.figure()
# plt.title('Contour length vs applied strain')
# plt.ylabel('Contour length')
# plt.xlabel('Applied Strain')
# for ii in range(n_paths):
#     plt.plot((np.array(inc))/L,np.array(cont_length[ii]), label = str(ii+1), c = 'k', lw = 5)

# plt.legend()
# plt.savefig(f'plots/n{n}/seed{seed}/contour.png')

# plt.figure()
# plt.title('Shortest path energy percent vs applied strain')
# plt.ylabel('Percent Stretch Energy')
# plt.xlabel('Applied Strain')

# for ii in range(n_paths):
#     plt.plot((np.array(inc))/L,np.array(shortest_path_energy_all[ii]), label = 'path' + str(ii+1), marker = 'o', fillstyle='none', ls = 'None')
# plt.plot((np.array(inc))/L, stretch_total_frac, label = 'Total Fiber Network')

# plt.legend()
# plt.savefig(f'plots/n{n}/seed{seed}/energy.png')

# plt.figure()
# plt.semilogy((np.array(inc))/L,cond)

# Save info for seperate plotting
# pathlib.Path(f'./plots/n{n}/seed{seed}/data').mkdir(parents=True, exist_ok=True)

# np.savetxt(f'./plots/n{n}/seed{seed}/data/strain.txt', np.array(inc)/L)
# np.savetxt(f'./plots/n{n}/seed{seed}/data/cont_length.txt', np.array(cont_length)/L)
# np.savetxt(f'./plots/n{n}/seed{seed}/data/force.txt', np.array(reaction_force)/float(ES))
plt.show()