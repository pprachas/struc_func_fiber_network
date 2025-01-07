from dolfin import *
import numpy as np
from ufl import Jacobian, diag, sign
import matplotlib.pyplot as plt
'''
This is the main script to fun FEniCS simulations of fiber networks with Simo-Reissner beams.
The network will be clamped (constrained displacement and rotation) at the top and bottom.
The top will have prescribed non-zero displacement.
'''

def read_mesh(f_name):
    '''
    Reads mesh and domain information

    Args:
        f_name: file name for mesh in XDMF format. The mesh should also contain domain information

    Returns:
        mesh: dolfin mesh 
        dx: marked domains for integration
    '''
    mesh = Mesh()

    mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()) # create class to store data out individual fibers

    with XDMFFile(f_name) as infile:
        infile.read(mesh) # read mesh 

    mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()) # create class to store data out individual fibers
    with XDMFFile(f_name) as infile:
        infile.read(mvc, "fibers") # read domain information

    mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    # Mark Domains
    dx = Measure("dx",domain=mesh, subdomain_data=mf)
    
    return mesh, mf, dx

def beam_function_space(mesh):
    '''
    Create mixed function space for beams

    Returns:
    V: function space
    v: displacement function space
    u: displacement trial function space
    u_: displacement test function space
    theta: angle trial function space
    theta_: angle test function space
    
    '''
    

    # Define function space
    Ue = VectorElement("CG", mesh.ufl_cell(), 2, dim=2) # displacement
    Te = FiniteElement("CG", mesh.ufl_cell(), 1) # rotation


    V = FunctionSpace(mesh, MixedElement([Ue, Te]))   

    v_ = TestFunction(V)
    u_, theta_ = split(v_)
    dv = TrialFunction(V)
    v = Function(V, name="Generalized displacement")
    u, theta = split(v)

    return V, v, v_, dv, u, u_, theta, theta_


def compute_tangents(mesh):
    '''
    Define tangent vectors and collect them into an initial rotation matrix
    Args:
        mesh: dolfin mesh

    Returns:
        R0: Initial rotation matrix
        g01: initial tangenet direction
    '''
    Jac = Jacobian(mesh)
    gdim = mesh.geometry().dim()
    Jac = as_vector([Jac[i, 0] for i in range(gdim)])
    g01 = Jac/sqrt(dot(Jac, Jac))
    g02 = as_vector([-g01[1],g01[0]])

    r01 = outer(g01,as_vector([1,0]))
    r02 = outer(g02, as_vector([0,1]))

    R0 = r02+r02


    return r01+r02, g01 

def rotation_matrix(theta): 
    '''
    Function to compute the current rotation matrix
    Args: 
        theta: current rotation function space
    Returns:
        current rotation matrix

    '''
    return as_tensor([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])


def define_bcs(L):
    '''
    used to define boundary conditions
    Args:
        L: domain size in the y-direction
    Returns:
        bc: Dirichlet BCs
    '''

def tgrad(u, g01):
    return dot(grad(u), g01)

def strains(R0,Rot,g01,u,theta):
    '''
    strain measures for beams (axial, shear, curvature)
    Args:
        u: displacement function
    Returns:
        defo: axial and shear strains
        curv: curvature strains
    '''
    # Strain measures
    defo = dot(R0.T,dot(Rot.T, g01 + tgrad(u, g01)) - g01)
    curv =  tgrad(theta, g01)
    return defo, curv

def postprocess_fiber(num_fibers, stretch, bend, shear, dx, translational = 0.0, rotational = 0.0, curvature = 0.0, compute_individual = False):
    '''
    Compute bending, stretch and shear energy of both the whole system and individual fibers
    Args:
        num_fibers: number of fibers
        stretch : stretching energy integrand
        bend: bending energy integrand
        shear: shearning energy integrand
        dx = integration domain seperated per fibers
        compute_individual: True if needed to compute energy individual fibers; False by default
    Returns:
        stretch_total: total stretch energy
        bend_total: total bend energy
        shear_total: total shear energy
        stretch_energy: stretch energy of each fiber
        bend_energy: bend energy of each fiber
        shear_energy: shear energy of each fiber
    '''
    stretch_total = assemble(stretch*dx)
    bend_total = assemble(bend*dx)
    shear_total = assemble(shear*dx)
    translational_total = assemble(translational*dx)
    rotational_total = assemble(rotational*dx)
    curvature_total = assemble(curvature*dx)

    stretch_energy = []
    bend_energy = []
    shear_energy = []
    translational_single = []
    rotational_single = []
    curvature_single = []

    if compute_individual:
        for ii in range(num_fibers):
            stretch_energy.append(assemble(stretch*dx(ii)))
            bend_energy.append(assemble(bend*dx(ii)))
            shear_energy.append(assemble(shear*dx(ii)))

            translational_single.append(assemble(translational*dx(ii)))
            rotational_single.append(assemble(rotational*dx(ii))/assemble(Constant(1.0)*dx(ii)))
            curvature_single.append(assemble(curvature*dx(ii)))
    
        return stretch_total, bend_total, shear_total, np.array(stretch_energy), np.array(bend_energy), np.array(shear_energy), \
                translational_total, rotational_total, curvature_total, np.array(translational_single), np.array(rotational_single), np.array(curvature_single)
    else:
        return stretch_total, bend_total, shear_total , translational_total, rotational_total, curvature_total
    
def compute_rxn_force(V, bc, residual):
    '''
    Calculates the reaction force at a given Dirichlet boundary condition
    Args:
        V: Functionspace of solution
        bc: boundary to evaluate the reaction forces at
        residual: total energy of the system
    Returns:
        reaction force
    '''
    v_reac = Function(V)
    bc.apply(v_reac.vector())

    return assemble(action(residual,v_reac))

def get_dof_coords(V):
    '''
    Gets the dof index and dof coordinates
    Args:
        V: function space to get dofs
    Returns:
        x_dofs: index of x dof
        y_dofs: index of y dof
        theta_dofs: index of rotation dof
        dof_coords: mesh coordinates of each dof
    '''
    x_dofs = V.sub(0).sub(0).dofmap().dofs()
    y_dofs = V.sub(0).sub(1).dofmap().dofs()
    theta_dofs = V.sub(1).dofmap().dofs()
    dofs = V.tabulate_dof_coordinates()
    dof_coords = dofs.reshape((-1, 2))

    return x_dofs, y_dofs, theta_dofs, dof_coords

def get_link_dof(V, link_points, tol=1e-6):
    '''
    Gets dof of crosslinks
    Args:
        V: function space to get dofs
        link_points: dof of crosslink
        tol: tolerance to get coordinates
    Returns:
        dofs of crosslinks 
    '''     
    x_dofs, y_dofs, theta_dofs, dof_coords = get_dof_coords(V)

    # Get x displacement at crosslinks
    x_disp_dof = []
    for ii in link_points:
        for jj in x_dofs:
            diff = np.sum((ii-dof_coords[jj])**2)
            if diff < tol:
                x_disp_dof.append(jj)
    y_disp_dof = []
    for ii in link_points:
        for jj in y_dofs:
            diff = np.sum((ii-dof_coords[jj])**2)
            if diff < tol:
                y_disp_dof.append(jj)
    
    return np.array([x_disp_dof, y_disp_dof])

def compute_damping_constant(c0, bcs, F_int, F_ext, v, dv):
    K = PETScMatrix()
    K = assemble(derivative(F_int-F_ext,v,dv), K)

    for bc in bcs:
        bc.apply(K)
    N = K.mat().getSize()[0]

    c = c0*(K.mat().getDiagonal().sum()/N)
    return c

def compute_damping_energy(damping_energy):
    return assemble(damping_energy)

def find_critical_disp(x0,x1,v0,v1,solver, disp, v, num_fibers, stretch,bend,shear,dx, tol):
    '''
    Finds the critical strain where the fiber switch from bending dominated to stretching dominated u
    using the bisection root finding method.
    Args:
        x0: displacemet right before transition
        x1 displacement right after transition
        solver: FEniCS Newton solver object
        disp: FEniCS expression object of applied nonzero displacement boundary condition
        v: solution vector
        num_fibers: number of fibers in simulation
        stretch: ufl expression to compute stretching energy
        bend: ufl expression to compute bending energy
        shear: ufl expression to compute shear energy
        dx: ufl integration domain
        tol: tolerance to stop convergence for bisection method

    Returns:
        x1: critical displacement of energy transition
    '''
    # Compute difference in energy at the intial guesses
    v.vector()[:]=v0
    stretch_total, bend_total, shear_total, _,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
    diff_energy0 = (stretch_total - bend_total)

    v.vector()[:]=v1
    stretch_total, bend_total, shear_total, _, _,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
    diff_energy1 = (stretch_total - bend_total)

    v.vector()[:] = v0 # better for initial guess for convergence

    x2 = (x0+x1)/2

    ii=0

    while(True):
        try:
            print(f'New Displacement: {x2}')
            ii+=1
            converged = False
            disp.t = x2
            iteration, converged = solver.solve()

        except RuntimeError:
            print('Newton solver failed during bisection process')
            x2 = np.nan
            break

        if converged:
            stretch_total, bend_total, shear_total,_,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            diff_energy2 = (stretch_total - bend_total)
            total_energy = stretch_total + bend_total + shear_total

            print(f'|Percent Energy difference:{diff_energy2/total_energy :.4e}| Displacement:{x2 :.4e}|')

            if diff_energy0*diff_energy2 <= 0:
                x1=x2
                diff_energy1 = diff_energy2
            elif diff_energy1*diff_energy2 < 0:
                x0=x2
                diff_energy0 = diff_energy2
            else:
                raise Exception('Boundary values do not contain root!')
            
            diff_x = x1-x0
            x2 = (x1+x0)/2
            v.vector()[:] = v0

            print(x0,x1)
            print(diff_energy0,diff_energy1)
            
            if np.abs(diff_energy2/total_energy) < tol or np.abs(diff_x/x2) < tol:
                break
     
        assert ii<=50, 'Root Finder Not Converging' # max 50 bisection iterations

    return x2, diff_energy2

def run_critical_strain(mesh_name, r, L, step_size, bc_type, c0 = None, rel_tol = 1e-8, abs_tol = 1e-9, relax = 1.0, quad = 3):

    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = quad
    parameters['reorder_dofs_serial'] = False

    mesh, mf, dx = read_mesh(mesh_name)

    init_step_size = step_size

    V,v,v_,dv,u,u_,theta,theta_ = beam_function_space(mesh)

    # define tangent vectors
    R0, g01 = compute_tangents(mesh)

    Rot = rotation_matrix(theta)

    # Boundary Conditions
    def top(x, on_boundary):
        return near(x[1],L,1e-6)

    def bot(x, on_boundary):
        return near(x[1],0,1e-6)
    

    # Case for boundary conditions:
    # bc_type[0] is the bottom BC and bc_type[1] is the top bc
    # 'f' is fixed, 'p' is pinned (prescibed displacement free rotation), 'h' is half-pinned (prescribed in y, free in x and rotation)

    # Bottom Dirichlet BC cases
    if bc_type[0] == 'f': 
        BC_bot = DirichletBC(V, Constant((0.0,0.0,0.0)), bot) # fixed displacement and rotation on the bottom
        bcs = [BC_bot]
    elif bc_type[0] == 'p':
        BC_bot = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot) # fixed displacement and rotation on the bottom
        bcs = [BC_bot]
    elif bc_type[0] == 'r': # for fiber networks only
        x_dofs = V.sub(0).sub(0).dofmap().dofs()
        dof_coords = V.tabulate_dof_coordinates().reshape((-1, 2))
        bot_node_coords = []
        for ii in x_dofs:
            if abs(dof_coords[ii,1]) <= 1e-6:
                    bot_node_coords.append(dof_coords[ii,0])
        
        bot_node_coords = np.array(bot_node_coords)

        def bot_x(x, on_boundary):
            return near(x[0],np.min(bot_node_coords)) and near(x[1],0,1e-6)

        BC_bot_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), bot_x) # fixed displacement and rotation on the bottom
        BC_bot_y = DirichletBC(V.sub(0).sub(1), Constant(0.0), bot) # fixed displacement and rotation on the bottom
        bcs = [BC_bot_x,BC_bot_y]
    else:
        raise Exception('Bottom BC type not supported!')

    # Top zero Dirichlet BC cases
    if bc_type[1] == 'f': 
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # fixed displacement and rotation on the bottom
        BC_top_rot = DirichletBC(V.sub(0).sub(1), Constant(0.0), top)
        bcs.append(BC_top)
        bcs.append(BC_top_rot)
    elif bc_type[1] == 'p':
        BC_bot = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # fixed displacement and rotation on the bottom
        bcs.append(BC_bot)
    elif bc_type[1] == 'r':
        pass
    else:
        raise Exception('Top BC type not supported!')

    # Applied displacement 
    disp = Expression("t",t=0,degree=0)
    BC_top_y = DirichletBC(V.sub(0).sub(1), disp, top)
    bcs.append(BC_top_y)

    # get strain measures:
    defo, curv = strains(R0,Rot,g01,u,theta)

    # Geometrical properties (circle cross section; 0 Poisson's ratio)
    r = Constant(r)
    S = pi*r**2
    I = pi*r**4/4
    kappa = Constant(6/7)

    # Stiffness moduli
    E = Constant(1e6)
    G = E/2

    # Beam stiffness parameters
    ES = E*S
    GS = kappa*G*S
    EI = E*I

    # Constitutive Equations
    C_N = diag(as_vector([ES, GS]))

    # Applied Load:
    F_max = Constant((0.0,0.0))
    M_max = Constant(0.0)

    elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

    F_int = derivative(elastic_energy, v, v_)
    F_ext = (M_max*theta_ + dot(F_max, u_)) * ds
    residual = F_int - F_ext
    # Add damping term if needed
    if c0 != None:
        damping_const = Expression("c", c=0.0, degree = 0) # Create expression to extend the top
        damping_energy = 0.5*damping_const*inner(v,v)*dx
        damping = derivative(damping_energy,v,v_)
        residual += damping
    tangent_form = derivative(residual, v, dv)

    #-------------------Define Newton Solver----------------------------------#
    problem = NonlinearVariationalProblem(residual, v, bcs, tangent_form)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters["newton_solver"]
    prm["linear_solver"] = "default"
    prm['relaxation_parameter'] = relax
    prm["absolute_tolerance"] = abs_tol*E # absolute tolerance should be with respect to E
    prm["relative_tolerance"] = rel_tol

    num_fibers = int(mf.array().max()+1) # number of fibers in the network

    # Energy functions to integrate
    stretch = 0.5 * defo[0]*dot(C_N, defo)[0]
    shear = 0.5 * defo[1]*dot(C_N, defo)[1]
    bend = 0.5*EI*curv**2

    reaction_force = []
    bcRy = DirichletBC(V.sub(0).sub(1), Constant(1.), top)

    dofs = V.tabulate_dof_coordinates()
    dof_coords = dofs.reshape((-1, 2))

    inc_cross = []
    v_cross = []

    # Iterate over load steps with Newton Solver
    while True:
        try:
            converged = False
            disp.t += step_size 

            if c0 != None:
                c = compute_damping_constant(c0, bcs, F_int, F_ext, v, dv)
                damping_const.c = c

            iteration,converged = solver.solve()
            v_converged = v.vector()[:]

            if c0 != None:
                damp = compute_damping_energy(damping_energy)
                percent_damp = damp/assemble(elastic_energy)
                
                if percent_damp > 0.02:
                    disp.t -= step_size
                    c0 = c0/10
                    converged = False
        except RuntimeError:
            print('Using smaller step:')
            disp.t -= step_size # go back to previously converged step
            step_size = step_size/2 # use smaller step 
            v.vector()[:] = v_converged # use previously converged solution as initial guess

            # Stop loop if step size is too small
            if init_step_size/step_size > 2**10:
                raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None
            
            continue
        if converged:
            #---------------------Compute energy------------------------------------------#
            stretch_total, bend_total, shear_total,_,_,_= postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            inc_cross.append(disp.t)
            v_cross.append(v.vector()[:])

            if stretch_total >= bend_total:
                print('Found Stretching dominated regime!')
                print('Starting Bisection Method')
                break

    tol = 1e-12 # convergence tolerance for secant method to find critical strain

    if len(inc_cross) == 1:
        while(True):
            v.vector().zero()
            step_size = step_size*1e-6
            disp.t = step_size

            solver.solve()
            stretch_total, bend_total, shear_total,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            total_energy = stretch_total + bend_total + shear_total
            print(stretch_total-bend_total)
            print(shear_total-stretch_total)

            if stretch_total <= bend_total:
                inc_cross.append(disp.t)
                v_cross.append(v.vector()[:])
                break
    
    if stretch_total == bend_total:
                diff_energy = 0
                crit_strain = 0
    else:
            crit_disp, diff_energy = find_critical_disp(inc_cross[-2],inc_cross[-1], v_cross[-2], v_cross[-1], solver, disp, v, num_fibers,stretch,bend,shear,dx, tol)
            crit_strain = crit_disp/L

    return crit_strain, diff_energy

def run_fea(mesh_name, f_name, r, L, step_size, bc_type, max_strain, c0 = None, compute_individual = False, rel_tol = 1e-8, abs_tol = 1e-9, relax = 1.0, quad = 3):

    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = quad
    parameters['reorder_dofs_serial'] = False

    mesh, mf, dx = read_mesh(mesh_name)

    init_step_size = step_size

    V,v,v_,dv,u,u_,theta,theta_ = beam_function_space(mesh)

    # define tangent vectors
    R0, g01 = compute_tangents(mesh)

    Rot = rotation_matrix(theta)

    # Boundary Conditions
    def top(x, on_boundary):
        return near(x[1],L,1e-6)

    def bot(x, on_boundary):
        return near(x[1],0,1e-6)

    # Case for boundary conditions:
    # bc_type[0] is the bottom BC and bc_type[1] is the top bc
    # 'f' is fixed, 'p' is pinned (prescibed displacement free rotation), 'h' is half-pinned (prescribed in y, free in x and rotation)

    # Bottom Dirichlet BC cases
    if bc_type[0] == 'f': 
        BC_bot = DirichletBC(V, Constant((0.0,0.0,0.0)), bot, method = 'pointwise') # fixed displacement and rotation on the bottom
        bcs = [BC_bot]
    elif bc_type[0] == 'p':
        BC_bot = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot, method = 'pointwise') # fixed displacement and rotation on the bottom
        bcs = [BC_bot]
    elif bc_type[0] == 'r': # for fiber networks only
        x_dofs = V.sub(0).sub(0).dofmap().dofs()

        dof_coords = V.tabulate_dof_coordinates().reshape((-1, 2))
        bot_node_coords = []
        for ii in x_dofs:
            if abs(dof_coords[ii,1]) <= 1e-6:
                    bot_node_coords.append(dof_coords[ii,0])
        
        bot_node_coords = np.array(bot_node_coords)

        def bot_x(x, on_boundary):
            return near(x[0],np.min(bot_node_coords)) and near(x[1],0,1e-6)

        BC_bot_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), bot_x, method = 'pointwise') # fixed displacement and rotation on the bottom
        BC_bot_y = DirichletBC(V.sub(0).sub(1), Constant(0.0), bot, method = 'pointwise') # fixed displacement and rotation on the bottom
        bcs = [BC_bot_x,BC_bot_y]
    else:
        raise Exception('Bottom BC type not supported!')

    # Top zero Dirichlet BC cases
    if bc_type[1] == 'f': 
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top, method = 'pointwise') # fixed displacement and rotation on the bottom
        BC_top_rot = DirichletBC(V.sub(0).sub(1), Constant(0.0), top, method = 'pointwise')
        bcs.append(BC_top)
        bcs.append(BC_top_rot)
    elif bc_type[1] == 'p':
        BC_bot = DirichletBC(V.sub(0).sub(0), Constant(0.0), top, method = 'pointwise') # fixed x displacement
        bcs.append(BC_bot)
    elif bc_type[1] == 'r':
        pass
    else:
        raise Exception('Top BC type not supported!')

    # Applied displacement 
    disp = Expression("t",t=0,degree=0)
    BC_top_y = DirichletBC(V.sub(0).sub(1), disp, top, method = 'pointwise')
    bcs.append(BC_top_y)

    # get strain measures:
    defo, curv = strains(R0,Rot,g01,u,theta)

    # Geometrical properties (circle cross section; 0 Poisson's ratio)
    r = Constant(r)
    S = pi*r**2
    I = pi*r**4/4
    kappa = Constant(6/7)

    # Stiffness moduli
    E = Constant(1e6)
    G = E/2

    # Beam stiffness parameters
    ES = E*S
    GS = kappa*G*S
    EI = E*I

    # Constitutive Equations
    C_N = diag(as_vector([ES, GS]))

    # Applied Load:
    F_max = Constant((0.0,0.0))
    M_max = Constant(0.0)

    elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

    F_int = derivative(elastic_energy, v, v_)
    F_ext = (M_max*theta_ + dot(F_max, u_)) * ds
    residual = F_int - F_ext

    # Add damping term if needed
    if c0 != None:
        damping_const = Expression("c", c=0.0, degree = 0) # Create expression to extend the top
        damping_energy = 0.5*damping_const*inner(v,v)*dx
        damping = derivative(damping_energy,v,v_)
        residual += damping

    tangent_form = derivative(residual, v, dv)

    #-------------------Define Newton Solver----------------------------------#
    problem = NonlinearVariationalProblem(residual, v, bcs, tangent_form)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters["newton_solver"]
    prm["linear_solver"] = "default"
    prm['relaxation_parameter'] = relax
    prm["absolute_tolerance"] = abs_tol*E # absolute tolerance should be with respect to E
    prm["relative_tolerance"] = rel_tol

    num_fibers = int(mf.array().max()+1) # number of fibers in the network

    # Energy functions to integrate
    stretch = 0.5 * defo[0]*dot(C_N, defo)[0]
    shear = 0.5 * defo[1]*dot(C_N, defo)[1]
    bend = 0.5*EI*curv**2

    reaction_force = []
    bcRy = DirichletBC(V.sub(0).sub(1), Constant(1.), top)

    dofs = V.tabulate_dof_coordinates()
    dof_coords = dofs.reshape((-1, 2))

    bend_all = []
    stretch_all = []
    shear_all = []
    bend_total_all = []
    stretch_total_all = []
    shear_total_all = []
    inc = []

    cont_all = []
    rot_all = []
    curv_all = []
    cont_total_all = []
    rot_total_all = []
    curv_total_all = []

    v_all = []

    damping_energy_all = []
    c0_all = []
    c_all = []
    stop = False
    # Iterate over load steps with Newton Solver
    while True:
        try:
            converged = False
            disp.t += step_size 
            if disp.t/L >= max_strain:
                disp.t = max_strain*L
                stop = True

            if c0 != None:
                c = compute_damping_constant(c0, bcs, F_int, F_ext, v, dv)
                damping_const.c = c

            iteration,converged = solver.solve()
            v_converged = v.vector()[:]

            if c0 != None:
                damp = compute_damping_energy(damping_energy)
                percent_damp = damp/assemble(elastic_energy)
                
                if percent_damp > 0.02: # make sure damping energy is below 2%
                    disp.t -= step_size
                    c0 = c0/10
                    converged = False
        except RuntimeError as err:
            print('Using smaller step:')
            disp.t -= step_size # go back to previously converged step
            step_size = step_size/2 # use smaller step 
            v.vector()[:] = v_converged # use previously converged solution as initial guess

            # Stop loop if step size is too small
            if init_step_size/step_size > 2**10:
                raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None
                
            continue
        if converged and stop == False:
            inc.append(disp.t)
            reaction_force.append(compute_rxn_force(V,bcRy, residual))
            v_all.append(v.vector()[:])
            if c0 != None:
                damping_energy_all.append(compute_damping_energy(damping_energy))
                c0_all.append(c0)
                c_all.append(c)
            #---------------------Compute energy------------------------------------------#
            if compute_individual:
                stretch_total, bend_total, shear_total, stretch_energy, bend_energy, shear_energy,\
                cont_total, rot_total, curv_total, cont, rot, curv_int = postprocess_fiber(num_fibers, stretch, bend, shear, dx,sqrt(dot(g01+tgrad(u,g01),g01+tgrad(u,g01))), theta, abs(curv), compute_individual)
        
                bend_all.append(bend_energy)
                stretch_all.append(stretch_energy)
                shear_all.append(shear_energy)

                bend_total_all.append(bend_total)
                stretch_total_all.append(stretch_total)
                shear_total_all.append(shear_total)

                cont_total_all.append(cont_total)
                rot_total_all.append(rot_total)
                curv_total_all.append(curv_total)
                cont_all.append(cont)
                rot_all.append(rot)
                curv_all.append(curv_int)       
            else:
                stretch_total, bend_total, shear_total, cont_total, rot_total, curv_total = postprocess_fiber(num_fibers, stretch, bend,  shear, dx, sqrt(dot(g01+tgrad(u,g01),g01+tgrad(u,g01))), theta, curv)

                bend_total_all.append(bend_total)
                stretch_total_all.append(stretch_total)
                shear_total_all.append(shear_total)

                cont_total_all.append(cont_total)
                rot_total_all.append(rot_total)
                curv_total_all.append(curv_total)
                
        #----------stop simulation once you hit the maximum displacement-------------#    
        elif converged and stop:

            inc.append(disp.t)
            reaction_force.append(compute_rxn_force(V,bcRy, residual))
            v_all.append(v.vector()[:])
            if c0 != None:
                damping_energy_all.append(compute_damping_energy(damping_energy))
                c0_all.append(c0)
                c_all.append(c)
            #---------------------Compute energy------------------------------------------#
            if compute_individual:
                stretch_total, bend_total, shear_total, stretch_energy, bend_energy, shear_energy,\
                cont_total, rot_total, curv_total, cont, rot, curv_int = postprocess_fiber(num_fibers, stretch, bend, shear, dx,sqrt(dot(g01+tgrad(u,g01),g01+tgrad(u,g01))), theta, abs(curv), compute_individual)
        
                bend_all.append(bend_energy)
                stretch_all.append(stretch_energy)
                shear_all.append(shear_energy)

                bend_total_all.append(bend_total)
                stretch_total_all.append(stretch_total)
                shear_total_all.append(shear_total)

                cont_total_all.append(cont_total)
                rot_total_all.append(rot_total)
                curv_total_all.append(curv_total)
                cont_all.append(cont)
                rot_all.append(rot)
                curv_all.append(curv_int)   
                print('Maximum applied strain achieved!')    
            else:
                stretch_total, bend_total, shear_total, cont_total, rot_total, curv_total = postprocess_fiber(num_fibers, stretch, bend,  shear, dx, sqrt(dot(g01+tgrad(u,g01),g01+tgrad(u,g01))), theta, curv)

                bend_total_all.append(bend_total)
                stretch_total_all.append(stretch_total)
                shear_total_all.append(shear_total)

                cont_total_all.append(cont_total)
                rot_total_all.append(rot_total)
                curv_total_all.append(curv_total)

                print('Maximum applied strain achieved!')
            break

    np.savetxt(f'{f_name}/disp.txt', np.array(inc))
    np.savetxt(f'{f_name}/force.txt', np.array(reaction_force))
    np.savetxt(f'{f_name}/stretch_total.txt', np.array(stretch_total_all))
    np.savetxt(f'{f_name}/bend_total.txt', np.array(bend_total_all))
    np.savetxt(f'{f_name}/shear_total.txt', np.array(shear_total_all))
    np.savetxt(f'{f_name}/cont_total.txt', np.array(cont_total_all))
    np.savetxt(f'{f_name}/rot_total.txt', np.array(rot_total_all))
    np.savetxt(f'{f_name}/curv_total.txt', np.array(curv_total_all))
    np.savetxt(f'{f_name}/v.txt', np.array(v_all))

    if c0 != None:
        np.savetxt(f'{f_name}/damping_energy.txt', np.array(damping_energy_all))
        np.savetxt(f'{f_name}/c0.txt', np.array(c0_all))
        np.savetxt(f'{f_name}/c.txt', np.array(c_all))

    if compute_individual:
        np.savetxt(f'{f_name}/stretch_energy.txt', np.array(stretch_all))
        np.savetxt(f'{f_name}/bend_energy.txt', np.array(bend_all))
        np.savetxt(f'{f_name}/shear_energy.txt', np.array(shear_all))
        np.savetxt(f'{f_name}/cont.txt', np.array(cont_all))
        np.savetxt(f'{f_name}/rot.txt', np.array(rot_all))
        np.savetxt(f'{f_name}/curv.txt', np.array(curv_all))
    else:
        pass