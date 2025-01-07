from dolfin import *
import numpy as np
from ufl import Jacobian, diag, sign
import matplotlib.pyplot as plt
import copy
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
    Create mixed function space for beams base on the Hellinger-Reissner Principle
    Args: 
        mesh: dolfin mesh file
    Returns:
        V: mixed function space
        dv: trial function space
        v: solution function
        u: displacement solution function
        u_: displacement test function space
        du: displacement trial function space
        theta : rotation solution function
        theta_: rotation test function space
        dtheta: rotation tr5ial function space
        n : force solution function
        n_: force test function space
        dn: force tr5ial function space
        m : moment solution function
        m_: moment test function space
        dm: moment tr5ial function space
    
    '''
    

    # Define function space
    Ue = VectorElement("CG", mesh.ufl_cell(), 2, dim=2) # displacement field
    Te = FiniteElement("CG", mesh.ufl_cell(), 1) # rotation field

    Ne = VectorElement("DG", mesh.ufl_cell(), 1, dim=2) # force field
    Me = FiniteElement("DG", mesh.ufl_cell(),0) # moment field


    V = FunctionSpace(mesh, MixedElement([Ue, Te, Ne, Me]))   

    v_ = TestFunction(V)
    u_, theta_, n_, m_ = split(v_)
    dv = TrialFunction(V)
    du,dtheta,dn,dm = split(dv)
    v = Function(V, name="Solution Function")
    u, theta, n, m = split(v)

    return V, v, v_, dv, u, u_,du,dtheta,dn,dm, theta, theta_, n, n_, m, m_


def compute_tangents(mesh):
    '''
    Define tangent vectors and collect them into an initial rotation matrix
    Args:
        mesh: dolfin mesh file

    Returns:
        R0: Initial rotation matrix
        g01: initial tangent direction
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
        current rotation tensir

    '''
    return as_tensor([[cos(theta), -sin(theta)],[sin(theta), cos(theta)]])

def tgrad(u, g01):
    return dot(grad(u), g01)

def strains(R0,Rot,g01,u,theta):
    '''
    strain measures for beams (axial, shear, curvature)
    Args:
        R0: initial rotation tensor (element to global rotation)
        Rot: current rotation tensor
        g01: initial tangent direction
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
        stretch : ufl stretching energy integrand
        bend: ufl bending energy integrand
        shear: ufl shearning energy integrand
        dx : integration domain (make sure its partioned by fibers)
        translational: ufl displacement integrand
        rotational: ufl rotation integrand
        curvature: ufl curvature integrand
        compute_individual: True if needed to compute energy individual fibers (will take longer computationally); False by default 
    Returns:
        stretch_total: total stretch energy
        bend_total: total bend energy
        shear_total: total shear energy
        translational_total: total length change(doesn't mean anything physically; used to check single calculations)
        rotational_total: total rotational (doesn't mean anything phyiscally; used to check single calculations)
        curvature_total: total curvature (doesn't mean anything phyiscally; used to check single calculation)
        stretch_energy: stretch energy of each fiber (only returns if compute_individual = True)
        bend_energy: bend energy of each fiber
        shear_energy: shear energy of each fiber
        translational_single: total length change of one fiber 
        rotational_single: total rotatation of one fiber normalized by fiber length
        curvature_single: total rotation of one fiber

         
    '''
    # conpute total values of all fibers
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

    # compute values per fiber
    if compute_individual:
        for ii in range(num_fibers):
            stretch_energy.append(assemble(stretch*dx(ii)))
            bend_energy.append(assemble(bend*dx(ii)))
            shear_energy.append(assemble(shear*dx(ii)))

            translational_single.append(assemble(translational*dx(ii)))
            rotational_single.append(assemble(rotational*dx(ii))/assemble(Constant(1.0)*dx(ii))) # make sure normalized by fiber length
            curvature_single.append(assemble(curvature*dx(ii)))
    
        return stretch_total, bend_total, shear_total, np.array(stretch_energy), np.array(bend_energy), np.array(shear_energy), \
                translational_total, rotational_total, curvature_total, np.array(translational_single), np.array(rotational_single), np.array(curvature_single)
    else:
        return stretch_total, bend_total, shear_total , translational_total, rotational_total, curvature_total
    
def compute_rxn_force(V, bc, residual):
    '''
    Calculates the reaction force at a given Dirichlet boundary condition
    Args:
        V: Function space of solution
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

def compute_damping_constant(c0, bcs, stretching, bending, u,theta,du, dtheta):
    '''
    Compute damping constant for Tokhonov regularization
    Args:
        c0: initial c0 (prescribed)
        bcs: boundary conditions 
        stretching: ufl form translational strain energy
        bending: ufl form bending strain energy
        u: displacement solution function
        theta: rotation solution function
        du: displacement trial function space
        dtheta: rotation trial solution
    Returns:
        c_N: damping constant for translation
        c_M: damping constant for rotation
    '''

    # Stretching component
    K_N = PETScMatrix()
    K_N = assemble(derivative(stretching,u,du), K_N)

    for bc in bcs:
        bc.apply(K_N)
    N = K_N.mat().getSize()[0]
    c_N = c0*(K_N.mat().getDiagonal().sum()/N)

    # Bending component
    K_M = PETScMatrix()
    K_M = assemble(derivative(bending,theta,dtheta), K_M)

    for bc in bcs:
        bc.apply(K_M)
    N = K_M.mat().getSize()[0]

    c_M = c0*(K_M.mat().getDiagonal().sum()/N)
    return c_N, c_M

def compute_damping_energy(damping_energy):
    '''
    Compute "energy" of Tikhonov regularization (addtional term in strain energy)
    Args:
        damping_energy: ufl form of damping energy
    Returns:
        damping energy value
    '''
    return assemble(damping_energy)

def compute_residual_norm(v, dv,  bcs, residual):
    res = assemble(residual)
    for bc in bcs:
        bc.apply(res)
    return res.norm('l2')

def find_critical_disp(x0,x1,v0,v1,solver, disp, V, v, v_, dv, bcRy, num_fibers, stretch_damp,bend_damp,stretch,bend,shear, residual, damping_energy, damping_const_N,damping_const_M, c0, max_damping_res, init_c0, bcs, u_p, theta_p,dx, u, du, theta,dtheta, tol):
    '''
    Finds the critical strain where the fiber switch from bending dominated to stretching dominated using the bisection root finding method.

    Args:
        x0: displacemet right before transition
        x1 displacement right after transition
        solver: FEniCS Newton solver object
        disp: FEniCS expression object of applied nonzero displacement boundary condition
        v: solution vector
        num_fibers: number of fibers in simulation
        stretch_damp: ufl exoression to compute stretching part of regularization
        bend_damp: ufl expression to compute bending part of regularization
        stretch: ufl expression to compute stretching energy
        bend: ufl expression to compute bending energy
        shear: ufl expression to compute shear energy
        residual: residual of fea problem
        damping_energy: total damping energy
        damping_const_N: damping constant of stretching part
        damping_cont_M: damping constant of bending part
        c0: damping constant
        max_damping_res: max residual 
        init_c0: initial damping constant
        bcs: Dirichlet bcs in the FEA problem
        u_p: previously converged displacement solution
        theta_p: previously converged rotational solution
        dx: ufl integration domain
        u: displacement solution
        du: trial function space
        theta: rotation function space
        dtheta: trial function space
        tol: tolerance to stop convergence for bisection method

    Returns:
        x2: critical displacement of energy transition
        diff_energy2: difference betrween bending and stretching energy
    '''
    # Compute difference in energy at the intial guesses
    v.vector()[:]=v0
    stretch_total, bend_total, shear_total, _,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
    diff_energy0 = (stretch_total - bend_total)

    v.vector()[:]=v1
    stretch_total, bend_total, shear_total, _, _,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
    diff_energy1 = (stretch_total - bend_total)

    x2 = (x0+x1)/2

    print(f'Iniial Bounds:{[x0,x1]}')
    ii=0

    v_converged = v0

    v.vector()[:]=v0

    u_p_init = v.sub(0,True).vector()[:] # for damping energy
    theta_p_init = v.sub(1,True).vector()[:] # for damping energy

    u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
    theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy

    flag = False
    if init_c0 != None:
        damping = derivative(damping_energy*dx,v,v_)
    else:
        damping = 0.0
    ii=0
    jj=0

    while(True):
        try:
            print(f'New Displacement: {x2}')
            converged = False
            disp.t = x2

            if c0 != None:
                c_N,c_M = compute_damping_constant(c0, bcs, stretch_damp, bend_damp, u,theta, du,dtheta)
                damping_const_N.c = c_N
                damping_const_M.c = c_M
            
            iteration, converged = solver.solve()

            current_disp = disp.t
            disp.t = 0
            computed_residual = compute_residual_norm(v,dv, bcs, residual-damping)
            f_reac = compute_rxn_force(V,bcRy, residual-damping)
            disp.t = current_disp

            # relaxing the damping parameter
            if init_c0 != None:
                damp = compute_damping_energy(damping_energy*dx)
                percent_damp = damp/assemble((stretch+bend+shear)*dx)
                if computed_residual/f_reac > max_damping_res:
                    print(f'|c0: {c0:.2e}|Residual Iteration: {ii}|Equilibrium Residual Ratio: {computed_residual/f_reac:.2e} (damped tol: {max_damping_res:.2e})|')
                    if ii > 30:
                        c0 = c0/2
                        ii=0
                    converged = False
                    u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
                    theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy
                    ii+=1

        except RuntimeError:
            if init_c0 != None:
                if c0 != 0.0:
                    c0 = 2*c0
                else:
                    c0 = init_c0
                v.vector()[:]=v_converged
                u_p.vector()[:] = u_p_init
                theta_p.vector()[:]=theta_p_init
                ii=0


        if converged:
            ii=0
            jj+=1
            stretch_total, bend_total, shear_total,_,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            diff_energy2 = (stretch_total - bend_total)
            total_energy = stretch_total + bend_total + shear_total

            diff_x = x1-x0
            
            print(f'Bisection tol: {tol :.4e} |Percent Energy difference:{diff_energy2/total_energy :.4e}| Percent Displacement Difference {np.abs(diff_x/x2)}|')

            if diff_energy0*diff_energy2 <= 0: # upper bound solution
                x1=x2
                diff_energy1 = diff_energy2
                v.vector()[:] = v_converged
                u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
                theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy
            elif diff_energy1*diff_energy2 < 0: # lower bound solution
                x0=x2
                diff_energy0 = diff_energy2
                v_converged = v.vector()[:]
                u_p.vector()[:] = u_p_init # for damping energy
                theta_p.vector()[:] = theta_p_init # for damping energy
    
            else:
                raise Exception('Boundary values do not contain root!')
            
            x2 = (x1+x0)/2

            print(x0,x1)
            print(diff_energy0,diff_energy1)
            
            if np.abs(diff_energy2/total_energy) < tol or np.abs(diff_x/x2) < tol:
                break
     
        assert jj<=50, 'Root Finder Not Converging' # max 50 bisection iterations

    return x2, diff_energy2

def run_critical_strain(mesh_name, r, L, step_size, bc_type, init_c0 = None, max_damping_res = 1e-4, rel_tol = 1e-8, abs_tol = 1e-9, relax = 1.0, quad = 3):
    '''
    Runs pipeline to compute the critical strain trnaisiton point

    Args:
        mesh_name: name of mesh
        r: radius of fibers
        L: initial domain size
        step_size: Incremental step for Newton iteration
        bc_type: type of bc choices enters as list of characters, choices as 'f' for clamped bcs, 'p' for pinned bcs, 'r' for roller supports (i.e. ['p','p'] for pinned-pinned supports)
        init_c0: initial damping constant 
        max_damping_res: equilibrium tolerance allowed
        rel_tol : relative tolerance for Newton solver
        abs_tol: absolute tolerance for Newton solver
        relax: Newton relaxation parameter
        quad: degree of Gaussian quadrature
    Returns:
        crit_strain: critical strain transition point
        diff_energy: difference between stretching and bending energy
    '''
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = quad
    parameters['reorder_dofs_serial'] = False

    mesh, mf, dx = read_mesh(mesh_name)

    init_step_size = step_size

    V, v, v_, dv, u, u_,du,dtheta,dn,dm, theta, theta_, n, n_, m, m_ = beam_function_space(mesh)

    # Set previous solutions for regularization
    Vu = V.sub(0).collapse()
    u_p = Function(Vu, name = 'previous displacements')

    Vt = V.sub(1).collapse()
    theta_p = Function(Vt, name = 'previous rotations')
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
        BC_bot_disp = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot) # fixed displacement
        BC_bot_rot = DirichletBC(V.sub(1), Constant(0.0), bot) # fixed rotation
        bcs = [BC_bot_disp, BC_bot_rot]
    elif bc_type[0] == 'p':
        BC_bot = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot) # fixed displacement free rotation
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
        BC_bot_y = DirichletBC(V.sub(0).sub(1), Constant(0.0), bot) # fixed displacement and rotation on the bottom
        bcs = [BC_bot_x,BC_bot_y]
    else:
        raise Exception('Bottom BC type not supported!')

    # Top zero Dirichlet BC cases
    if bc_type[1] == 'f': 
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # prescribed x displacement
        BC_top_rot = DirichletBC(V.sub(0).sub(1), Constant(0.0), top) # fixed rotation
        bcs.append(BC_top)
        bcs.append(BC_top_rot)
    elif bc_type[1] == 'p':
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # fixed x displacement
        bcs.append(BC_top)
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
    E = Constant(1.0)
    G = E/2

    # Beam stiffness parameters
    ES = E*S
    GS = kappa*G*S
    EI = E*I

    # Constitutive Equations
    C_N = diag(as_vector([ES, GS]))
    S_N = diag(as_vector([1/ES, 1/GS]))

    # Applied Load (None in this case):
    F_max = Constant((0.0,0.0))
    M_max = Constant(0.0)

    # Rotate force field to the same frame
    N = dot(R0.T,dot(Rot.T,n))

    HR_axial_int = dot(N,defo) - 0.5*(dot(N,dot(S_N,N))) # axial component of HR functional ; equivalent to strain energy
    HR_rot_int = m*curv - 0.5*(1/EI)*m**2 # rotational component of HR functional
    HR_int = (HR_axial_int + HR_rot_int)*dx

    elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

    F_int = derivative(HR_int, v, v_)
    F_ext = (M_max*theta_ + dot(F_max, u_)) * ds
    residual = F_int - F_ext
    # Add damping term if needed
    if init_c0 != None:
        c0 = 0.0 # No damping unless needed
        damping_const_N = Expression("c", c=0.0, degree = 0)
        damping_const_M = Expression("c", c=0.0, degree = 0)
        damping_energy = 0.5*(damping_const_N*(inner(u-u_p,u-u_p))+damping_const_M*inner(theta-theta_p,theta-theta_p))
        damping = derivative(damping_energy*dx,v,v_)
        residual += damping
    else:
        c0 = None
        damping_const = None
        damping_energy = None
        damping_const_N = None
        damping_const_M = None
    tangent_form = derivative(residual, v, dv)

    #-------------------Define Newton Solver----------------------------------#
    problem = NonlinearVariationalProblem(residual, v, bcs, tangent_form)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters["newton_solver"]
    prm["linear_solver"] = "default"
    prm['relaxation_parameter'] = relax
    prm["absolute_tolerance"] = abs_tol
    prm["relative_tolerance"] = rel_tol
    prm["maximum_iterations"] = 20

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

    prev_converged = True
    prev_damp = None
    # Iterate over load steps with Newton Solver
    v_converged = v.vector()[:] # Make sure first guess is from a converged solution
    ii=0
    prev_computed_residual = 1.0
    while True:
        try:
            converged = False

            disp.t += step_size 

            print(f'|Step size: {step_size}| Displacement: {disp.t}|')

            if c0 != None:
                c_N,c_M = compute_damping_constant(c0, bcs, derivative(0.5*dot(defo, dot(C_N, defo))*dx,u,u_), derivative(0.5*EI*curv**2,theta,theta_)*dx,u,theta, du,dtheta)
                damping_const_N.c = c_N
                damping_const_M.c = c_M

            iteration,converged = solver.solve()
            prev_converged = converged

            if init_c0 != None and c0 != 0.0:
                damp = compute_damping_energy(damping_energy*dx)
                percent_damp = damp/assemble(elastic_energy)
        
                current_disp = disp.t
                disp.t = 0
                computed_residual = compute_residual_norm(v,dv, bcs, residual-damping)
                f_reac = compute_rxn_force(V,bcRy, residual-damping)
                disp.t = current_disp

                if computed_residual/f_reac > max_damping_res: # case where equilibrium criteria is not met
                    print(f'|c0: {c0:.2e}|Residual Iteration: {ii}|Equilibrium Residual Ratio: {computed_residual/f_reac:.2e} (damped tol: {max_damping_res:.2e})|')
                    if ii > 30:
                        c0 = c0/2
                        ii = 0
                    disp.t -= step_size
                    converged = False
                    u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
                    theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy
                    ii+=1



        except RuntimeError:

            computed_residual = compute_residual_norm(v,dv, bcs, residual-damping)
            f_reac = compute_rxn_force(V,bcRy, residual-damping)

            if c0 != None and np.abs(computed_residual/f_reac) < max_damping_res:
                continue


            disp.t -= step_size # go back to previously converged step

            # for case with damping
            if init_c0 != None:
                if init_step_size/step_size < 2**10:
                    if c0 == 0.0:
                        step_size=step_size/2
                    elif c0 != 0:
                        step_size=step_size/2
                        c0 = 2*c0
                else:
                    if c0 == 0.0: # start damping if needed
                        c0=init_c0
                        step_size=init_step_size
                    else:
                        raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None

            elif init_c0 == None:
                if init_step_size/step_size < 2**10:
                    step_size=step_size/2
                else:
                    raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None

            v.vector()[:] = v_converged # use previously converged solution as initial guess

            u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
            theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy

            prev_converged = False
            ii=0
            
            continue
        if converged:
            ii=0
            v_converged = v.vector()[:] # For Newton convergnece -- not related to damping
            #---------------------Compute energy------------------------------------------#
            stretch_total, bend_total, shear_total,_,_,_= postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            inc_cross.append(disp.t)
            v_cross.append(v.vector()[:])

            u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
            theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy

            if stretch_total >= bend_total:
                print('Found Stretching dominated regime!')
                print('Starting Bisection Method')
                break
            if prev_converged:
                step_size = min(step_size*2, init_step_size)
            
            prev_converged = True


    tol = 1e-12 # convergence tolerance for secant method to find critical strain

    if len(inc_cross) == 1:
        print('Finding Lower Bound:')
        while(True):
            v.vector().zero()
            step_size = step_size/2
            disp.t = step_size

            solver.solve()
            stretch_total, bend_total, shear_total,_,_,_ = postprocess_fiber(num_fibers, stretch, bend, shear, dx)
            total_energy = stretch_total + bend_total + shear_total
            print(stretch_total-bend_total)

            if stretch_total <= bend_total:
                inc_cross = [disp.t] + inc_cross
                v_cross = [v.vector()[:]] + v_cross
                break
    
    crit_disp, diff_energy = find_critical_disp(inc_cross[-2],inc_cross[-1], v_cross[-2], v_cross[-1], solver, disp, V, v, v_, dv, bcRy, num_fibers,derivative(0.5*dot(defo, dot(C_N, defo))*dx,u,u_), derivative(0.5*EI*curv**2,theta,theta_)*dx, stretch,bend,shear,residual, damping_energy, damping_const_N,damping_const_M, c0, max_damping_res, init_c0, bcs, u_p, theta_p,dx,u, du,theta,dtheta, tol)
    crit_strain = crit_disp/L

    return crit_strain, diff_energy

def run_fea(mesh_name, f_name, r, L, step_size, bc_type, max_strain, init_c0 = None, max_damping_res = 5e-4, compute_individual = False, compute_residual = False, rel_tol = 1e-8, abs_tol = 1e-9, relax = 1.0, quad = 3):
    '''
    Runs pipeline to compute the critical strain trnaisiton point

    Args:
        mesh_name: name of mesh
        f_name: folder to save results
        r: radius of fibers
        L: initial domain size
        step_size: Incremental step for Newton iteration
        bc_type: type of bc choices enters as list of characters, choices as 'f' for clamped bcs, 'p' for pinned bcs, 'r' for roller supports (i.e. ['p','p'] for pinned-pinned supports)
        init_c0: initial damping constant 
        max_damping_res: equilibrium tolerance allowed
        rel_tol : relative tolerance for Newton solver
        abs_tol: absolute tolerance for Newton solver
        relax: Newton relaxation parameter
        quad: degree of Gaussian quadrature
    Returns:
        None; all the results are saved in the f_name folder
        results obtained from this function are:
        disp.txt: applied displacement
        force.txt: applied force
        stretch_total.txt: total stretching energy
        bend_total.txt: total bending energy
        shear_total.txt: total shear energy
        cont_total.txt: total change in contour length of all fibers (no real physical meaning) 
        rot_total.txt: total rotation in all fibers(no physical meaning)
        curv_total.txt: total curvature on all fibers(no physical meaning)
        v.txt: displacement solution vector
    Optional Returns:
        damping_energy.txt: damping energy from regularization
        stretch_energy.txt: stretching energy per fiber
        bend_energy.txt: bending energy per fiber
        shear_energy.txt: shear energy per fiber
        cont.txt: contour length per fiber
        rot.txt: rotation per fiber
        curv.txt: curvature per fiber

    '''
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = quad
    parameters['reorder_dofs_serial'] = False

    mesh, mf, dx = read_mesh(mesh_name)

    init_step_size = step_size

    V, v, v_, dv, u, u_,du,dtheta,dn,dm, theta, theta_, n, n_, m, m_ = beam_function_space(mesh)

    # Set previous solutions for regularization
    Vu = V.sub(0).collapse()
    u_p = Function(Vu, name = 'previous displacements')

    Vt = V.sub(1).collapse()
    theta_p = Function(Vt, name = 'previous rotations')
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
        BC_bot_disp = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot) # fixed displacement
        BC_bot_rot = DirichletBC(V.sub(1), Constant(0.0), bot) # fixed rotation
        bcs = [BC_bot_disp, BC_bot_rot]
    elif bc_type[0] == 'p':
        BC_bot = DirichletBC(V.sub(0), Constant((0.0,0.0)), bot) # fixed displacement free rotation
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
        BC_bot_y = DirichletBC(V.sub(0).sub(1), Constant(0.0), bot) # fixed displacement and rotation on the bottom
        bcs = [BC_bot_x,BC_bot_y]
    else:
        raise Exception('Bottom BC type not supported!')

    # Top zero Dirichlet BC cases
    if bc_type[1] == 'f': 
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # prescribed x displacement
        BC_top_rot = DirichletBC(V.sub(0).sub(1), Constant(0.0), top) # fixed rotation
        bcs.append(BC_top)
        bcs.append(BC_top_rot)
    elif bc_type[1] == 'p':
        BC_top = DirichletBC(V.sub(0).sub(0), Constant(0.0), top) # fixed x displacement
        bcs.append(BC_top)
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
    E = Constant(1.0)
    G = E/2

    # Beam stiffness parameters
    ES = E*S
    GS = kappa*G*S
    EI = E*I

    # Constitutive Equations
    C_N = diag(as_vector([ES, GS]))
    S_N = diag(as_vector([1/ES, 1/GS]))

    # Applied Load (None in this case):
    F_max = Constant((0.0,0.0))
    M_max = Constant(0.0)

    # Rotate force field to the same frame
    N = dot(R0.T,dot(Rot.T,n))

    HR_axial_int = dot(N,defo) - 0.5*(dot(N,dot(S_N,N))) # axial component of HR functional ; equivalent to strain energy
    HR_rot_int = m*curv - 0.5*(1/EI)*m**2 # rotational component of HR functional
    HR_int = (HR_axial_int + HR_rot_int)*dx

    elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

    F_int = derivative(HR_int, v, v_)
    F_ext = (M_max*theta_ + dot(F_max, u_)) * ds
    residual = F_int - F_ext
    # Add damping term if needed
    if init_c0 != None:
        c0 = 0.0 # No damping unless needed
        damping_const_N = Expression("c", c=0.0, degree = 0)
        damping_const_M = Expression("c", c=0.0, degree = 0)
        damping_energy = 0.5*(damping_const_N*(inner(u-u_p,u-u_p))+damping_const_M*inner(theta-theta_p,theta-theta_p))
        damping = derivative(damping_energy*dx,v,v_)
        residual += damping
    else:
        c0 = None
        damping = 0.0
    tangent_form = derivative(residual, v, dv)

    #-------------------Define Newton Solver----------------------------------#
    problem = NonlinearVariationalProblem(residual, v, bcs, tangent_form)
    solver = NonlinearVariationalSolver(problem)
    prm = solver.parameters["newton_solver"]
    prm["linear_solver"] = "default"
    prm['relaxation_parameter'] = relax
    prm["absolute_tolerance"] = abs_tol # absolute tolerance should be with respect to E
    prm["relative_tolerance"] = rel_tol
    prm["maximum_iterations"] = 20

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
    residual_all = []

    prev_converged = True
    prev_damp = None
    stop = False
    # Iterate over load steps with Newton Solver
    v_converged = v.vector()[:] # make sure initial guess is from a converged solution
    while True:
        try:
            converged = False

            disp.t += step_size
            if disp.t/L >= max_strain:
                stop = True
                disp.t = max_strain*L

            if c0 != None:
                print(f'|Step size: {step_size}| Displacement: {disp.t}| Max displacement: {max_strain*L}|c0: {c0:.2e}|')
                c_N,c_M = compute_damping_constant(c0, bcs, derivative(0.5*dot(defo, dot(C_N, defo))*dx,u,u_), derivative(0.5*EI*curv**2,theta,theta_)*dx,u,theta, du,dtheta)
                damping_const_N.c = c_N
                damping_const_M.c = c_M

            iteration,converged = solver.solve()
            v_converged = v.vector()[:]
            prev_converged = converged

            if init_c0 != None and c0 != 0.0:
                damp = compute_damping_energy(damping_energy*dx)
                percent_damp = damp/assemble(elastic_energy)
        
                current_disp = disp.t
                disp.t = 0
                computed_residual = compute_residual_norm(v,dv, bcs, residual-damping)
                f_reac = compute_rxn_force(V,bcRy, residual-damping)
                disp.t = current_disp

                if computed_residual/f_reac > max_damping_res: # case where equilibrium criteria is not met
                    print(f'|c0: {c0:.2e}|Residual Iteration: {ii}|Equilibrium Residual Ratio: {computed_residual/f_reac:.2e} (damped tol: {max_damping_res:.2e})|')
                    if ii > 30:
                        c0 = c0/2
                        ii = 0
                    disp.t -= step_size
                    converged = False
                    u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
                    theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy
                    ii+=1

        except RuntimeError:
            computed_residual = compute_residual_norm(v,dv, bcs, residual-damping)
            f_reac = compute_rxn_force(V,bcRy, residual-damping)

            if c0 != None and np.abs(computed_residual/f_reac) < max_damping_res:
                continue


            disp.t -= step_size # go back to previously converged step

            # case with damping
            if init_c0 != None:
                if init_step_size/step_size < 2**10:
                    if c0 == 0.0:
                        step_size=step_size/2
                    elif c0 != 0:
                        step_size=step_size/2
                        c0 = 2*c0
                else:
                    if c0 == 0.0: # start damping if none initially
                        c0=init_c0
                        step_size=init_step_size
                    else:
                        raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None

            elif init_c0 == None:
                if init_step_size/step_size < 2**10:
                    step_size=step_size/2
                else:
                    raise RuntimeError(f'Newton solver not converging with step size {step_size}') from None

            v.vector()[:] = v_converged # use previously converged solution as initial guess

            u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
            theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy

            prev_converged = False
            ii=0
            
            continue
        if converged and stop == False:
            ii=0
            v_converged = v.vector()[:] # For Newton convergence -- not related to damping

            u_p.vector()[:] = v.sub(0,True).vector()[:] # for damping energy
            theta_p.vector()[:] = v.sub(1,True).vector()[:] # for damping energy

            inc.append(disp.t)
            reaction_force.append(compute_rxn_force(V,bcRy, residual-damping))
            v_all.append(v.vector()[:])
            
            if prev_converged:
                step_size = min(step_size*2, init_step_size)
            
            if compute_residual:
                current_disp = disp.t
                disp.t = 0 # homogenize BCs to conform to test function
                residual_val = compute_residual_norm(v,dv, bcs, residual-damping)
                residual_all.append(residual_val)
                disp.t = current_disp
            
            prev_converged = True

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
            reaction_force.append(compute_rxn_force(V,bcRy, residual-damping))
            v_all.append(v.sub(0).vector()[:]) # only store displacement
            if c0 != None:
                damping_energy_all.append(compute_damping_energy(damping_energy*dx))
            
            if compute_residual:
                current_disp = disp.t
                disp.t = 0 # homogenize BCs for correct test function
                residual_val = compute_residual_norm(v,dv, bcs, residual-damping)
                residual_all.append(residual_val)
                disp.t = current_disp

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
    #----------------Save Solution-------------------------------------#
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

    if compute_individual:
        np.savetxt(f'{f_name}/stretch_energy.txt', np.array(stretch_all))
        np.savetxt(f'{f_name}/bend_energy.txt', np.array(bend_all))
        np.savetxt(f'{f_name}/shear_energy.txt', np.array(shear_all))
        np.savetxt(f'{f_name}/cont.txt', np.array(cont_all))
        np.savetxt(f'{f_name}/rot.txt', np.array(rot_all))
        np.savetxt(f'{f_name}/curv.txt', np.array(curv_all))
    if residual:
        np.savetxt(f'{f_name}/residual.txt', np.array(residual_all))