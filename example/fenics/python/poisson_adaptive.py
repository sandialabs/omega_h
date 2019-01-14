# Solving the Poisson problem adaptively using Omega_h

from dolfin import *;
from mshr import *;
# In case omega_h is not in the default path, append the path where the PyOmega_h shared library is located
#import sys
#sys.path.append("/usr/local/lib/python3/dist-packages/")
import PyOmega_h as omega_h;

#Init dolfin mesh

# Build mesh with omega_h & initialize metrics
comm_osh = omega_h.world()
mesh_osh = omega_h.build_box(comm_osh, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 32, 32, 0)
mesh_osh.balance()
mesh_osh.set_parting(omega_h.GHOSTED, 0)
omega_h.add_implied_metric_tag(mesh_osh)
mesh_osh.set_parting(omega_h.ELEM_BASED, 0)

maxiter = 20
i = 0

def boundary(x):
	return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

while(i < maxiter):
	# Import mesh from omega_h to dolfin 
    mesh=Mesh()
    omega_h.mesh_to_dolfin(mesh, mesh_osh)

    # Create function space
    V = FunctionSpace(mesh, "Lagrange", 1)
	   
    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
    g = Expression("sin(5*x[0])", degree=2)
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx + g*v*ds

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    # Save solution in VTK format
    file = File("poisson_u_"+ str(i) +".pvd")
    file << u

    # Import u from dolfin to omega_h
    omega_h.function_from_dolfin(mesh_osh, u._cpp_object, "u")
 
    # Set up metric, adaptivity parameters
    mesh_osh.set_parting(omega_h.GHOSTED, 0);
    metric_input = omega_h.MetricInput()
    source = omega_h.MetricSource(omega_h.VARIATION, 2e-3, "u")
    metric_input.add_source(source)
    metric_input.should_limit_lengths = True
    metric_input.max_length = 1.0 / 2.0
    metric_input.should_limit_gradation = True
    omega_h.generate_target_metric_tag(mesh_osh, metric_input) 
    opts = omega_h.AdaptOpts(mesh_osh)
    opts.verbosity = omega_h.EXTRA_STATS
		
    # Adapt mesh
    while(omega_h.approach_metric(mesh_osh, opts)):
        omega_h.adapt(mesh_osh, opts)

    i+=1

