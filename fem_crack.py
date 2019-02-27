from fenics import *

import os, tempfile, shutil

# gmshloc = "/Applications/Gmsh.app/Contents/MacOS/gmsh"
gmshloc = "gmsh"
dbloc = "femdata.db"
geo = os.path.abspath("crack_gen.geo")

# Some global properties
E = 1.0
nu= 0.25

def sim(x1,y1, x2,y2, clscale):
    global dx, ds
    # Make a temporary directory to make a new set of meshes
    tmpdir = tempfile.mkdtemp()
    os.system("cp "+geo+" "+tmpdir+"/crack.geo")
    arg_pass_string = " ".join(["-setnumber {0} {1}".format(a,b)
                                    for a,b in
                                    [("x1",x1),("y1",y1),("x2",x2),("y2",y2),] ])
    os.system(gmshloc+" "+tmpdir+"/"+"crack.geo -clscale {0} {1} - ".format(clscale,arg_pass_string))
    os.system("dolfin-convert {0}/crack2.msh {0}/crack2.xml".format(tmpdir))
    # Open up the mesh files and define integrals
    mesh = Mesh(tmpdir+"/crack2.xml")
    cellids = MeshFunction("size_t",mesh,tmpdir+"/crack2_physical_region.xml")
    facetids = MeshFunction( "size_t", mesh,tmpdir+"/crack2_facet_region.xml")
    dx = dx(subdomain_data=cellids)
    ds = ds(subdomain_data=facetids)
    # Define the function space and solution/test/trial functions
    V = VectorFunctionSpace(mesh,"CG",1)
    u,tu,Du = Function(V,name='u'),TestFunction(V),TrialFunction(V)
    # Define the weak form
    n = FacetNormal(mesh)
    K = E/(3.0-6.0*nu)
    G = E/(2.0+2.0*nu)
    P = 1
    Dsigma = (K-2.0/3.0*G)*tr(sym(grad(Du)))*Identity(2) + 2.0*G*sym(grad(Du))
    a = inner(grad(tu),Dsigma)*dx
    L = - inner(tu, n*P)*ds(2) - inner(tu, n*P)*ds(4)
    bcs = [
        DirichletBC(V, Constant((0.0,0.0)), facetids,3)
    ]
    # Solve the problem
    solve(a==L, u, bcs=bcs)

    W = assemble( inner(grad(u),Dsigma)*dx )
    G = assemble( inner(tu,n*P)*ds(2) + inner(tu,n*P)*ds(4) )
    # Delete the meshes
    shutil.rmtree(tmpdir)
    return {
        'W':W,
        'G':G,
        'x':mesh.coordinates(),
        'u':u.vector().get_local(),
    }
    
sim(0.1,0.1,0.2,0.3,1.0)
