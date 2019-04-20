from fenics import *
from SimDataDB import SimDataDB
from draw_fracture import draw_fracture
import os, tempfile, shutil
import numpy as np

# gmshloc = "/Applications/Gmsh.app/Contents/MacOS/gmsh"
gmshloc = "gmsh"
dbloc = "femdata.db"
geo = os.path.abspath("crack_gen.geo")

# Some global properties
E = 1.0
nu= 0.25
P = 1

pixel_res = 32
xs = np.linspace(-1,1,pixel_res)
# pix_xy = np.empty(pixel_res**2,2)
# for x in xs:
#     for y in xs:
#         pix_xy[i,:]=x,y
sdb = SimDataDB('pix2.sqlite')

@sdb.Decorate('static',
             [('x1','float'),('y1','float'),('x2','float'),('y2','float'),('clscale','float')],
             [('W','float'),('W_P','float'),('G_c','float'),('V','float'),('pic','array')])
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
    eps = sym(grad(u))
    w = (K-2.0/3.0*G)/2.0*tr(eps)**2 + G * inner(eps,eps)
    Dsigma = (K-2.0/3.0*G)*tr(sym(grad(Du)))*Identity(2) + 2.0*G*sym(grad(Du))
    a = inner(grad(tu),Dsigma)*dx
    L = - inner(tu, n*P)*ds(6) - inner(tu, n*P)*ds(7)
    # TODO I want the left side to be roller
    bcs = [
        DirichletBC(V, Constant((0.0,0.0)), facetids,2),
        DirichletBC(V, Constant((0.0,0.0)), facetids,3),
        DirichletBC(V, Constant((0.0,0.0)), facetids,4),               
        DirichletBC(V.sub(0), Constant(0.0), facetids,5),
    ]

    # Solve the problem
    solve(a==L, u, bcs=bcs)

    # Post process
    # Draw a picture
    pic = np.empty((pixel_res,pixel_res,3),dtype=np.float32)
    frac_pic = draw_fracture(x1,y1,x2,y2, (pixel_res,pixel_res))
    for i,x in enumerate(xs):
        for j,y in enumerate(xs):
            uN = u(x,y)
            pic[j,i,0] = uN[0]
            pic[j,i,1] = uN[1]
            pic[j,i,2] = frac_pic.getpixel((i,j))
    # Strain energy
    W = assemble( w*dx )
    # Fracture energy, equal to length
    G_c = assemble( dot(n,n)*ds(6) )
    # Virtual work on fracture
    W_P = assemble( inner(u,n*P)*ds(6) + inner(u,n*P)*ds(7) ) 
    Volume = assemble(inner(u,n)*ds(6)+inner(u,n)*ds(7))
    
    # Delete the meshes
    shutil.rmtree(tmpdir)
    
    # Return the database
    return {
        'W':W,
        'G_c':G_c,
        'W_P':W_P,
        'V':Volume,
        'pic':pic
#         'x':mesh.coordinates(),
#         'u':u.vector().get_local(),
    }
    
for tip_x in np.linspace(-0.7,0.9,6):
    for tip_y in np.linspace(-0.5,0.5,6):
        for crook_x in np.linspace(-0.9, tip_x-0.1, 6):
            for crook_y in np.linspace( -0.5,0.5, 6):
                try:
                    sim(crook_x,crook_y, tip_x, tip_y,1.0)
                except Exception as e:
                    print(e)

