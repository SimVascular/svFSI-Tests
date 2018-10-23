# Header
import sympy as sp
import numpy as np
import vtk as vtk
import os as os
sp.init_printing(use_unicode = True)

class bcolors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    ULINE   = '\033[4m'

# Define symbols
d        = sp.symbols('d')        # dimension, d=2, 3
x, y, z  = sp.symbols('x y z')    # coordinate axes
t, T0    = sp.symbols('t T0')     # time, time period
omega, c = sp.symbols('omega, c') # angular frequency (displacement)
rho      = sp.symbols('rho')      # density of the structure
k        = sp.symbols('k')        # bulk modulus / penalty parameter
mu       = sp.symbols('mu')       # shear modulus (Neo-Hookean model)
C1, C2   = sp.symbols('C1 C2')    # material constants (Mooney-Rivlin model)
M0, L0   = sp.symbols('M0 L0')    # mass, length scales
beta     = sp.symbols('beta')     # angular frequency (pressure)

# Input parameters
nt      = 21                  # num. time steps
mshFile = 'mesh/N004/mesh-complete.mesh.vtu'

#=======================================================================
# Define displacement field
u1 = c * L0 * (t/T0)**2 * sp.sin(omega*y) * sp.sin(omega*z)
u2 = 0
u3 = 0
u  = sp.Matrix([u1, u2, u3])

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Displacement field, u: " + bcolors.ENDC)
print("")
sp.pprint(u)
print("")

#=======================================================================
# Compute the deformation gradient tensor
F = sp.eye(3) + sp.Matrix([[sp.diff(u1, x), sp.diff(u1, y), sp.diff(u1, z)], \
                           [sp.diff(u2, x), sp.diff(u2, y), sp.diff(u2, z)], \
                           [sp.diff(u3, x), sp.diff(u3, y), sp.diff(u3, z)]])

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Deformation tensor, F: " + bcolors.ENDC)
print("")
sp.pprint(F)
print("")

# Jacobian
J = F.det()
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Jacobian, J: " + bcolors.ENDC)
print("")
sp.pprint(J)
print("")

# Cauchy-Green deformation tensor and its inverse
C  = F.T * F
Ci = C**-1

#=======================================================================
# Now construct 1st & 2nd Piola-Kirchhoff Stress tensor
# Volumetric 2nd Piola-Kirchhoff stress tensor:

p = (M0/L0) * (t**2/T0**4) * sp.sin(beta*x) * sp.sin(beta*y) * sp.sin(beta*z)

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Pressure, p: " + bcolors.ENDC)
print("")
sp.pprint(p)
print("")

Svol = p*J*Ci

# Isochoric 2nd Piola-Kirhhoff stress tensor:
#    Neo-Hookean
Sbar = mu*sp.eye(3)

# TrS := S:C
TrS  = sp.tensorcontraction(sp.tensorcontraction(sp.tensorproduct(Sbar, C), (1,2)), (0,1))

Siso = J**(-2/d) * (Sbar - (1/d * TrS) * Ci)

PK1 = F*(Svol + Siso)

#=======================================================================
# Structure internal force, Div.(PK1)
DivP = sp.diff(PK1.col(0), x)
DivP = DivP.col_insert(1, sp.diff(PK1.col(1),y))
DivP = DivP.col_insert(2, sp.diff(PK1.col(2),z))
fs = sp.Matrix([sum(DivP.row(0)), sum(DivP.row(1)), sum(DivP.row(2))]) 

#=======================================================================
# Inertial force, ddot(u)
fi = sp.Matrix([sp.diff(u1, t, t), sp.diff(u2, t, t), sp.diff(u3, t, t)])

#=======================================================================
# Body force, fb
fb  = fi - (1/rho)*fs

fb_fn = sp.lambdify([t, T0, M0, L0, omega, beta, c, d, rho, mu, x, y, z], \
    fb, "numpy")

#=======================================================================
# Define time, domain and material parameters
# Set time parameters: t, T0
T0    = 0.001
tr    = np.linspace(0, T0, nt)

M0    = 1.0
L0    = 1.0

# Set parameters: d, omega
d     = 3
c     = 1.0
omega = 0.1*np.pi
beta  = 0.2*np.pi

# Set material parameters:
# Density, g/cm3
rho   = 1.0

# Values chosen based on E = 3.0E4 dyn/cm2, nu = 0.5
# Isochoric: NeoHookean (mu, dyn/cm2), Mooney-Rivilin (C1, C2)
mu    = 10**4

#=======================================================================
# Read VTU mesh file and load point coordinates
mshReader = vtk.vtkXMLUnstructuredGridReader()
mshReader.SetFileName(mshFile)
mshReader.Update()

msh = vtk.vtkUnstructuredGrid()
msh = mshReader.GetOutput()

msh_npts = msh.GetNumberOfPoints()
msh_pts  = np.zeros((msh_npts,d))
for ipt in range(0, msh_npts):
    msh_pts[ipt,:] = msh.GetPoint(ipt)

#=======================================================================
# Evaluate body force over the domain and write to file
bff = vtk.vtkDoubleArray()
bff.SetNumberOfComponents(3)
bff.Allocate(msh_npts,128)
bff.SetNumberOfTuples(msh_npts*3)
bff.SetName("FB")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Writing csv and vtu files from t=%.4f to t=%.4f" % \
    (tr[0], tr[nt-1]) + bcolors.ENDC)
print("")
os.system('mkdir -p csv  vtu')
i = 0
for time in tr:
    fbg = fb_fn(time, T0, M0, L0, omega, beta, c, d, rho, mu, \
        msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
    i = i + 1
    for j in range(0, d):
        fname = "csv/fb%d_%02d.csv" % (j+1, i)
        np.savetxt(fname, fbg[j,0,:], delimiter=',')
    for j in range(0, msh_npts):
        rtmp = fbg[:,0,j]
        bff.SetTuple3(j, rtmp[0], rtmp[1], rtmp[2])
    fname = "vtu/fb_%02d.vtu" % (i)
    msh.GetPointData().AddArray(bff)
    mshWrite = vtk.vtkXMLUnstructuredGridWriter()
    mshWrite.SetInputData(msh)
    mshWrite.SetFileName(fname)
    mshWrite.Write()
print(bcolors.HEADER + "="*80 + bcolors.ENDC)

#=======================================================================
# EOF

