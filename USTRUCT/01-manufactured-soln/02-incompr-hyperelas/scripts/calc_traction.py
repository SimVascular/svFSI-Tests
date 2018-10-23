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

def genFaceBC(FileName, nv, flag) :

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

    #=======================================================================
    # Define displacement field
    u1 = c * L0 * (t/T0)**2 * sp.sin(omega*y) * sp.sin(omega*z)
    u2 = 0
    u3 = 0
    u  = sp.Matrix([u1, u2, u3])

    if flag:
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

    if flag:
        print(bcolors.HEADER + "="*80 + bcolors.ENDC)
        print("")
        print(bcolors.OKBLUE + "Deformation tensor, F: " + bcolors.ENDC)
        print("")
        sp.pprint(F)
        print("")

    # Jacobian
    J = F.det()
    if flag:
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
    # Now construct 1st and 2nd Piola-Kirchhoff stress tensors
    # Volumetric 2nd Piola-Kirchhoff stress tensor:

    p = (M0/L0) * (t**2/T0**4) * sp.sin(beta*x) * sp.sin(beta*y) * sp.sin(beta*z)

    if flag:
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

    PK1 = F * (Svol + Siso)

    #=======================================================================
    # Traction vector
    hv = PK1 * nv

#    h  = hv[0]*nv[0] + hv[1]*nv[1] + hv[2]*nv[2]
    h_fn = sp.lambdify([t, T0, M0, L0, omega, beta, c, d, mu, x, y, z], \
        hv, "numpy")

    #=======================================================================
    # Define time, domain and material parameters
    # Set time parameters: t, T0
    nt    = 21
    T0    = 0.001
    tr    = np.linspace(0, T0, nt)

    M0    = 1.0
    L0    = 1.0

    # Set parameters: d, c, omega, gamma
    d     = 3
    c     = 1.0
    omega = 0.1*np.pi
    beta  = 0.2*np.pi

    # Values chosen based on E = 3.0E4 dyn/cm2, nu = 0.5
    # Isochoric: NeoHookean (mu, dyn/cm2), Mooney-Rivilin (C1, C2)
    mu    = 10**4

    #=======================================================================
    vtpReader = vtk.vtkXMLPolyDataReader()
    vtpReader.SetFileName(FileName)
    vtpReader.Update()

    pdata = vtk.vtkPolyData()
    pdata = vtpReader.GetOutput()

    pdata_npts = pdata.GetNumberOfPoints()
    pdata_pts  = np.zeros((pdata_npts,d))
    pdata_nid  = np.zeros((pdata_npts,1),int)
    pdata_nid  = pdata.GetPointData().GetArray('GlobalNodeID')
    for ipt in range(0, pdata_npts):
        pdata_pts[ipt,:] = pdata.GetPoint(ipt)

    path, fname = os.path.split(FileName)
    fhdr, ext   = os.path.splitext(fname)

    # Write face nodes' global node Id to a file
    fname = "csv/bc_%s_nodeid.csv" % (fhdr)
    np.savetxt(fname, pdata_nid, delimiter=',', fmt='%d')

    # Evaluate traction over the face and write to file
    print(bcolors.HEADER + "="*80 + bcolors.ENDC)
    print("")
    print(bcolors.OKBLUE + "Writing BC data for face %s from t=%.4f to t=%.4f" % \
        (fhdr, tr[0], tr[nt-1]) + bcolors.ENDC)
    print("")
    i = 0
    for time in tr:
        hg = h_fn(time, T0, M0, L0, omega, beta, c, d, mu, \
            pdata_pts[:,0], pdata_pts[:,1], pdata_pts[:,2])        
        i = i + 1
#        for j in range(0, np.size(hg,0)):
#            fname = "csv/bc_%s_h%d_%02d.csv" % (fhdr, j+1, i)
#            np.savetxt(fname, hg[j,0,:])
        fname = "csv/bc_%s_h1_%02d.csv" % (fhdr, i)
        np.savetxt(fname, hg[0,0])
#        fname = "csv/bc_%s_h2_%02d.csv" % (fhdr, i)
#        np.savetxt(fname, hg[1,0])
        fname = "csv/bc_%s_h3_%02d.csv" % (fhdr, i)
        np.savetxt(fname, hg[2,0])


if __name__ == '__main__':
#    nv = sp.Matrix([-1,0,0])
#    genFaceBC("mesh/N004/mesh-surfaces/X0.vtp", nv, True)
    
#    nv = sp.Matrix([1,0,0])
#    genFaceBC("mesh/N004/mesh-surfaces/X1.vtp", nv, False)
    
#    nv = sp.Matrix([0,1,0])
#    genFaceBC("mesh/N004/mesh-surfaces/Y0.vtp", nv, True)
    
#    nv = sp.Matrix([0,1,0])
#    genFaceBC("mesh/N004/mesh-surfaces/Y1.vtp", nv, False)
    
    nv = sp.Matrix([0,0,-1])
    genFaceBC("mesh/N004/mesh-surfaces/Z0.vtp", nv, True)
    
    nv = sp.Matrix([0,0,1])
    genFaceBC("mesh/N004/mesh-surfaces/Z1.vtp", nv, False)

