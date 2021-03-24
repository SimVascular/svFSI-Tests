# Header
import numpy as np
import vtk
import os

class bcolors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    ULINE   = '\033[4m'

def setDomainID(FileName):
    
    print(bcolors.HEADER + "="*60 + bcolors.ENDC)
    print("")
    print(bcolors.OKBLUE + "Loading mesh file <--- %s" % (FileName) + bcolors.ENDC)
    print("")

    mshReader = vtk.vtkXMLUnstructuredGridReader()
    mshReader.SetFileName(FileName)
    mshReader.Update()

    msh = vtk.vtkUnstructuredGrid()
    msh = mshReader.GetOutput()

    msh_npts = msh.GetNumberOfPoints()
    msh_ncel = msh.GetNumberOfCells()

    domain_id = np.ones((msh_ncel,1))
    for icel in range(0, msh_ncel):
        cur_cel = msh.GetCell(icel)
        cur_cel_npts = cur_cel.GetNumberOfPoints()
        pts_cel = cur_cel.GetPointIds()
        cnt = 0
        for i in range(0, cur_cel_npts):
            xp = msh.GetPoint(pts_cel.GetId(i))
            if xp[1] < 1.001:
                cnt = cnt + 1
        if cnt == cur_cel_npts:
            domain_id[icel] = 2

    path, fname = os.path.split(FileName)

    # Write domain ID to file
    fname = "%s/domain_info.dat" % (path)
    
    print("")
    print(bcolors.OKBLUE + "Writing domain file ---> %s" % (fname) + bcolors.ENDC)
    print("")
    print(bcolors.HEADER + "="*60 + bcolors.ENDC)
    np.savetxt(fname, domain_id, delimiter='', fmt='%d')

if __name__ == '__main__':
    setDomainID("h1.00/mesh-complete.mesh.vtu")

