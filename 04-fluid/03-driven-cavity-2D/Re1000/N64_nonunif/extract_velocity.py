import numpy as np 
import vtk
from vtk.util import numpy_support as vtknp

def readvtu(fname):
    # read the ASCII version
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()


    output = reader.GetOutput()

    # number of points
    NoP = output.GetNumberOfPoints()

    # number of cells
    NoC = output.GetNumberOfCells()

    # coordinates of points
    temp_vtk_array = output.GetPoints().GetData()
    coord = vtknp.vtk_to_numpy(temp_vtk_array)

    # displacement vector defined on each point
    temp_vtk_array = output.GetPointData().GetVectors('Velocity')
    u = vtknp.vtk_to_numpy(temp_vtk_array)

    conn = [ ]
    # cell connectivity
    for i in range(NoC):
        cell = output.GetCell(i)
        npts = cell.GetNumberOfPoints()
        connt = [cell.GetPointId(j) for j in range(npts)]
        conn.append(connt)

    conn = np.asarray(conn)

    return u, coord, conn 

fname = "./16-procs/result_2000.vtu"
vel, coord, conn = readvtu(fname)

indx = np.where(coord[:,1]==0.5)[0]
xc = coord[indx,0]
indx = indx[np.argsort(xc)]

fid = open('v_x0.5_Re1000.dat','w')
fid.write("variables = x, v \n")
for i in indx:
    fid.write("{}  {}\n".format(coord[i,0],vel[i,1]))

fid.close()

indx = np.where(coord[:,0]==0.5)[0]
yc = coord[indx,1]
indx = indx[np.argsort(yc)]

fid = open('u_y0.5_Re1000.dat','w')
fid.write("variables = y, u \n")
for i in indx:
    fid.write("{}  {}\n".format(coord[i,1],vel[i,0]))

fid.close()
