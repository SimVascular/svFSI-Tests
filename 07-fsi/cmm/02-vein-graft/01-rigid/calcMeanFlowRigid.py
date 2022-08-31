# Header
import numpy as np
import vtk

# User-defined inputs
RESULTS_FOLDER  = "48-procs"
START_TIME_STEP = 3000
END_TIME_STEP   = 4000
DIFF_TIME_STEP  = 25
TIME_STEP_SIZE  = 0.000441
WALLS_FILE      = "../mesh/mesh-surfaces/walls_combined.vtp"

#----------------------------------------------------------------------
def loadVTU(fileName):

    print ("   Loading vtu file   <---   {}".format(fileName))
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()
    vtuMesh = vtuReader.GetOutput()

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getMeshData(msh, keyWord):
    msh_data = vtk.vtkDoubleArray()
    msh_data = msh.GetPointData().GetArray(keyWord)

    return msh_data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getWallNodes(FileName):
    modelReader = vtk.vtkXMLPolyDataReader()
    modelReader.SetFileName(FileName)
    modelReader.Update()

    model = vtk.vtkPolyData()
    model = modelReader.GetOutput()

    model_npts = model.GetNumberOfPoints()
    model_ID = model.GetPointData().GetArray('GlobalNodeID')
    id_list = np.zeros((model_npts,1), dtype='int32')
    for i in range(0, model_npts):
        id_list[i] = int(model_ID.GetTuple1(i))
    return id_list
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getSurfaceData(msh, srf_ids, keyWord):
    msh_data = getMeshData(msh, keyWord)
    num_comp = msh_data.GetNumberOfComponents()

    srf_nno = np.size(srf_ids)
    if num_comp == 1:
        srf_data = np.zeros((srf_nno,))
    else:
        srf_data = np.zeros((srf_nno,num_comp))
    for ipt in range(0, srf_nno):
        h = msh_data.GetTuple(srf_ids[ipt]-1)
        if num_comp == 1:
            srf_data[ipt] = h[0]
        else:
            for j in range(0, num_comp):
                srf_data[ipt,j] = h[j]

    return srf_data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def createVTKDataArray(dType, numComp, numTupls, dName):

    if dType == "double":
        D = vtk.vtkDoubleArray()
    elif dType == "int":
        D = vtk.vtkIntArray()

    D.SetNumberOfComponents(numComp)
    D.Allocate(numTupls)
    D.SetNumberOfTuples(numTupls)
    D.SetName(dName)

    return D
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def vtpInitWriter(fwall):
    modelReader = vtk.vtkXMLPolyDataReader()
    modelReader.SetFileName(fwall)
    modelReader.Update()

    model = vtk.vtkPolyData()
    model = modelReader.GetOutput()

    return model
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def vtpAddScalarArray(model, dArr, dName):
    numTupls = np.size(dArr)
    vtkD = createVTKDataArray("double", 1, numTupls, dName)
    for i in range(0, numTupls):
        vtkD.SetTuple1(i, dArr[i])

    model.GetPointData().AddArray(vtkD)

    return model
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def vtpAddVectorArray(model, dArr, dName):
    numTupls = np.shape(dArr)[0]
    vtkD = createVTKDataArray("double", 3, numTupls, dName)
    for i in range(0, numTupls):
        vtkD.SetTuple3(i, dArr[i,0], dArr[i,1], dArr[i,2])

    model.GetPointData().AddArray(vtkD)

    return model
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def vtpWriteToFile(model, fout):
    modelWrite = vtk.vtkXMLPolyDataWriter()
    modelWrite.SetInputData(model)
    modelWrite.SetFileName(fout)
    modelWrite.Write()

    return
#----------------------------------------------------------------------

if __name__ == '__main__':
    srcdir = RESULTS_FOLDER
    nstart = START_TIME_STEP
    nend   = END_TIME_STEP
    nfreq  = DIFF_TIME_STEP
    dt     = TIME_STEP_SIZE
    fwall  = WALLS_FILE

    wall_ids = getWallNodes(fwall)
    wall_nno = np.size(wall_ids)
    mean_P   = np.zeros((wall_nno), dtype='float64')
    mean_H   = np.zeros((wall_nno,3), dtype='float64')
    mean_WSS = np.zeros((wall_nno,3), dtype='float64')

    nframe = int((nend - nstart)/nfreq) + 1
    for i in range(0, nframe):
        ntime = nstart + i*nfreq
        time  = float(ntime)*dt
        if (ntime < 100):
            fname = "%s/result_%03d.vtu" %(srcdir, ntime)
        else:
            fname = "%s/result_%d.vtu" %(srcdir, ntime)
        print ("Reading file    <-----   {}".format(fname))
        vtuData  = loadVTU(fname)
        mean_P   = mean_P   + getSurfaceData(vtuData, wall_ids, 'Pressure')
        mean_H   = mean_H   + getSurfaceData(vtuData, wall_ids, 'Traction')
        mean_WSS = mean_WSS + getSurfaceData(vtuData, wall_ids, 'WSS')

    mean_P   = mean_P   / float(nframe)
    mean_H   = mean_H   / float(nframe)
    mean_WSS = mean_WSS / float(nframe)

    model = vtpInitWriter(fwall)
    vtpAddScalarArray(model, mean_P, "Pressure")
    vtpAddVectorArray(model, mean_H, "Traction")
    vtpAddVectorArray(model, mean_WSS, "WSS")

    fout = 'avg_flow_rigid.vtp'
    print ("Writing average flow file    ----->   {}".format(fout))
    vtpWriteToFile(model, fout)

#=======================================================================
# EOF


