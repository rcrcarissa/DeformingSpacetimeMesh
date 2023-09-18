import numpy as np
import vtk
from vtkmodules.util import numpy_support
import pickle


def vtuWriteTriangleAttr(fpath: str, pts: np.array, cls: np.array, scalar_fields: dict = {}, vector_fields: dict = {}):
    vtk_position = numpy_support.numpy_to_vtk(pts)
    points = vtk.vtkPoints()
    triangle = vtk.vtkTriangle()
    cells = vtk.vtkCellArray()

    for cell in cls:
        for i in range(3):
            triangle.GetPointIds().SetId(i, cell[i])
        cells.InsertNextCell(triangle)

    vtu = vtk.vtkUnstructuredGrid()
    points.SetData(vtk_position)
    vtu.SetPoints(points)
    vtu.SetCells(vtk.VTK_TRIANGLE, cells)

    # setup cell values
    pd = vtu.GetCellData()
    for i, (k, v) in enumerate(scalar_fields.items()):
        vtk_array = numpy_support.numpy_to_vtk(v)
        vtk_array.SetName(k)
        if i == 0:
            pd.SetScalars(vtk_array)
        else:
            pd.AddArray(vtk_array)

    for i, (k, v) in enumerate(vector_fields.items()):
        vtk_array = numpy_support.numpy_to_vtk(v)
        vtk_array.SetName(k)
        if i == 0:
            pd.SetVectors(vtk_array)
        else:
            pd.AddArray(vtk_array)

    # # setup point values
    # pd = vtu.GetPointData()
    # for i, (k, v) in enumerate(scalar_fields.items()):
    #     vtk_array = numpy_support.numpy_to_vtk(v)
    #     vtk_array.SetName(k)
    #     if i == 0:
    #         pd.SetScalars(vtk_array)
    #     else:
    #         pd.AddArray(vtk_array)
    # 
    # for i, (k, v) in enumerate(vector_fields.items()):
    #     vtk_array = numpy_support.numpy_to_vtk(v)
    #     vtk_array.SetName(k)
    #     if i == 0:
    #         pd.SetVectors(vtk_array)
    #     else:
    #         pd.AddArray(vtk_array)

    writer = vtk.vtkXMLDataSetWriter()
    writer.SetFileName(fpath)
    writer.SetInputData(vtu)
    writer.Write()


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("independent_patches_cellID.pickle", "rb") as f:
        indSys = pickle.load(f)
    with open("comp_lower.pickle", "rb") as f:
        comp_lower = pickle.load(f)
    with open("comp_upper.pickle", "rb") as f:
        comp_upper = pickle.load(f)
    cell_lower_segIDs = np.zeros(cells.shape[0])
    cell_upper_segIDs = np.zeros(cells.shape[0])
    for i, sys in enumerate(indSys):
        for tri in sys[0]:
            cell_lower_segIDs[tri] = i
        comp_lower.remove(sys[0])
        for tri in sys[1]:
            cell_upper_segIDs[tri] = i
        comp_upper.remove(sys[1])
    curr_len = len(indSys)
    for i, comp in enumerate(comp_lower):
        for tri in comp:
            cell_lower_segIDs[tri] = curr_len + i
    curr_len = len(indSys) + len(comp_lower)
    for i, comp in enumerate(comp_upper):
        for tri in comp:
            cell_upper_segIDs[tri] = curr_len + i
    vtuWriteTriangleAttr("SyntheticMeshSegID.vtu",
                         np.concatenate([nodes_coor, np.array([np.zeros(nodes_coor.shape[0])]).T], axis=1), cells,
                         {"cell_lower_segIDs": cell_lower_segIDs, "cell_upper_segIDs": cell_upper_segIDs})
