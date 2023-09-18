'''
This code generates a vtu file for tetrahedral mesh.
Author: Congrong Ren
Last Modified: Mar 10, 2023
'''
import pickle
import numpy as np
import vtk
from vtkmodules.util import numpy_support

upper_height = 0.1


def vtuWriteTetra(fpath: str, pts: np.array, cls: np.array, vector_fields: dict = {}):
    points = vtk.vtkPoints()
    for pt in pts:
        points.InsertNextPoint(pt[0], pt[1], pt[2])

    data_save = vtk.vtkUnstructuredGrid()
    data_save.SetPoints(points)
    for c in cls:
        cellId = vtk.vtkIdList()
        for i in c:
            cellId.InsertNextId(i)
        data_save.InsertNextCell(vtk.VTK_TETRA, cellId)

    pd = data_save.GetCellData()
    for (k, v) in vector_fields.items():
        vtk_array = numpy_support.numpy_to_vtk(v)
        vtk_array.SetName(k)
        pd.AddArray(vtk_array)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(data_save)
    writer.SetFileName(fpath)
    writer.Write()


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    vtu_fname = "all.vtu"
    with open("tetrahedra.pickle", "rb") as f:
        tetrahedralization = pickle.load(f)
    with open("mergedTri_tetrahedra.pickle", "rb") as f:
        tetrahedralization_mergedTri = pickle.load(f)
    # with open("boundaryTris_tetrahedra.pickle", "rb") as f:
    #     tetrahedralization_boundaryTri = pickle.load(f)
    # with open("boundaryTris_tetrahedra_upper.pickle", "rb") as f:
    #     tetrahedralization_boundaryTri_upper = pickle.load(f)
    tetrahedra = []
    for tets in tetrahedralization.values():
        tetrahedra += tets
    for tets in tetrahedralization_mergedTri.values():
        tetrahedra += tets
    # for tets in tetrahedralization_boundaryTri.values():
    #     tetrahedra += tets
    # for tets in tetrahedralization_boundaryTri_upper.values():
    #     tetrahedra += tets
    # with open("all_tets.pickle", "rb") as f:
    #     tetrahedra = pickle.load(f)
    # tetrahedra = [[(923, 'l'), (924, 'u'), (925, 'u'), (926, 'u')],[(922, 'l'), (923, 'l'), (924, 'u'), (925, 'u')]]
    i = 0
    related_nodes_coor = np.concatenate((nodes_coor, np.zeros((nodes_coor.shape[0], 1))), axis=1)
    temp = np.concatenate((nodes_coor, upper_height * np.ones((nodes_coor.shape[0], 1))), axis=1)
    related_nodes_coor = np.concatenate((related_nodes_coor, temp), axis=0)
    related_cells = []
    related_node_idx = []
    related_node_layer = []
    for tet in tetrahedra:
        cell = []
        node_layers = []
        node_idx = []
        for n in tet:
            if n[1] == 'l':
                cell.append(n[0])
                node_layers.append(0)
            else:
                cell.append(n[0] + nodes_coor.shape[0])
                node_layers.append(1)
            node_idx.append(n[0])
        related_cells.append(cell)
        related_node_layer.append(node_layers)
        related_node_idx.append(node_idx)

    vtuWriteTetra(vtu_fname, related_nodes_coor, related_cells,
                  {'node_idx': np.array(related_node_idx), 'node_layer': np.array(related_node_layer)})
