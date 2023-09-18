'''
This code draws one independent system and its correlated real face(s).
Author: Congrong Ren
Last Modified: Dec 29, 2022
'''
import pickle
import numpy as np
import vtk
from vtkmodules.util import numpy_support

nphi = 16
upper_height = 0.1


def orientTet(tet, nodes_coor):
    p0 = nodes_coor[tet[1]] - nodes_coor[tet[0]]
    p1 = nodes_coor[tet[2]] - nodes_coor[tet[0]]
    p2 = nodes_coor[tet[3]] - nodes_coor[tet[0]]
    sign = np.dot(np.cross(p0, p1), p2)
    if sign > 0:
        return tet
    return [tet[0], tet[2], tet[1], tet[3]]


def vtuWriteTetra(fpath: str, pts: np.array, cls: np.array, scalar_fields: dict = {}):
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

    pd = data_save.GetPointData()
    for i, (k, v) in enumerate(scalar_fields.items()):
        vtk_array = numpy_support.numpy_to_vtk(v)
        vtk_array.SetName(k)
        if i == 0:
            pd.SetScalars(vtk_array)
        else:
            pd.AddArray(vtk_array)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(data_save)
    writer.SetFileName(fpath)
    writer.Write()


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    with open("dataGen/values.pickle", "rb") as f:
        values = pickle.load(f)
    with open("dataGen/values_denoise.pickle", "rb") as f:
        values_denoise = pickle.load(f)

    vtu_fname = "straightMesh.vtu"
    tet_nodes_coor = np.zeros((2 * num_nodes, 3))
    tet_nodes_coor[:num_nodes, :-1] = nodes_coor
    tet_nodes_coor[num_nodes:, :-1] = nodes_coor
    tet_nodes_coor[num_nodes:, -1] = np.ones(num_nodes) * upper_height
    tetrahedra = []
    for c in cells:
        sorted_c = np.sort(c)
        t = [sorted_c[0], sorted_c[1], sorted_c[2], sorted_c[2] + num_nodes]
        tetrahedra.append(orientTet(t, tet_nodes_coor))
        t = [sorted_c[0], sorted_c[0] + num_nodes, sorted_c[1] + num_nodes, sorted_c[2] + num_nodes]
        tetrahedra.append(orientTet(t, tet_nodes_coor))
        t = [sorted_c[0], sorted_c[1], sorted_c[1] + num_nodes, sorted_c[2] + num_nodes]
        tetrahedra.append(orientTet(t, tet_nodes_coor))
    vtuWriteTetra(vtu_fname, tet_nodes_coor, tetrahedra)

    vtu_fname = "straightVolRen.vtu"
    tet_nodes_coor = np.zeros(((nphi + 1) * num_nodes, 3))
    for phi in range(nphi + 1):
        tet_nodes_coor[num_nodes * phi:num_nodes * (phi + 1), :-1] = nodes_coor
        tet_nodes_coor[num_nodes * phi:num_nodes * (phi + 1), -1] = np.ones(num_nodes) * upper_height * phi
    tetrahedra = []
    for phi in range(nphi):
        for c in cells:
            sorted_c = np.sort(c)
            t = [sorted_c[0] + phi * num_nodes, sorted_c[1] + phi * num_nodes, sorted_c[2] + phi * num_nodes,
                 sorted_c[2] + (phi + 1) * num_nodes]
            tetrahedra.append(orientTet(t, tet_nodes_coor))
            t = [sorted_c[0] + phi * num_nodes, sorted_c[0] + (phi + 1) * num_nodes,
                 sorted_c[1] + (phi + 1) * num_nodes, sorted_c[2] + (phi + 1) * num_nodes]
            tetrahedra.append(orientTet(t, tet_nodes_coor))
            t = [sorted_c[0] + phi * num_nodes, sorted_c[1] + phi * num_nodes, sorted_c[1] + (phi + 1) * num_nodes,
                 sorted_c[2] + (phi + 1) * num_nodes]
            tetrahedra.append(orientTet(t, tet_nodes_coor))
    vol_values = np.zeros((nphi + 1) * num_nodes)
    vol_values_denoise = np.zeros((nphi + 1) * num_nodes)
    for phi in range(nphi + 1):
        vol_values[phi * num_nodes:(phi + 1) * num_nodes] = values[:, phi % nphi]
        vol_values_denoise[phi * num_nodes:(phi + 1) * num_nodes] = values_denoise[:, phi % nphi]
    vtuWriteTetra(vtu_fname, tet_nodes_coor, np.array(tetrahedra),
                  {'values': vol_values, 'values_denoise': vol_values_denoise})
