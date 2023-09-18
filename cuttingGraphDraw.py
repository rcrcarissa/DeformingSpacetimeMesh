import pickle
import numpy as np
import vtk
from vtkmodules.util import numpy_support


def vtuWriteLines(fpath: str, pts: np.array, lns: np.array):
    points = vtk.vtkPoints()
    for pt in pts:
        points.InsertNextPoint(pt[0], pt[1], pt[2])

    lines = vtk.vtkCellArray()
    for l in lns:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, l[0])
        line.GetPointIds().SetId(1, l[1])
        lines.InsertNextCell(line)

    data_save = vtk.vtkPolyData()
    data_save.SetPoints(points)
    data_save.SetLines(lines)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(data_save)
    writer.SetFileName(fpath)
    writer.Write()


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("real_faces.pickle", "rb") as f:
        real_faces = pickle.load(f)
    num_cutting_edges = len(real_faces)
    cutting_graph_lower = np.zeros((num_cutting_edges, 2), dtype=int)
    cutting_graph_upper = np.zeros((num_cutting_edges, 2), dtype=int)
    for i in range(num_cutting_edges):
        cutting_graph_lower[i, 0] = real_faces[i][0]
        cutting_graph_lower[i, 1] = real_faces[i][1]
        cutting_graph_upper[i, 0] = real_faces[i][2]
        cutting_graph_upper[i, 1] = real_faces[i][3]

    vtuWriteLines('cutting_graph_lower.vtp',
                  np.concatenate([nodes_coor, np.array([np.zeros(nodes_coor.shape[0])]).T], axis=1),
                  cutting_graph_lower)
    vtuWriteLines('cutting_graph_upper.vtp',
                  np.concatenate([nodes_coor, np.array([np.zeros(nodes_coor.shape[0])]).T], axis=1),
                  cutting_graph_upper)
