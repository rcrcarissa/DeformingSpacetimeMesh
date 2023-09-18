import pickle
import numpy as np



def orientTet(tet, nodes_coor):
    tet_n0 = np.append(nodes_coor[tet[0][0]], int(tet[0][1] == 'u'))
    tet_n1 = np.append(nodes_coor[tet[1][0]], int(tet[1][1] == 'u'))
    tet_n2 = np.append(nodes_coor[tet[2][0]], int(tet[2][1] == 'u'))
    tet_n3 = np.append(nodes_coor[tet[3][0]], int(tet[3][1] == 'u'))
    p0 = tet_n1 - tet_n0
    p1 = tet_n2 - tet_n0
    p2 = tet_n3 - tet_n0
    sign = np.dot(np.cross(p0, p1), p2)
    if sign > 0:
        return tet
    return [tet[0], tet[2], tet[1], tet[3]]


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)

    with open("tetrahedra.pickle", "rb") as f:
        tetrahedralization = pickle.load(f)
    with open("mergedTri_tetrahedra.pickle", "rb") as f:
        tetrahedralization_mergedTri = pickle.load(f)
    with open("boundaryTris_tetrahedra.pickle", "rb") as f:
        tetrahedralization_boundaryTri = pickle.load(f)
    with open("boundaryTris_tetrahedra_upper.pickle", "rb") as f:
        tetrahedralization_boundaryTri_upper = pickle.load(f)
    tetrahedra = []
    for tets in tetrahedralization.values():
        for tet in tets:
            tetrahedra.append(orientTet(tet, nodes_coor))
    for tets in tetrahedralization_mergedTri.values():
        for tet in tets:
            tetrahedra.append(orientTet(tet, nodes_coor))
    for tets in tetrahedralization_boundaryTri.values():
        for tet in tets:
            tetrahedra.append(orientTet(tet, nodes_coor))
    for tets in tetrahedralization_boundaryTri_upper.values():
        for tet in tets:
            tetrahedra.append(orientTet(tet, nodes_coor))
    with open("all_tets.pickle", "wb") as f:
        pickle.dump(tetrahedra, f)
