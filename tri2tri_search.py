import numpy as np
import pickle

if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    # node_to_triangle hash and node_to_edge hash
    triangles_of_node = {}
    edges_of_node = {}
    for i in range(num_nodes):
        triangles_of_node[i] = []
        edges_of_node[i] = []
    for cl in cells:
        for pt in cl:
            triangles_of_node[pt].append(cl)
            adjacent_nodes = cl[cl != pt]
            edges_of_node[pt].append([pt, adjacent_nodes[0]])
            edges_of_node[pt].append([pt, adjacent_nodes[1]])
            # Note: "edges_of_node" may contain duplicate edges for a node
    tri2tri_pairs = []
    tri2edge_idx = []
    tri2path_idx = []
    tri2node_idx = []
    for i, cl in enumerate(cells):
        next_pts = [next_nodes[pt] for pt in cl]
        if len(set(next_pts)) == 3:
            find_cell = [idx for idx, cell in enumerate(triangles_of_node[next_pts[2]]) if set(next_pts) == set(cell)]
            if len(find_cell) > 0:
                idx = find_cell[0]
                cell = triangles_of_node[next_pts[2]][idx]
                idx = np.where((cells == cell).all(axis=1))[0][0]
                tri2tri_pairs.append((i, idx))
        elif len(set(next_pts)) == 2:
            find_edge = [idx for idx, edge in enumerate(edges_of_node[next_pts[0]]) if set(next_pts) == set(edge)]
            if len(find_edge) > 0:
                tri2edge_idx.append(i)
            else:
                tri2path_idx.append(i)
        elif len(set(next_pts)) == 1:
            tri2node_idx.append(i)
    with open("tri2tri_pairs.pickle", "wb") as f:
        pickle.dump(tri2tri_pairs, f)
    with open("tri2edge_idx.pickle", "wb") as f:
        pickle.dump(tri2edge_idx, f)
    with open("tri2path_idx.pickle", "wb") as f:
        pickle.dump(tri2path_idx, f)
    with open("tri2node_idx.pickle", "wb") as f:
        pickle.dump(tri2node_idx, f)
