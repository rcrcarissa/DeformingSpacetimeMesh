import numpy as np
import pickle
import copy

if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)
    with open("tri2edge_idx.pickle", "rb") as f:
        tri2edge_idx = pickle.load(f)
    with open("tri2node_idx.pickle", "rb") as f:
        tri2node_idx = pickle.load(f)

    all_tetrahedra = {}
    for idx in tri2node_idx:
        all_tetrahedra[idx] = [
            [(cells[idx][0], 'l'), (cells[idx][1], 'l'), (cells[idx][2], 'l'), (next_nodes[cells[idx][0]], 'u')]]
    for idx in tri2edge_idx:
        cell = copy.deepcopy(cells[idx])
        cell = np.sort(cell)
        if next_nodes[cell[0]] == next_nodes[cell[1]]:
            all_tetrahedra[idx] = [
                [(cell[0], 'l'), (cell[1], 'l'), (cell[2], 'l'), (next_nodes[cell[2]], 'u')],
                [(cell[0], 'l'), (cell[1], 'l'), (next_nodes[cell[0]], 'u'), (next_nodes[cell[2]], 'u')]]
        elif next_nodes[cell[0]] == next_nodes[cell[2]]:
            all_tetrahedra[idx] = [
                [(cell[0], 'l'), (cell[1], 'l'), (cell[2], 'l'), (next_nodes[cell[0]], 'u')],
                [(cell[0], 'l'), (cell[1], 'l'), (next_nodes[cell[0]], 'u'), (next_nodes[cell[1]], 'u')]]
        else:
            all_tetrahedra[idx] = [
                [(cell[0], 'l'), (cell[1], 'l'), (cell[2], 'l'), (next_nodes[cell[1]], 'u')]]
    with open("mergedTri_tetrahedra-1.pickle", "wb") as f:
        pickle.dump(all_tetrahedra, f)
