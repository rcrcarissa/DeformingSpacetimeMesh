import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

cells_to_draw = [209, 276, 389, 536, 567, 573, 817, 1022, 1816]

if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    mesh = Triangulation(nodes_coor[:, 0], nodes_coor[:, 1], cells)
    fig, ax = plt.subplots()
    ax.triplot(mesh, color='grey', linewidth=0.5)
    edges_to_draw = []
    for c in cells_to_draw:
        nodes = np.sort(cells[c])
        if (nodes[0], nodes[1]) not in edges_to_draw:
            edges_to_draw.append((nodes[0], nodes[1]))
        if (nodes[0], nodes[2]) not in edges_to_draw:
            edges_to_draw.append((nodes[0], nodes[2]))
        if (nodes[1], nodes[2]) not in edges_to_draw:
            edges_to_draw.append((nodes[1], nodes[2]))
    for e in edges_to_draw:
        ax.plot([nodes_coor[e[0]][0], nodes_coor[e[1]][0]], [nodes_coor[e[0]][1], nodes_coor[e[1]][1]], color="b",
                linewidth=1)
    plt.show()
