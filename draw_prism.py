'''
This code draws one independent system and its correlated real face(s).
Author: Congrong Ren
Last Modified: Dec 29, 2022
'''
import pickle
import numpy as np
import matplotlib.pyplot as plt

# indSys_indices = range(1375)
indSys_idx = 51

upper_height = 0.5

if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    with open("real_faces.pickle", "rb") as f:
        real_faces = pickle.load(f)
    with open("half_real_faces.pickle", "rb") as f:
        half_real_faces = pickle.load(f)
    with open("comp_lower.pickle", "rb") as f:
        comp_lower = pickle.load(f)
    with open("comp_upper.pickle", "rb") as f:
        comp_upper = pickle.load(f)
    ind_patch = [[54], [527]]
    edges_lower = set()
    edges_upper = set()
    for c_i in ind_patch[0]:
        c = list(cells[c_i])
        c.sort()
        edges_lower.add((c[0], c[1]))
        edges_lower.add((c[0], c[2]))
        edges_lower.add((c[1], c[2]))
    for c_i in ind_patch[1]:
        c = list(cells[c_i])
        c.sort()
        edges_upper.add((c[0], c[1]))
        edges_upper.add((c[0], c[2]))
        edges_upper.add((c[1], c[2]))
    related_realFaces = []
    for e in edges_lower:
        if (e[0], e[1], next_nodes[e[1]], next_nodes[e[0]]) in real_faces or (
                e[0], e[1], next_nodes[e[0]]) in half_real_faces:
            related_realFaces.append((e[0], e[1], next_nodes[e[1]], next_nodes[e[0]]))
    ax = plt.figure().add_subplot(projection='3d')
    for e in edges_lower:
        ax.plot(nodes_coor[list(e)][:, 0], nodes_coor[list(e)][:, 1], [0, 0], color='gray', linewidth=0.8)
    for e in edges_upper:
        ax.plot(nodes_coor[list(e)][:, 0], nodes_coor[list(e)][:, 1], [upper_height, upper_height], color='gray',
                linewidth=0.8)
    for f in related_realFaces:
        x = np.array([nodes_coor[[f[0], f[3]]][:, 0], nodes_coor[[f[1], f[2]]][:, 0]])
        y = np.array([nodes_coor[[f[0], f[3]]][:, 1], nodes_coor[[f[1], f[2]]][:, 1]])
        z = np.array([[0, upper_height], [0, upper_height]])
        ax.plot(x[0], y[0], z[0], color='b', linewidth=0.8)
        ax.plot(x[1], y[1], z[1], color='b', linewidth=0.8)
        ax.plot_surface(x, y, z, color='b', alpha=0.3)
        if f[0] < f[1]:
            ax.plot(nodes_coor[[f[0], f[2]]][:, 0], nodes_coor[[f[0], f[2]]][:, 1], [0, upper_height], color='b',
                    linewidth=0.5)
        else:
            ax.plot(nodes_coor[[f[1], f[3]]][:, 0], nodes_coor[[f[1], f[3]]][:, 1], [0, upper_height], color='b',
                    linewidth=0.5)
    plt.show()
    # for indSys_idx in indSys_indices:
    #     ind_patch = patches_cellID[indSys_idx]
    #     if len(ind_patch[0]) < 3 or len(ind_patch[1]) < 3:
    #         continue
    #     edges_lower = set()
    #     edges_upper = set()
    #     for c_i in ind_patch[0]:
    #         c = list(cells[c_i])
    #         c.sort()
    #         edges_lower.add((c[0], c[1]))
    #         edges_lower.add((c[0], c[2]))
    #         edges_lower.add((c[1], c[2]))
    #     for c_i in ind_patch[1]:
    #         c = list(cells[c_i])
    #         c.sort()
    #         edges_upper.add((c[0], c[1]))
    #         edges_upper.add((c[0], c[2]))
    #         edges_upper.add((c[1], c[2]))
    #     related_realFaces = []
    #     for e in edges_lower:
    #         if (e[0], e[1], next_nodes[e[1]], next_nodes[e[0]]) in real_faces or (
    #                 e[0], e[1], next_nodes[e[0]]) in half_real_faces:
    #             related_realFaces.append((e[0], e[1], next_nodes[e[1]], next_nodes[e[0]]))
    #     ax = plt.figure().add_subplot(projection='3d')
    #     for e in edges_lower:
    #         ax.plot(nodes_coor[list(e)][:, 0], nodes_coor[list(e)][:, 1], [0, 0], color='gray', linewidth=0.8)
    #     for e in edges_upper:
    #         ax.plot(nodes_coor[list(e)][:, 0], nodes_coor[list(e)][:, 1], [upper_height, upper_height], color='gray',
    #                 linewidth=0.8)
    #     for f in related_realFaces:
    #         x = np.array([nodes_coor[[f[0], f[3]]][:, 0], nodes_coor[[f[1], f[2]]][:, 0]])
    #         y = np.array([nodes_coor[[f[0], f[3]]][:, 1], nodes_coor[[f[1], f[2]]][:, 1]])
    #         z = np.array([[0, upper_height], [0, upper_height]])
    #         ax.plot(x[0], y[0], z[0], color='b', linewidth=0.8)
    #         ax.plot(x[1], y[1], z[1], color='b', linewidth=0.8)
    #         ax.plot_surface(x, y, z, color='b', alpha=0.3)
    #         if f[0] < f[1]:
    #             ax.plot(nodes_coor[[f[0], f[2]]][:, 0], nodes_coor[[f[0], f[2]]][:, 1], [0, upper_height], color='b',
    #                     linewidth=0.5)
    #         else:
    #             ax.plot(nodes_coor[[f[1], f[3]]][:, 0], nodes_coor[[f[1], f[3]]][:, 1], [0, upper_height], color='b',
    #                     linewidth=0.5)
    #     plt.savefig('indSys/' + str(indSys_idx) + ".png")
