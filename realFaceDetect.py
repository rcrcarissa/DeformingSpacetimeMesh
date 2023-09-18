'''
This code detects all real faces and dumps them into a file.
Author: Congrong Ren
Date: Nov 28, 2022
'''

import numpy as np
import pickle

if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    num_cells = cells.shape[0]
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    # adj_lists: {node:[all nodes adjacent to it]}
    adj_lists = {}
    edges_list = []
    for i in range(num_nodes):
        adj_lists[i] = []  # initialization
    for i, cl in enumerate(cells):
        cl = np.sort(cl)
        adj_lists[cl[0]] += [cl[1], cl[2]]
        adj_lists[cl[1]] += [cl[0], cl[2]]
        adj_lists[cl[2]] += [cl[0], cl[1]]
        if [cl[0], cl[1]] not in edges_list:
            edges_list.append([cl[0], cl[1]])
        if [cl[0], cl[2]] not in edges_list:
            edges_list.append([cl[0], cl[2]])
        if [cl[1], cl[2]] not in edges_list:
            edges_list.append([cl[1], cl[2]])

    # find all real faces: a list of (i, j, j', i'). i < j
    #      and half real faces: a list of (i, j, i'). i < j
    real_faces = []
    half_real_faces = []
    for i, edge in enumerate(edges_list):
        node1 = edge[0]
        node2 = edge[1]
        next_node1 = next_nodes[node1]
        next_node2 = next_nodes[node2]
        if next_node2 in adj_lists[next_node1]:
            real_faces.append((node1, node2, next_node2, next_node1))
        elif next_node1 == next_node2:
            half_real_faces.append((node1, node2, next_node1))
        if i > 0 and i % 10000 == 0:
            print(str(i) + "/" + str(len(edges_list)) + " edges have been processed!")
    with open("real_faces.pickle", "wb") as f:
        pickle.dump(real_faces, f)
    with open("half_real_faces.pickle", "wb") as f:
        pickle.dump(half_real_faces, f)
