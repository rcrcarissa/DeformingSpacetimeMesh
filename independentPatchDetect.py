'''
This code detects all independent systems devided by real faces and dumps them into a file.
Finding independent patches is equivalent to finding connected components on a graph. The algorithm is based on graph traversal.
Author: Congrong Ren
Date: Dec 11, 2022
'''

import numpy as np
import networkx as nx
import pickle


def boundary_of_triangle_mesh(cells_idx, nodes_of_cells):
    '''
    Args:
        cells_idx:
        Indices of the cells composing the mesh.
        nodes_of_cells:
        Indices of nodes on the cells.
    Returns:
        A list of unordered edges representing the boundary of the mesh.
    '''
    # find all edges that appear only once
    single_edges = []
    for c in cells_idx:
        nodes_of_c = np.sort(nodes_of_cells[c])
        if (nodes_of_c[0], nodes_of_c[1]) in single_edges:
            single_edges.remove((nodes_of_c[0], nodes_of_c[1]))
        else:
            single_edges.append((nodes_of_c[0], nodes_of_c[1]))
        if (nodes_of_c[0], nodes_of_c[2]) in single_edges:
            single_edges.remove((nodes_of_c[0], nodes_of_c[2]))
        else:
            single_edges.append((nodes_of_c[0], nodes_of_c[2]))
        if (nodes_of_c[1], nodes_of_c[2]) in single_edges:
            single_edges.remove((nodes_of_c[1], nodes_of_c[2]))
        else:
            single_edges.append((nodes_of_c[1], nodes_of_c[2]))
    return single_edges


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    num_cells = cells.shape[0]
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    lower_boundary = []
    upper_boundary = []
    with open("real_faces.pickle", "rb") as f:
        real_faces = pickle.load(f)
    with open("half_real_faces.pickle", "rb") as f:
        half_real_faces = pickle.load(f)
    for face in real_faces:
        lower_boundary.append((face[0], face[1]))
        if face[3] < face[2]:
            upper_boundary.append((face[3], face[2]))
        else:
            upper_boundary.append((face[2], face[3]))
    for face in half_real_faces:
        lower_boundary.append((face[0], face[1]))

    # dictionary: {edge: 1 (boundary) or 2 adjacent triangles}
    edge2tri_dict = {}
    for i, cl in enumerate(cells):
        c = np.sort(cl)
        try:
            edge2tri_dict[(c[0], c[1])].append(i)
        except:
            edge2tri_dict[(c[0], c[1])] = [i]
        try:
            edge2tri_dict[(c[0], c[2])].append(i)
        except:
            edge2tri_dict[(c[0], c[2])] = [i]
        try:
            edge2tri_dict[(c[1], c[2])].append(i)
        except:
            edge2tri_dict[(c[1], c[2])] = [i]

    # construct dual graphs for both upper layer and lower layer
    g_lower = nx.Graph()
    g_lower.add_nodes_from(range(num_cells))
    g_upper = nx.Graph()
    g_upper.add_nodes_from(range(num_cells))
    i = 0
    for e, cs in edge2tri_dict.items():
        if len(cs) == 2:
            if e not in lower_boundary:
                g_lower.add_edge(cs[0], cs[1])
            if e not in upper_boundary:
                g_upper.add_edge(cs[0], cs[1])
        i += 1

    # get all connected components
    comp_lower = nx.connected_components(g_lower)
    comp_lower = list(comp_lower)
    comp_upper = nx.connected_components(g_upper)
    comp_upper = list(comp_upper)

    with open("comp_lower.pickle", "wb") as f:
        pickle.dump(comp_lower, f)
    with open("comp_upper.pickle", "wb") as f:
        pickle.dump(comp_upper, f)

    # with open("comp_lower.pickle", "rb") as f:
    #     comp_lower = pickle.load(f)
    # with open("comp_upper.pickle", "rb") as f:
    #     comp_upper = pickle.load(f)

    # match the patches in upper and lower layers
    boundary_nodes_lower = []
    boundary_nodes_upper = []
    for i, comp in enumerate(comp_lower):
        comp_boundary = boundary_of_triangle_mesh(comp, cells)
        comp_bound_in_realFace = []
        for edge in comp_boundary:
            if edge in lower_boundary:
                comp_bound_in_realFace.append(edge)
        boundary_nodes_next_layer = [next_nodes[edge[0]] for edge in comp_bound_in_realFace]
        boundary_nodes_next_layer += [next_nodes[edge[1]] for edge in comp_bound_in_realFace]
        boundary_nodes_lower.append(set(boundary_nodes_next_layer))
    for i, comp in enumerate(comp_upper):
        comp_boundary = boundary_of_triangle_mesh(comp, cells)
        comp_bound_in_realFace = []
        for edge in comp_boundary:
            if edge in upper_boundary:
                comp_bound_in_realFace.append(edge)
        boundary_nodes = [edge[0] for edge in comp_bound_in_realFace]
        boundary_nodes += [edge[1] for edge in comp_bound_in_realFace]
        boundary_nodes_upper.append(set(boundary_nodes))
    with open("boundary_nodes_lower.pickle", "wb") as f:
        pickle.dump(boundary_nodes_lower, f)
    with open("boundary_nodes_upper.pickle", "wb") as f:
        pickle.dump(boundary_nodes_upper, f)

    # with open("boundary_nodes_lower.pickle", "rb") as f:
    #     boundary_nodes_lower = pickle.load(f)
    # with open("boundary_nodes_upper.pickle", "rb") as f:
    #     boundary_nodes_upper = pickle.load(f)

    independent_sys = []
    unpaired_lower = []
    unpaired_upper = list(range(len(comp_upper)))
    for i, nodes in enumerate(boundary_nodes_lower):
        # if nodes in boundary_nodes_upper:
        try:
            idx = boundary_nodes_upper.index(nodes)
            independent_sys.append([comp_lower[i], comp_upper[idx]])
            unpaired_upper.remove(idx)
        # else:
        except:
            unpaired_lower.append(i)
    with open("independent_patches_cellID.pickle", "wb") as f:
        pickle.dump(independent_sys, f)
    with open("unpaired_lower.pickle", "wb") as f:
        pickle.dump(unpaired_lower, f)
    with open("unpaired_upper.pickle", "wb") as f:
        pickle.dump(unpaired_upper, f)
