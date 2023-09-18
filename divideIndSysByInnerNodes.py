'''
This code splits independent systems whose upper meshes have inner nodes..

Author: Congrong Ren
Date: Mar 05, 2023
'''
import pickle
import numpy as np
import networkx as nx


def cellFormatConvert(ind_sys, cells):
    '''
    Args:
        ind_sys: [[all lower cells' indices],[all upper cells' indices]]
        cells: [[three nodes for all cells]]
    Returns:
        adjacency lists for upper and lower meshes
        dual graphs for upper and lower meshes
        inner nodes in upper meshes
    '''
    lower_adjacency_list = {}
    upper_adjacency_list = {}
    lower_edge_to_triangles = {}
    upper_edge_to_triangles = {}
    lower_dual_graph = nx.Graph()
    upper_dual_graph = nx.Graph()
    upper_nodes = []
    upper_boundary_nodes = []
    for cell_idx in ind_sys[0]:
        cell = list(cells[cell_idx])
        cell.sort()
        try:
            lower_adjacency_list[cell[0]].append(cell[1])
            lower_adjacency_list[cell[0]].append(cell[2])
        except:
            lower_adjacency_list[cell[0]] = []
            lower_adjacency_list[cell[0]].append(cell[1])
            lower_adjacency_list[cell[0]].append(cell[2])
        try:
            lower_adjacency_list[cell[1]].append(cell[0])
            lower_adjacency_list[cell[1]].append(cell[2])
        except:
            lower_adjacency_list[cell[1]] = []
            lower_adjacency_list[cell[1]].append(cell[0])
            lower_adjacency_list[cell[1]].append(cell[2])
        try:
            lower_adjacency_list[cell[2]].append(cell[0])
            lower_adjacency_list[cell[2]].append(cell[1])
        except:
            lower_adjacency_list[cell[2]] = []
            lower_adjacency_list[cell[2]].append(cell[0])
            lower_adjacency_list[cell[2]].append(cell[1])
        try:
            lower_edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        except:
            lower_edge_to_triangles[(cell[0], cell[1])] = []
            lower_edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        try:
            lower_edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        except:
            lower_edge_to_triangles[(cell[0], cell[2])] = []
            lower_edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        try:
            lower_edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
        except:
            lower_edge_to_triangles[(cell[1], cell[2])] = []
            lower_edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
    for node, neighbors in lower_adjacency_list.items():
        lower_adjacency_list[node] = list(set(neighbors))
    for edge, triangles in lower_edge_to_triangles.items():
        if len(triangles) == 2:
            lower_dual_graph.add_edge(triangles[0], triangles[1])
    for cell_idx in ind_sys[1]:
        cell = list(cells[cell_idx])
        cell.sort()
        upper_nodes += cell
        try:
            upper_adjacency_list[cell[0]].append(cell[1])
            upper_adjacency_list[cell[0]].append(cell[2])
        except:
            upper_adjacency_list[cell[0]] = []
            upper_adjacency_list[cell[0]].append(cell[1])
            upper_adjacency_list[cell[0]].append(cell[2])
        try:
            upper_adjacency_list[cell[1]].append(cell[0])
            upper_adjacency_list[cell[1]].append(cell[2])
        except:
            upper_adjacency_list[cell[1]] = []
            upper_adjacency_list[cell[1]].append(cell[0])
            upper_adjacency_list[cell[1]].append(cell[2])
        try:
            upper_adjacency_list[cell[2]].append(cell[0])
            upper_adjacency_list[cell[2]].append(cell[1])
        except:
            upper_adjacency_list[cell[2]] = []
            upper_adjacency_list[cell[2]].append(cell[0])
            upper_adjacency_list[cell[2]].append(cell[1])
        try:
            upper_edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        except:
            upper_edge_to_triangles[(cell[0], cell[1])] = []
            upper_edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        try:
            upper_edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        except:
            upper_edge_to_triangles[(cell[0], cell[2])] = []
            upper_edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        try:
            upper_edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
        except:
            upper_edge_to_triangles[(cell[1], cell[2])] = []
            upper_edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
    for node, neighbors in upper_adjacency_list.items():
        upper_adjacency_list[node] = list(set(neighbors))
    for edge, triangles in upper_edge_to_triangles.items():
        if len(triangles) == 2:
            upper_dual_graph.add_edge(triangles[0], triangles[1])
        elif len(triangles) == 1:
            upper_boundary_nodes += list(edge)
    upper_inner_nodes = list(set(upper_nodes) - set(upper_boundary_nodes))
    return lower_adjacency_list, upper_adjacency_list, lower_edge_to_triangles, upper_edge_to_triangles, lower_dual_graph, upper_dual_graph, upper_inner_nodes


def getMeshBoundary(cell_indices, cells):
    edge_to_triangles = {}
    for cell_idx in cell_indices:
        cell = list(cells[cell_idx])
        cell.sort()
        try:
            edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        except:
            edge_to_triangles[(cell[0], cell[1])] = []
            edge_to_triangles[(cell[0], cell[1])].append(cell_idx)
        try:
            edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        except:
            edge_to_triangles[(cell[0], cell[2])] = []
            edge_to_triangles[(cell[0], cell[2])].append(cell_idx)
        try:
            edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
        except:
            edge_to_triangles[(cell[1], cell[2])] = []
            edge_to_triangles[(cell[1], cell[2])].append(cell_idx)
    boundary_nodes = []
    for edge, triangles in edge_to_triangles.items():
        if len(triangles) == 1:
            boundary_nodes += list(edge)
    return list(set(boundary_nodes))


def findCuttingEdges(lower_adj_list, upper_adj_list, lower_edge_to_triangles, inner_node, previous_nodes, node_coors):
    neighbors = upper_adj_list[inner_node]
    cutting_pairs_lower = []
    cutting_pairs_upper = []
    cutting_pairs_angle_cos = []
    for i, n1 in enumerate(neighbors):
        for n2 in neighbors[i + 1:]:
            if n1 in previous_nodes.keys() and n2 in previous_nodes.keys():
                prev_n1 = list(set(previous_nodes[n1]).intersection(set(lower_adj_list.keys())))[0]
                prev_n2 = list(set(previous_nodes[n2]).intersection(set(lower_adj_list.keys())))[0]
                if prev_n1 in lower_adj_list[prev_n2]:
                    if prev_n1 < prev_n2:
                        if len(lower_edge_to_triangles[(prev_n1, prev_n2)]) == 2:
                            cutting_pairs_lower.append([prev_n1, prev_n2])
                            cutting_pairs_upper.append([n1, inner_node, n2])
                            cutting_pairs_angle_cos.append(
                                cosOfAngle(node_coors[inner_node], node_coors[n1], node_coors[n2]))
                    else:
                        if len(lower_edge_to_triangles[(prev_n2, prev_n1)]) == 2:
                            cutting_pairs_lower.append([prev_n2, prev_n1])
                            cutting_pairs_upper.append([n2, inner_node, n1])
                            cutting_pairs_angle_cos.append(
                                cosOfAngle(node_coors[inner_node], node_coors[n1], node_coors[n2]))
    if len(cutting_pairs_angle_cos) > 0:
        idx = cutting_pairs_angle_cos.index(min(cutting_pairs_angle_cos))
        return cutting_pairs_lower[idx], cutting_pairs_upper[idx]
    else:
        return None, None


def cosOfAngle(o, p1, p2):
    v1 = np.array(p1 - o)
    v2 = np.array(p2 - o)
    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    num_cells = cells.shape[0]
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    previous_nodes = {}
    for i, node in enumerate(next_nodes):
        try:
            previous_nodes[node].append(i)
        except:
            previous_nodes[node] = []
            previous_nodes[node].append(i)

    with open("independent_patches_cellID.pickle", "rb") as f:
        independent_patches_cellID = pickle.load(f)
    independent_subpatches_cellID = []
    # for ind_sys_idx, indSys in enumerate(independent_patches_cellID):
    for ind_sys_idx, indSys in enumerate(independent_patches_cellID):
        lower_adj_list, upper_adj_list, lower_edge_to_triangles, upper_edge_to_triangles, lower_dual, upper_dual, upper_inner_nodes = cellFormatConvert(
            indSys, cells)
        if len(upper_inner_nodes) > 0 and ind_sys_idx != 37:
            for inner_n in upper_inner_nodes:
                lower_edge, upper_edges = findCuttingEdges(lower_adj_list, upper_adj_list, lower_edge_to_triangles,
                                                           inner_n, previous_nodes, nodes_coor)
                if lower_edge is None:
                    continue
                lower_triangles = lower_edge_to_triangles[(lower_edge[0], lower_edge[1])]
                upper_edge1 = upper_edges[:-1]
                upper_edge1.sort()
                upper_edge2 = upper_edges[1:]
                upper_edge2.sort()
                upper_triangles1 = upper_edge_to_triangles[(upper_edge1[0], upper_edge1[1])]
                upper_triangles2 = upper_edge_to_triangles[(upper_edge2[0], upper_edge2[1])]
                if (lower_triangles[0], lower_triangles[1]) in lower_dual.edges and (
                        upper_triangles1[0], upper_triangles1[1]) in upper_dual.edges and (
                        upper_triangles2[0], upper_triangles2[1]) in upper_dual.edges:
                    lower_dual.remove_edge(lower_triangles[0], lower_triangles[1])
                    upper_dual.remove_edge(upper_triangles1[0], upper_triangles1[1])
                    upper_dual.remove_edge(upper_triangles2[0], upper_triangles2[1])
            comps_lower = nx.connected_components(lower_dual)
            comps_lower = list(comps_lower)
            comps_upper = nx.connected_components(upper_dual)
            comps_upper = list(comps_upper)
            comps_boundary_upper = []
            for comp in comps_upper:
                boundary_nodes = getMeshBoundary(comp, cells)
                boundary_nodes = set(boundary_nodes) - set(upper_inner_nodes)
                comps_boundary_upper.append(boundary_nodes)
            for lower_idx, comp in enumerate(comps_lower):
                boundary_nodes = getMeshBoundary(comp, cells)
                boundary_nodes_upper_matching = set([next_nodes[node] for node in boundary_nodes])
                upper_idx = [idx for idx, comps_boundary in enumerate(comps_boundary_upper) if
                             boundary_nodes_upper_matching.intersection(
                                 comps_boundary) == boundary_nodes_upper_matching][0]
                independent_subpatches_cellID.append([comps_lower[lower_idx], comps_upper[upper_idx]])
        else:
            independent_subpatches_cellID.append(indSys)

    with open("subpatches_cellID.pickle", "wb") as f:
        pickle.dump(independent_subpatches_cellID, f)
