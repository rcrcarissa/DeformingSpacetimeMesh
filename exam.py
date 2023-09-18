import numpy as np
import matplotlib.pyplot as plt
import pickle


def is_four_points_coplanar(p, q, r, s):
    """
    Args: Coordinates of four points
    Return: True/False
    """
    A = q - p
    B = r - p
    C = s - p
    return np.abs(np.dot(A, np.cross(B, C))) < 1e-7


def node2coor(node, nodes_coor):
    node_coor_3d = np.zeros(3)
    node_coor_3d[:-1] = nodes_coor[node[0]]
    node_coor_3d[-1] = int(node[1] == 'u')
    return node_coor_3d


def exportAdjacencyList(cells, nodes_of_cells):
    node_list = []
    edge_list = []
    for cell in cells:
        cell = nodes_of_cells[cell]
        node_list += [(cell[0], 'u'), (cell[1], 'u'), (cell[2], 'u')]
        sorted_nodes = np.sort(cell)
        edge_list.append(((sorted_nodes[0], 'u'), (sorted_nodes[1], 'u')))
        edge_list.append(((sorted_nodes[0], 'u'), (sorted_nodes[2], 'u')))
        edge_list.append(((sorted_nodes[1], 'u'), (sorted_nodes[2], 'u')))
    node_list = list(set(node_list))
    edge_list = list(set(edge_list))
    adjacency_list = {}
    for node in node_list:
        adjacency_list[node] = []
    for edge in edge_list:
        adjacency_list[edge[0]].append(edge[1])
        adjacency_list[edge[1]].append(edge[0])
    return edge_list, adjacency_list


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)
    edges_to_draw = [((695, 'l'), (778, 'l')), ((959, 'l'), (960, 'l')), ((867, 'l'), (960, 'l')),
                     ((866, 'l'), (867, 'l')), ((779, 'l'), (867, 'l')), ((866, 'l'), (959, 'l')),
                     ((865, 'l'), (959, 'l')), ((865, 'l'), (866, 'l')), ((695, 'l'), (779, 'l')),
                     ((778, 'l'), (779, 'l')), ((694, 'l'), (778, 'l')), ((867, 'l'), (959, 'l')),
                     ((694, 'l'), (695, 'l')), ((778, 'l'), (866, 'l')), ((778, 'l'), (867, 'l')),
                     ((783, 'u'), (871, 'u')), ((783, 'u'), (872, 'u')), ((871, 'u'), (965, 'u')),
                     ((872, 'u'), (873, 'u')), ((873, 'u'), (966, 'u')), ((871, 'u'), (872, 'u')),
                     ((699, 'u'), (700, 'u')), ((872, 'u'), (966, 'u')), ((784, 'u'), (872, 'u')),
                     ((872, 'u'), (965, 'u')), ((699, 'u'), (783, 'u')), ((699, 'u'), (784, 'u')),
                     ((783, 'u'), (784, 'u')), ((965, 'u'), (966, 'u')), ((784, 'u'), (873, 'u')),
                     ((700, 'u'), (784, 'u')), ((694, 'l'), (699, 'u')), ((695, 'l'), (700, 'u')),
                     ((778, 'l'), (783, 'u')), ((779, 'l'), (873, 'u')), ((865, 'l'), (871, 'u')),
                     ((866, 'l'), (783, 'u')), ((867, 'l'), (873, 'u')), ((959, 'l'), (965, 'u')),
                     ((960, 'l'), (966, 'u')), ((694, 'l'), (700, 'u')), ((695, 'l'), (784, 'u')),
                     ((779, 'l'), (784, 'u')), ((867, 'l'), (966, 'u')), ((959, 'l'), (966, 'u')),
                     ((865, 'l'), (965, 'u')), ((866, 'l'), (871, 'u')), ((694, 'l'), (783, 'u'))]

    # cells_to_draw = [504, 505, 29, 142, 520, 521, 28, 140]
    # e, a = exportAdjacencyList(cells_to_draw, cells)
    # print(e)
    ax = plt.figure().add_subplot(projection='3d')
    for e in edges_to_draw:
        ax.plot([nodes_coor[e[0][0]][0], nodes_coor[e[1][0]][0]], [nodes_coor[e[0][0]][1], nodes_coor[e[1][0]][1]],
                [int(e[0][1] == 'u'), int(e[1][1] == 'u')], color='gray', linewidth=0.8)
    for e in []:
        ax.plot([nodes_coor[e[0][0]][0], nodes_coor[e[1][0]][0]], [nodes_coor[e[0][0]][1], nodes_coor[e[1][0]][1]],
                [int(e[0][1] == 'u'), int(e[1][1] == 'u')], color='blue', linewidth=0.8)
    # ax.plot(nodes_coor[54672][0], nodes_coor[54672][1], 0, 'ro')
    # ax.plot(nodes_coor[next_nodes[54672]][0], nodes_coor[next_nodes[54672]][1], 1, 'ro')
    plt.show()
    # print(is_four_points_coplanar(node2coor((434, 'l'), nodes_coor), node2coor((500, 'l'), nodes_coor),
    #                               node2coor((437, 'u'), nodes_coor), node2coor((500, 'u'), nodes_coor)))
    # normal = np.cross(node2coor((782, 'u'), nodes_coor) - node2coor((871, 'u'), nodes_coor),
    #                   node2coor((777, 'l'), nodes_coor) - node2coor((871, 'u'), nodes_coor))
    # print(np.dot(normal, node2coor((776, 'l'), nodes_coor) - node2coor((871, 'u'), nodes_coor)))
