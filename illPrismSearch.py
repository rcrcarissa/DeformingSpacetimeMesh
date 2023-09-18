'''
This code finds all tri2tri cases with ill twist.

Author: Congrong Ren
Date: Mar 23, 2023
'''
import numpy as np
import pickle


def intersect_line_plane(pt0, pt1, pl0, pl1, pl2):
    '''
    Arguments:
    pt0, pt1: Define the line.
    pl0, pl1, pl2: Define the plane.
    Return:
    A point where the line intersects the plane.
    '''
    p_norm = np.cross(pl1 - pl0, pl2 - pl0)  # normal vector to the plane
    u = pt1 - pt0
    dot = np.dot(p_norm, u)
    if dot == 0:
        return None
    w = pl0 - pt0
    fra = np.dot(p_norm, w) / dot
    u = u * fra
    return pt0 + u


def is_point_in_triangle(pt, v0, v1, v2):
    '''
    Arguments:
    pt: Define the point. It is in the same plane of the triangle.
    v0, v1, v2: Define the triangle.
    Return:
    A boolean variable indicating if pt is in the triangle.
    '''
    area = np.linalg.norm(np.cross(v1 - v0, v2 - v0))
    alpha = np.linalg.norm(np.cross(v0 - pt, v1 - pt)) / area
    beta = np.linalg.norm(np.cross(v0 - pt, v2 - pt)) / area
    gamma = np.linalg.norm(np.cross(v1 - pt, v2 - pt)) / area
    if np.abs(alpha + beta + gamma - 1) < 1e-7:
        return True
    else:
        return False


def is_point_on_segment(pt, v0, v1):
    '''
    Args:
        pt: Define the point.
        v0, v1: Define the segment.
    Returns:
        True or False.
    '''
    cross_prod = np.cross(pt - v0, v1 - v0)
    if np.linalg.norm(cross_prod) < 1e-7:
        ratio = (pt[0] - v0[0]) / (v1[0] - v0[0])
        if 0 <= ratio <= 1:
            return True
    return False


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    with open("dataGen/next_nodes.pickle", "rb") as f:
        next_nodes = pickle.load(f)

    with open("tri2tri_pairs.pickle", "rb") as f:
        tri2tri_pairs = pickle.load(f)

    ill_twist_sys = []
    for i, indSys in enumerate(tri2tri_pairs):
        sorted_lower_nodes = np.sort(cells[indSys[0]])
        n1 = np.append(nodes_coor[sorted_lower_nodes[0]], 0)
        n2 = np.append(nodes_coor[sorted_lower_nodes[1]], 0)
        n1_prime = np.append(nodes_coor[next_nodes[sorted_lower_nodes[0]]], 1)
        n2_prime = np.append(nodes_coor[next_nodes[sorted_lower_nodes[1]]], 1)
        n3_prime = np.append(nodes_coor[next_nodes[sorted_lower_nodes[2]]], 1)
        intersect = intersect_line_plane(n2, n3_prime, n1, n1_prime, n2_prime)
        if intersect is not None and is_point_in_triangle(intersect, n1, n1_prime, n2_prime):
            if is_point_on_segment(intersect, n2, n3_prime):
                ill_twist_sys.append(i)
    print(len(ill_twist_sys))
    with open("ill_twist_sys.pickle", "wb") as f:
        pickle.dump(ill_twist_sys, f)
