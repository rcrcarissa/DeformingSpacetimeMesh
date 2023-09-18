'''
This code draws one independent system and its correlated real face(s).
Author: Congrong Ren
Last Modified: Dec 29, 2022
'''
import pickle
import numpy as np
import vtk
import matplotlib.pyplot as plt

nphi = 16
upper_height = 0.1
num_upsampled = 1024
signal_mag = 5
R = 1
signal_R = 0.3 * R
signal_cr = R / 2


def fwhm2sigma(fwhm):
    return fwhm / np.sqrt(8 * np.log(2))


def BivariateGaussian(point, center, sigma_x, sigma_y, rho):
    return np.exp(-(
            ((point[0] - center[0]) / sigma_x) ** 2 - 2 * rho * (point[0] - center[0]) * (point[1] - center[1]) / (
            sigma_x * sigma_y) + ((point[1] - center[1]) / sigma_y) ** 2) / (2 * (1 - rho ** 2))) / (
            2 * np.pi * sigma_x * sigma_y * np.sqrt(1 - rho ** 2))


def tet_volume(p, q, r, s):
    v0 = q - p
    v1 = r - p
    v2 = s - p
    return np.abs(np.dot(np.cross(v0, v1), v2)) / 6


def tet_barycentric_interp(point_coor, tet_coors, values):
    total_volume = tet_volume(tet_coors[0], tet_coors[1], tet_coors[2], tet_coors[3])
    a0 = tet_volume(point_coor, tet_coors[1], tet_coors[2], tet_coors[3]) / total_volume
    a1 = tet_volume(point_coor, tet_coors[0], tet_coors[2], tet_coors[3]) / total_volume
    a2 = tet_volume(point_coor, tet_coors[0], tet_coors[1], tet_coors[3]) / total_volume
    a3 = tet_volume(point_coor, tet_coors[0], tet_coors[1], tet_coors[2]) / total_volume
    return a0 * values[0] + a1 * values[1] + a2 * values[2] + a3 * values[3]


def tet_cellLocator(point, cellLocator):
    closest_pt = [0.0, 0.0, 0.0]
    cellId = vtk.reference(0)
    subId = vtk.reference(0)
    d = vtk.reference(np.float64())
    cellLocator.FindClosestPoint(point, closest_pt, cellId, subId, d)
    if d < 1e-7:
        return int(cellId)
    return closest_pt


if __name__ == "__main__":
    with open("dataGen/syn_mesh.pickle", "rb") as f:
        nodes_coor, cells = pickle.load(f)
    num_nodes = nodes_coor.shape[0]
    with open("dataGen/values.pickle", "rb") as f:
        values = pickle.load(f)
    with open("dataGen/values_denoise.pickle", "rb") as f:
        values_denoise = pickle.load(f)

    straight_vtu_fname = "straightMesh.vtu"
    deformed_vtu_fname = "all.vtu"

    reader0 = vtk.vtkXMLUnstructuredGridReader()
    reader0.SetFileName(straight_vtu_fname)
    reader0.Update()
    data0 = reader0.GetOutput()
    points_coor0 = np.array(data0.GetPoints().GetData())
    tets0 = np.array(data0.GetCells().GetData())
    cell_locator0 = vtk.vtkCellLocator()
    cell_locator0.SetDataSet(data0)
    cell_locator0.BuildLocator()

    reader1 = vtk.vtkXMLUnstructuredGridReader()
    reader1.SetFileName(deformed_vtu_fname)
    reader1.Update()
    data1 = reader1.GetOutput()
    points_coor1 = np.array(data1.GetPoints().GetData())
    tets1 = np.array(data1.GetCells().GetData())
    cell_locator1 = vtk.vtkCellLocator()
    cell_locator1.SetDataSet(data1)
    cell_locator1.BuildLocator()

    straight_interp_values = np.zeros((num_nodes, num_upsampled))
    deformed_interp_values = np.zeros((num_nodes, num_upsampled))
    straight_interp_values_denoise = np.zeros((num_nodes, num_upsampled))
    deformed_interp_values_denoise = np.zeros((num_nodes, num_upsampled))
    for i in range(num_upsampled):
        start_slot = int(nphi / num_upsampled * i)
        height = (nphi / num_upsampled * i - start_slot) * upper_height
        for n in range(num_nodes):
            pt_coor = [nodes_coor[n][0], nodes_coor[n][1], height]

            cell_id = tet_cellLocator(pt_coor, cell_locator0)
            if not isinstance(cell_id, int):
                cell_id = tet_cellLocator(cell_id, cell_locator0)
            tet = tets0[5 * cell_id + 1:5 * cell_id + 5]
            tet_coor = [points_coor0[node] for node in tet]
            value = []
            for node in tet:
                if node < num_nodes:
                    value.append(values[node, start_slot])
                else:
                    value.append(values[node - num_nodes, (start_slot + 1) % nphi])
            straight_interp_values[n, i] = tet_barycentric_interp(pt_coor, tet_coor, value)
            value = []
            for node in tet:
                if node < num_nodes:
                    value.append(values_denoise[node, start_slot])
                else:
                    value.append(values_denoise[node - num_nodes, (start_slot + 1) % nphi])
            straight_interp_values_denoise[n, i] = tet_barycentric_interp(pt_coor, tet_coor, value)

            cell_id = tet_cellLocator(pt_coor, cell_locator1)
            if not isinstance(cell_id, int):
                cell_id = tet_cellLocator(cell_id, cell_locator1)
            tet = tets1[5 * cell_id + 1:5 * cell_id + 5]
            tet_coor = [points_coor1[node] for node in tet]
            value = []
            for node in tet:
                if node < num_nodes:
                    value.append(values[node, start_slot])
                else:
                    value.append(values[node - num_nodes, (start_slot + 1) % nphi])
            deformed_interp_values[n, i] = tet_barycentric_interp(pt_coor, tet_coor, value)
            value = []
            for node in tet:
                if node < num_nodes:
                    value.append(values_denoise[node, start_slot])
                else:
                    value.append(values_denoise[node - num_nodes, (start_slot + 1) % nphi])
            deformed_interp_values_denoise[n, i] = tet_barycentric_interp(pt_coor, tet_coor, value)
        if i % 20 == 0:
            print(str(i) + "/" + str(num_upsampled))

    with open("straight_interp_values.pickle", "wb") as f:
        pickle.dump(straight_interp_values, f)
    with open("deformed_interp_values.pickle", "wb") as f:
        pickle.dump(deformed_interp_values, f)
    with open("straight_interp_values_denoise.pickle", "wb") as f:
        pickle.dump(straight_interp_values_denoise, f)
    with open("deformed_interp_values_denoise.pickle", "wb") as f:
        pickle.dump(deformed_interp_values_denoise, f)

    # with open("straight_interp_values.pickle", "rb") as f:
    #     straight_interp_values = pickle.load(f)
    # with open("deformed_interp_values.pickle", "rb") as f:
    #     deformed_interp_values = pickle.load(f)
    # with open("straight_interp_values_denoise.pickle", "rb") as f:
    #     straight_interp_values_denoise = pickle.load(f)
    # with open("deformed_interp_values_denoise.pickle", "rb") as f:
    #     deformed_interp_values_denoise = pickle.load(f)

    values_gt = np.zeros((num_nodes, num_upsampled))
    sigma = fwhm2sigma(signal_R)
    for i in range(num_upsampled):
        angle0 = 2 * np.pi / num_upsampled * i
        angle1 = 2 * np.pi / num_upsampled * i + np.pi
        signal_center0 = np.array([signal_cr * np.cos(angle0), signal_cr * np.sin(angle0)])
        signal_center1 = np.array([signal_cr * np.cos(angle1), signal_cr * np.sin(angle1)])
        for n in range(num_nodes):
            value = signal_mag * BivariateGaussian(nodes_coor[n], signal_center0, sigma, sigma, 0)
            values_gt[n, i] = value + signal_mag * BivariateGaussian(nodes_coor[n], signal_center1, sigma, sigma, 0)

    psnr_straight = []
    psnr_deformed = []
    psnr_straight_denoise = []
    psnr_deformed_denoise = []
    for i in range(num_upsampled):
        max_signal = max(values_gt[:, i])
        mse = np.sum((straight_interp_values[:, i] - values_gt[:, i]) ** 2) / num_nodes
        psnr = 20 * np.log10(max_signal) - 10 * np.log10(mse)
        psnr_straight.append(psnr)
        mse = np.sum((deformed_interp_values[:, i] - values_gt[:, i]) ** 2) / num_nodes
        psnr = 20 * np.log10(max_signal) - 10 * np.log10(mse)
        psnr_deformed.append(psnr)
        mse = np.sum((straight_interp_values_denoise[:, i] - values_gt[:, i]) ** 2) / num_nodes
        psnr = 20 * np.log10(max_signal) - 10 * np.log10(mse)
        psnr_straight_denoise.append(psnr)
        mse = np.sum((deformed_interp_values_denoise[:, i] - values_gt[:, i]) ** 2) / num_nodes
        psnr = 20 * np.log10(max_signal) - 10 * np.log10(mse)
        psnr_deformed_denoise.append(psnr)

    with open("psnr_straight.pickle", "wb") as f:
        pickle.dump(psnr_straight, f)
    with open("psnr_deformed.pickle", "wb") as f:
        pickle.dump(psnr_deformed, f)
    with open("psnr_straight_denoise.pickle", "wb") as f:
        pickle.dump(psnr_straight_denoise, f)
    with open("psnr_deformed_denoise.pickle", "wb") as f:
        pickle.dump(psnr_deformed_denoise, f)

    # with open("psnr_straight.pickle", "rb") as f:
    #     psnr_straight = pickle.load(f)
    # with open("psnr_deformed.pickle", "rb") as f:
    #     psnr_deformed = pickle.load(f)
    # with open("psnr_straight_denoise.pickle", "rb") as f:
    #     psnr_straight_denoise = pickle.load(f)
    # with open("psnr_deformed_denoise.pickle", "rb") as f:
    #     psnr_deformed_denoise = pickle.load(f)

    plt.subplot(2, 2, 1)
    plt.plot(range(1024), psnr_straight)
    plt.title("SL (noised)")
    plt.xticks([0, 250, 500, 750, 1000], ["0", "pi/2", "pi", "3*pi/2", "2*pi"])
    plt.ylim(0, 45)
    plt.xlabel("phi")
    plt.ylabel("PSNR")

    plt.subplot(2, 2, 2)
    plt.plot(range(1024), psnr_deformed)
    plt.title("PL (noised)")
    plt.xticks([0, 250, 500, 750, 1000], ["0", "pi/2", "pi", "3*pi/2", "2*pi"])
    plt.ylim(0, 45)
    plt.xlabel("phi")
    plt.ylabel("PSNR")

    plt.subplot(2, 2, 3)
    plt.plot(range(1024), psnr_straight_denoise)
    plt.title("SL")
    plt.xticks([0, 250, 500, 750, 1000], ["0", "pi/2", "pi", "3*pi/2", "2*pi"])
    plt.ylim(0, 80)
    plt.xlabel("phi")
    plt.ylabel("PSNR")

    plt.subplot(2, 2, 4)
    plt.plot(range(1024), psnr_deformed_denoise)
    plt.title("PL")
    plt.xticks([0, 250, 500, 750, 1000], ["0", "pi/2", "pi", "3*pi/2", "2*pi"])
    plt.ylim(0, 80)
    plt.xlabel("phi")
    plt.ylabel("PSNR")
    plt.show()
