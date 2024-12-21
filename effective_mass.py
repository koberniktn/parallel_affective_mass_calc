import copy
import os
import sys
import re
from ortogonal_vertors import basis_vectors
from matplotlib.pyplot import legend
import numpy as np
from obtain_k_array import generate_all_points
import matplotlib.pyplot as plt


st3 = []
st3.append([0.0, 0.0, 0.0]) # 0
st3.append([-1.0, 0.0, 0.0]);  st3.append([1.0, 0.0, 0.0])  # dx  1-2
st3.append([0.0, -1.0, 0.0]);  st3.append([0.0, 1.0, 0.0])  # dy  3-4
st3.append([0.0, 0.0, -1.0]);  st3.append([0.0, 0.0, 1.0])  # dz  5-6
st3.append([-1.0, -1.0, 0.0]); st3.append([1.0, 1.0, 0.0]); st3.append([1.0, -1.0, 0.0]); st3.append([-1.0, 1.0, 0.0]) # dxdy 7-10
st3.append([-1.0, 0.0, -1.0]); st3.append([1.0, 0.0, 1.0]); st3.append([1.0, 0.0, -1.0]); st3.append([-1.0, 0.0, 1.0]) # dxdz 11-14
st3.append([0.0, -1.0, -1.0]); st3.append([0.0, 1.0, 1.0]); st3.append([0.0, 1.0, -1.0]); st3.append([0.0, -1.0, 1.0]) # dydz 15-18


def create_rotation_matrix(a, b, c):
    a_unit = a / np.linalg.norm(a)
    b_unit = b / np.linalg.norm(b)
    c_unit = c / np.linalg.norm(c)

    return np.column_stack([a_unit, b_unit, c_unit])


def generate_kpoints(kpt_frac, h, flag="cartesian", points=None): #generate k-points to mass calculations
    kpt_rec = kpt_frac
    kpoints = []

    if flag=="cartesian": #mode for calc tensor in usial cartesial coordinates
        for i in range(len(st3)):
            k_c = [ kpt_rec[j] + st3[i][j]*h for j in range(3) ]
            kpoints.append([k_c[0], k_c[1], k_c[2]])

    elif flag=="direction": #mode for calc tensor in the particular diretcion
        if not points:
            print("FROM generate_kpoints:\n"
                  "You specified the \"direction\" "
                  "flag, but did not specify the direction for "
                    "calculating the effective mass. You must give "
                    "the coordinates of the start and end point")
            sys.exit(1)

        else:
            basis = basis_vectors(points[1], points[0])
            rotation_matrix = create_rotation_matrix(basis[0], basis[1], basis[2])
            for i in range(len(st)):
                v = np.dot(rotation_matrix, st[i])
                k_c = [kpt_rec[j] + v[j] * h for j in range(3)]
                kpoints.append([k_c[0], k_c[1], k_c[2]])

    return kpath_to_siman(kpoints)


def kpath_to_siman(path_list): #obtain k-point in siman view
    k_path = [len(path_list)]
    for point in path_list:
        elem = (1, *point)
        k_path.append(elem)
    return k_path


def calc_eff_mass_st3(bs, flag, band):
    if flag == 'st3':
        kpoints_indices = range(1, len(st3)+1)
    else:
        sys.exit(1)

    spin = list(bs.bands.keys())[0]  # using spin up
    kis = kpoints_indices

    kp = [bs.kpoints[i - 1].cart_coords for i in kis]
    e = [bs.bands[spin][band - 1][i - 1] / 27.1 for i in kis]
    # print([bs.bands[spin][band - 1][i - 1] for i in kis])

    dh_A = np.linalg.norm(kp[1] - kp[0])  # displacement in k-space
    h = dh_A * 0.529177  # converted to Bohr-1 units

    m = np.zeros([3, 3])

    m[0][0] = (e[1] - 2.0 * e[0] + e[2]) / h ** 2
    m[1][1] = (e[3] - 2.0 * e[0] + e[4]) / h ** 2
    m[2][2] = (e[5] - 2.0 * e[0] + e[6]) / h ** 2

    m[0][1] = (e[7] + e[8] - e[9] - e[10]) / (4.0 * h ** 2)
    m[0][2] = (e[11] + e[12] - e[13] - e[14]) / (4.0 * h ** 2)
    m[1][2] = (e[15] + e[16] - e[17] - e[18]) / (4.0 * h ** 2)

    # symmetrize
    m[1][0] = m[0][1]
    m[2][0] = m[0][2]
    m[2][1] = m[1][2]

    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if m[i][j] != 0:
                m[i][j] = 1 / m[i][j]

    print(m)
    return m


def uniq_elements(x, y):
    import pandas as pd
    df = pd.DataFrame({
        'x': x,
        'y': y
    })

    df = df.drop_duplicates(subset=['x'])
    return df['x'].tolist(), df['y'].tolist()


def parse_EIGENVAL_VASP(path, band, diff2_size=None, debug=False):
    with open(path, 'r') as eigenval_fh:
        for _ in range(5):
            eigenval_fh.readline()

        nelec, nkpt, nband = [int(s) for s in eigenval_fh.readline().split()]
        if debug: print(f'From EIGENVAL: Number of the valence band is {nelec/2} from {nband}\nCurrent band is {band}')
        if band > nband:
            print(f'Requested band {band} is larger than total number of the calculated bands {nband}!')
            sys.exit(1)

        energies = []
        k_points = []

        while eigenval_fh.readline(): # empty line
            k_points.append(eigenval_fh.readline().split()[:3]) # k point coordinates
            for j in range(1, nband+1):
                line = eigenval_fh.readline()
                if band == j:
                    energies.append(float(line.split()[1]))

        if diff2_size and len(energies) != diff2_size:
            print("ATTENTION!! The number of calculated points in EIGENVAL does "
                                                        "not correspond to the selected differential scheme\n"
                                                        "From parse_EIGENVAL_VASP")
            sys.exit(1)

        for i in range(len(k_points)):
            coords = k_points[i]
            k_points[i] = tuple([float(coord) for coord in coords])

    return energies, k_points


def find_extr(y, x, smoothing_window=2, debag=False):
    y_old = copy.copy(y)
    y = np.convolve(y, np.ones(smoothing_window) / smoothing_window, mode='same')
    dy = np.gradient(y) #first central differences derivative (dy[0] = y[1]-y[0]; dy[1]=(y[2]-y[0])/2)
    dy_old = np.gradient(y_old)
    dy[:4], dy[-4:] = dy_old[:4], dy_old[-4:]
    ddy_sign = np.diff(np.sign(dy))

    i = 1
    extr_id = []
    while i < len(x) - 1:
        if ddy_sign[i] != 0:
            if ddy_sign[i - 1] == 0:
                extr_id.append(i)
            if ddy_sign[i] == ddy_sign[i - 1] and extr_id[-1] == i - 1:
                extr_id[-1] = i
        i += 1
    extr_id.extend([0, len(x) - 1])

    extr_x = [tuple(round(coord, 5)  for coord in x[i]) for i in extr_id]
    extr_y = [round(y_old[i], 5) for i in extr_id]

    if debag:
        plt.plot(range(len(y)), y_old)
        plt.scatter(extr_id, extr_y)
        plt.grid()
        plt.show()

    return uniq_elements(extr_x, extr_y)
     







