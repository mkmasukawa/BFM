import numpy as np
import itertools
from mayavi import mlab
import random

def generates_vectors():
    vector_list = [[2, 0, 0], [2, 1, 0], [2, 1, 1], [2, 2, 1], [3, 0, 0], [3, 1, 0]]
    signal_permutation_list = [[1, 1, 1], [1, -1, -1], [1, 1, -1], [1, -1, 1], [-1, 1, 1], [-1, -1, -1], [-1, 1, -1],
                               [-1, -1, 1]]
    complete_vector_list = []

    for i in vector_list:
        for j in list(itertools.permutations(i)):
            for k in signal_permutation_list:
                complete_vector_list.append([a * b for a, b in zip(j, k)])

    vector_list_no_duplicate = []

    for i in complete_vector_list:
        if i not in vector_list_no_duplicate:
            vector_list_no_duplicate.append(i)

    return vector_list_no_duplicate

VECTORS = generates_vectors()

DICT_INT_TO_VECTOR = dict(zip(range(len(VECTORS)), [tuple(l) for l in VECTORS]))

DICT_VECTOR_TO_INT = {v: k for k, v in DICT_INT_TO_VECTOR.items()}

#https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

#function that creates matrix
def plots_bonds(bond_matrix):
    x_len = len(bond_matrix)
    y_len = len(bond_matrix[0])
    z_len = len(bond_matrix[0][0])
    x = np.zeros(x_len)
    y = np.zeros(y_len)
    z = np.zeros(z_len)
    scalars, connections = [], []
    x = []
    for i in range(x_len):
        for j in range(y_len):
            for k in range(z_len):
                a = bond_matrix[i][j][k]
                if a != 0:
                    x.append((i, j, k))
                    scalars.append(1)
                b = [int(x) for x in str(a)]
                split_bonds = chunks(b[:-1], 4)
                for bond in split_bonds:
                    # connections.append([(i, j, k), (i+bond[1], j+bond[2], k+bond[3])])
                    strings = [str(integer) for integer in bond[1:]]
                    a_string = "".join(strings)
                    an_integer = int(a_string)
                    bond_coordinates = DICT_INT_TO_VECTOR[an_integer]
                    connections.append([(i, j, k), (i+bond_coordinates[0], j+bond_coordinates[1], k+bond_coordinates[2])])

    connection_index = []
    for c in connections:
        connection_index.append([x.index(c[0]), x.index(c[1])])

    u, v, w = np.array(x).T
    scalars = np.array(scalars)

    mlab.clf()

    pts = mlab.points3d(u*10, v*10, w*10, 2 * scalars.max() - scalars,
                        scale_factor=10, resolution=7)
    pts.mlab_source.dataset.lines = np.array(connection_index)

    tube = mlab.pipeline.tube(pts, tube_radius=0.0001, tube_sides=4)
    tube.filter.radius_factor = 1
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

    mlab.show()


def update_system(bond_matrix):
    x = list(range(len(bond_matrix)))
    y = list(range(len(bond_matrix[0])))
    z = list(range(len(bond_matrix[0][0])))
    x_random = random.choice(x)
    y_random = random.choice(y)
    z_random = random.choice(z)
    for i in range(200):
        if bond_matrix[x_random, y_random, z_random] != 0:
            break
        else:
            x_random = random.choice(x)
            y_random = random.choice(y)
            z_random = random.choice(z)
    update_position(bond_matrix, x_random, y_random, z_random)

def update_position(bond_matrix, i, j, k):
    if bond_matrix[i, j, k] != 0:
        # choose random vector
        available_positions = []
        #see which of the movements are available
        for v in VECTORS:
            if bond_matrix[i+v[0], j+v[1], k+v[2]] == 0:
                available_positions.append(v)
        chosen_v = random.choice(available_positions)
        #check if the position is available
        #if it is available remake the bonds
        print(available_positions)
    else:
        print("value is zero")


def creates_box(M, N):
    test_m = np.zeros((M, M, M), dtype=int)

    def particle(i, j , k):
        # if 10 < i < M and 10 < j < M and 10 < k < M:
        if 6 < i < M-6 and 6 < j < M-6 and 6 < k < M-6:
            test_m[i, j, k] = 1002106110601
            test_m[i - 2, j + 2, k - 1] = 105710601
            test_m[i - 4, j + 4, k - 2] = 105710601
            test_m[i - 6, j + 6, k - 3] = 10571
            test_m[i - 2, j - 2, k + 1] = 105610611
            test_m[i - 4, j - 4, k + 2] = 105610611
            test_m[i - 6, j - 6, k + 3] = 10561
            test_m[i, j + 2, k] = 100310021
            test_m[i, j + 4, k] = 100310021
            test_m[i, j + 6, k] = 1
        # else:
        #     raise Exception("Values outside matrix")

    x = list(range(len(test_m)))
    y = list(range(len(test_m[0])))
    z = list(range(len(test_m[0][0])))

    for i in range(N):

        x_random = random.choice(x)
        y_random = random.choice(y)
        z_random = random.choice(z)
        particle(x_random, y_random, z_random)
    # print(DICT_VECTOR_TO_INT)
    # print(DICT_VECTOR_TO_INT[(0, -2, 0)])

    # test_m[3, 3, 3] = 12001
    # test_m[5, 3, 3] = 1
    return test_m

    # for i in range(M):
    #     test_m[i]