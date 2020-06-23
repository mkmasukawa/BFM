import numpy as np
import itertools
from mayavi import mlab
import random

random.seed(300)

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

EXCLUSION_VOLUME_VECTORS = []
for a in range(-1, 2):
    for b in range(-1, 2):
        for c in range(-1, 2):
            EXCLUSION_VOLUME_VECTORS.append((a, b , c))

DICT_INT_TO_VECTOR = dict(zip(range(len(VECTORS)), [tuple(l) for l in VECTORS]))

DICT_VECTOR_TO_INT = {v: k for k, v in DICT_INT_TO_VECTOR.items()}




# https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# def vector_to_integer(vector):


# the information of the bonds is stored in the matrix values in the form of an int
# this returns a dictionary
def split_info(a, i, j, k):
    b = [int(x) for x in str(a)]
    monomer_dict = {"type": b[-1]}
    split_bonds = chunks(b[:-1], 4)

    for bond in split_bonds:
        # connections.append([(i, j, k), (i+bond[1], j+bond[2], k+bond[3])])
        strings = [str(integer) for integer in bond[1:]]
        a_string = "".join(strings)
        an_integer = int(a_string)
        bond_coordinates = DICT_INT_TO_VECTOR[an_integer]
        # key is the kind of bond, value is the bond coordinate
        monomer_dict[bond_coordinates] = bond[0]

    return monomer_dict


# function that creates matrix
def plots_bonds(bond_matrix):
    x_len = len(bond_matrix)
    y_len = len(bond_matrix[0])
    z_len = len(bond_matrix[0][0])
    scalars, connections = [], []
    x = []
    for i in range(x_len):
        for j in range(y_len):
            for k in range(z_len):
                a = bond_matrix[i][j][k]
                if a != 0:
                    m_d = split_info(a, i, j, k)
                    m_d.pop("type")
                    for m in m_d:
                        # if (i, j, k) not in x:
                        x.append((i, j, k))
                        x.append((i+m[0], j+m[1], k+m[2]))
                        scalars.append(1)
                        scalars.append(1)

    if x:
        connection_index = []
        for i in range(int(len(x) / 2)):
            connection_index.append([2 * i, 2 * i + 1])

        # for c in connections:
        #     connection_index.append([x.index(c[0]), x.index(c[1])])

        u, v, w = np.array(x).T
        scalars = np.array(scalars)

        mlab.clf()

        pts = mlab.points3d(u * 10, v * 10, w * 10, 2 * scalars.max() - scalars,
                            scale_factor=10, resolution=7)

        pts.mlab_source.dataset.lines = np.array(connection_index)
        tube = mlab.pipeline.tube(pts, tube_radius=0.0001, tube_sides=4)
        tube.filter.radius_factor = 1
        tube.filter.vary_radius = 'vary_radius_by_scalar'
        mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

        mlab.show()

#updates the system within a given spot
#because it looks only at a slice of the matrix, this is why the coordinates need to be relative
def update_system(bond_matrix):
    x = list(range(len(bond_matrix)))
    y = list(range(len(bond_matrix[0])))
    z = list(range(len(bond_matrix[0][0])))
    x_random = random.choice(x)
    y_random = random.choice(y)
    z_random = random.choice(z)
    #tries to find a spot that has a monomer
    for i in range(200):
        if bond_matrix[x_random, y_random, z_random] != 0:
            break
        else:
            x_random = random.choice(x)
            y_random = random.choice(y)
            z_random = random.choice(z)
    update_position(bond_matrix, x_random, y_random, z_random)

#https://stackoverflow.com/questions/20887555/dead-simple-example-of-using-multiprocessing-queue-pool-and-locking
#multiprocessing

# checks the possible movements and returns the movement that is allowed
def chooses_allowed_movement(bond_matrix, i, j, k):
    available_positions = []
    for a in range(-1, 2):
        for b in range(-1, 2):
            for c in range(-1, 2):
                #the new position must be within boundaries
                if 0 <= i+a < len(bond_matrix) and 0 <= j+b < len(bond_matrix[0]) and 0 <= k+c < len(bond_matrix[0][0]):
                    #position must be free
                    if bond_matrix[i + a, j + b, k + c] == 0:
                        invades_exclusion = False
                        # overstretches = False
                        for d in range(-1, 2):
                            for e in range(-1, 2):
                                for f in range(-1, 2):
                                    if 0 <= i + a  + d < len(bond_matrix) and 0 <= j + b + e < len(
                                            bond_matrix[0]) and 0 <= k + c + f < len(bond_matrix[0][0]):
                                        #not consider the monomer being moved
                                        if (a + d + i, b + e + j, c + f + k) != (i, j, k):
                                            #if the neighbor regions are occupied
                                            if bond_matrix[a + d + i, b + e + j, c + f + k] != 0:
                                                invades_exclusion = True
                        if invades_exclusion == False:
                            # if overstretches == False:
                            #needs to check if it overstretches the bond
                            #needs to check if there is bond crossing
                            available_positions.append((a, b, c))
    # remove the possibility that any bond is overstretched
    old_mon_dict = split_info(bond_matrix[i, j, k], i, j, k)
    old_mon_dict.pop("type")
    to_remove = []
    for a_p in available_positions:
        for m in old_mon_dict:
            if (m[0]-a_p[0], m[1]-a_p[1], m[2]-a_p[2] ) not in DICT_VECTOR_TO_INT:
                if a_p in available_positions:
                    to_remove.append(a_p)
    corrected = []
    for a in available_positions:
        if a not in to_remove:
            corrected.append(a)
    if corrected:
        chosen_movement = random.choice(corrected)
        #check if this movement does not overstretch the bonds
        return chosen_movement
    else:
        return False


# updates the position of the matrix returning true if it did, false if no change was made
def update_position(bond_matrix, i, j, k):
    # only positions that have monomers can be modified
    if bond_matrix[i, j, k] != 0:
        # see how many positions are allowed
        move_to = chooses_allowed_movement(bond_matrix, i, j, k)
        if move_to:
            # check all the atoms that are connected to the bonds and modify them too

            #change the atom first
            old_mon_dict = split_info(bond_matrix[i, j, k], i, j, k)
            new_mon = str(old_mon_dict["type"])
            old_mon_dict.pop("type")
            for m in old_mon_dict:
                #redo the bond in the monomer
                new_bond_int = DICT_VECTOR_TO_INT[(m[0] - move_to[0], m[1] - move_to[1], m[2] - move_to[2])]
                new_mon = str(old_mon_dict[m])+ str(new_bond_int).zfill(3) + new_mon
                #redo the bonds in the neighbors
                #neighbor bonds
                # print(bond_matrix[i + m[0], j + m[1], k + m[2]])
                neighbor_dict = split_info(bond_matrix[i + m[0], j + m[1], k + m[2]], i + m[0], j + m[1], k + m[2])
                # print("neighbour_dict", neighbor_dict)
                to_be_changed = (-m[0], -m[1], -m[2])
                # print(DICT_VECTOR_TO_INT[to_be_changed])
                # print(DICT_VECTOR_TO_INT[-m[0]+move_to[0], -m[1]+move_to[1], -m[2]+move_to[2]])
                neighbor_dict[(-m[0]+move_to[0], -m[1]+move_to[1], -m[2]+move_to[2])] = old_mon_dict[m]

                del neighbor_dict[to_be_changed]
                # print(neighbor_dict)
                c = str(neighbor_dict["type"])
                del neighbor_dict["type"]
                #join the bonds again
                for n in neighbor_dict:
                    c = str(neighbor_dict[n]) + str(DICT_VECTOR_TO_INT[n]).zfill(3) + c
                #convert dict back to int
                # print(c)

                # copy_value_neighbor = str(bond_matrix[i + m[0], j + m[1], k + m[2]])
                # to_substitute = str(old_mon_dict[m])+ str(DICT_VECTOR_TO_INT[(-m[0], -m[1], -m[2])]).zfill(3)
                # substitute = str(old_mon_dict[m])+ str(DICT_VECTOR_TO_INT[(-m[0] + move_to[0], -m[1] + move_to[1], -m[2] + move_to[2])]).zfill(3)
                # c = copy_value_neighbor.replace(to_substitute, substitute)
                bond_matrix[i + m[0], j + m[1], k + m[2]] = int(c)
                # print(c)

            bond_matrix[i, j, k] = 0
            bond_matrix[i + move_to[0], j + move_to[1], k + move_to[2]] = new_mon

            return True
        else:
            return False
        # modify position
    else:
        return False


def creates_box(m, n):
    test_m = np.zeros((m, m, m), dtype=int)

    def particle(i, j, k):
        # if 10 < i < M and 10 < j < M and 10 < k < M:
        if 6 < i < m - 6 and 6 < j < m - 6 and 6 < k < m - 6:
            test_m[i, j, k] = 1002106110601
            test_m[i - 2, j + 2, k - 1] = 105710601
            test_m[i - 4, j + 4, k - 2] = 105710601
            test_m[i - 6, j + 6, k - 3] = 10571
            test_m[i - 2, j - 2, k + 1] = 105610611
            test_m[i - 4, j - 4, k + 2] = 105610611
            test_m[i - 6, j - 6, k + 3] = 10561
            test_m[i, j + 2, k] = 100310021
            test_m[i, j + 4, k] = 100310021
            test_m[i, j + 6, k] = 10031
        # else:
        #     raise Exception("Values outside matrix")

    x = list(range(len(test_m)))
    y = list(range(len(test_m[0])))
    z = list(range(len(test_m[0][0])))

    for number in range(n):
        x_random = random.choice(x)
        y_random = random.choice(y)
        z_random = random.choice(z)
        if test_m[x_random, y_random, z_random] == 0:
            particle(x_random, y_random, z_random)
    # print(DICT_VECTOR_TO_INT)
    # print(DICT_VECTOR_TO_INT[(0, -2, 0)])
    #
    # test_m[3, 3, 3] = 12001
    # test_m[5, 3, 3] = 1
    return test_m
