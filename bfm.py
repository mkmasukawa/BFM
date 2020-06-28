import numpy as np
import itertools
from mayavi import mlab
import random
import gc
import shapely.geometry as sg

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
            EXCLUSION_VOLUME_VECTORS.append((a, b, c))

DICT_INT_TO_VECTOR = dict(zip(range(len(VECTORS)), [tuple(l) for l in VECTORS]))

DICT_VECTOR_TO_INT = {v: k for k, v in DICT_INT_TO_VECTOR.items()}


# https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


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


def particle_coordinates(bond_matrix):
    x_len = len(bond_matrix)
    y_len = len(bond_matrix[0])
    z_len = len(bond_matrix[0][0])
    scalars, connections = [], []
    x = []
    colors = []
    for i in range(x_len):
        for j in range(y_len):
            for k in range(z_len):
                cell = bond_matrix[i][j][k]
                if cell != 0:
                    m_d = split_info(cell, i, j, k)
                    monomer_type = m_d["type"]
                    m_d.pop("type")
                    for m in m_d:
                        # if (i, j, k) not in x:
                        x.append((i, j, k))
                        x.append((i + m[0], j + m[1], k + m[2]))
                        # scalars.append(1)
                        if monomer_type == 1:
                            # colors.append((1, 0, 0))
                            scalars.append(0.25)
                        if monomer_type == 2:
                            # colors.append((0, 1, 0))
                            scalars.append(0.5)
                        if monomer_type == 3:
                            scalars.append(0.75)
                            # colors.append((0, 0, 1))
                        scalars.append(0.5)

    if x:
        connection_index = []
        for i in range(int(len(x) / 2)):
            connection_index.append([2 * i, 2 * i + 1])

        # for c in connections:
        #     connection_index.append([x.index(c[0]), x.index(c[1])])

        u, v, w = np.array(x).T
        scalars = np.array(scalars)

        return u, v, w, scalars, connection_index


def animates(bond_matrix, max_steps):
    # View it.
    x, y, z, scalars, connections = particle_coordinates(bond_matrix)
    l = mlab.points3d(x * 10, y * 10, z * 10, 2 * scalars.max() - scalars,
                      scale_factor=10, resolution=7, colormap="copper")
    l.mlab_source.dataset.lines = np.array(connections)
    tube = mlab.pipeline.tube(l, tube_radius=0.0001, tube_sides=4)
    tube.filter.radius_factor = 1
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

    # Now animate the data.
    ms = l.mlab_source

    @mlab.animate(delay=10)
    def anim():
        # f = mlab.gcf()
        mlab.gcf()
        # for i in range(duration):
        i = 0
        # while(True):
        while (i < max_steps):
            print(i)
            x, y, z, s, c = particle_coordinates(bond_matrix)
            ms.reset(x=10 * x, y=10 * y, z=10 * z, scalars=s)
            ms.dataset.lines = np.array(c)

            # l = mlab.points3d(x * 10, y * 10, z * 10, 2 * s.max() - s,
            #                   scale_factor=10, resolution=7, colormap="copper")
            # l.mlab_source.dataset.lines = np.array(c)
            # tube = mlab.pipeline.tube(l, tube_radius=0.0001, tube_sides=4)
            # tube.filter.radius_factor = 1
            # tube.filter.vary_radius = 'vary_radius_by_scalar'
            # mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))
            # f.scene.render()

            # update_system(bond_matrix)
            update_system_with_bonds(bond_matrix)
            i += 1
            # yield chunk.unlink()
            yield
            gc.collect(generation=1)

    anim()
    mlab.show()


# function that creates matrix
def plots_bonds(bond_matrix):
    x_len = len(bond_matrix)
    y_len = len(bond_matrix[0])
    z_len = len(bond_matrix[0][0])
    scalars, connections = [], []
    x = []
    colors = []
    for i in range(x_len):
        for j in range(y_len):
            for k in range(z_len):
                cell = bond_matrix[i][j][k]
                if cell != 0:
                    m_d = split_info(cell, i, j, k)
                    monomer_type = m_d["type"]
                    m_d.pop("type")
                    for m in m_d:
                        # if (i, j, k) not in x:
                        x.append((i, j, k))
                        x.append((i + m[0], j + m[1], k + m[2]))
                        # scalars.append(1)
                        if monomer_type == 1:
                            # colors.append((1, 0, 0))
                            scalars.append(0.25)
                        if monomer_type == 2:
                            # colors.append((0, 1, 0))
                            scalars.append(0.5)
                        if monomer_type == 3:
                            scalars.append(0.75)
                            # colors.append((0, 0, 1))
                        scalars.append(1)

    if x:
        connection_index = []
        for i in range(int(len(x) / 2)):
            connection_index.append([2 * i, 2 * i + 1])

        # for c in connections:
        #     connection_index.append([x.index(c[0]), x.index(c[1])])

        u, v, w = np.array(x).T
        scalars = np.array(scalars)

        # pts = mlab.points3d(u * 10, v * 10, w * 10, 2 * scalars.max() - scalars,
        #                     scale_factor=10,
        #                     resolution=7, color=colors)

        pts = mlab.points3d(u * 10, v * 10, w * 10, 2 * scalars.max() - scalars,
                            scale_factor=10, resolution=7, colormap="copper")

        pts.mlab_source.dataset.lines = np.array(connection_index)
        tube = mlab.pipeline.tube(pts, tube_radius=0.0001, tube_sides=4)
        tube.filter.radius_factor = 1
        tube.filter.vary_radius = 'vary_radius_by_scalar'
        mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

        # print(type(mlab.screenshot()))
        mlab.show()
        # return mlab.screenshot()
        # mlab.savefig("equilibrated" + '.png')
        # mlab.savefig("test.vtk")
        # pts.module_manager.source.save_output('output.vtk')


# TODO: make simpified plot of the whole system if too heavy to be loaded by mayavi


# updates the system within a given spot
# because it looks only at a slice of the matrix, this is why the coordinates need to be relative
def update_system(bond_matrix):
    filled_cells = np.nonzero(bond_matrix)
    if filled_cells:
        chosen = random.choice(np.transpose(filled_cells))
        update_position(bond_matrix, chosen[0], chosen[1], chosen[2])
    return bond_matrix


# converts representations of monomers from dictionary to int (opposite of split function)
def dict_to_int(monomer_dict):
    try:
        s = str(monomer_dict["type"])
        del monomer_dict["type"]
        if monomer_dict:
            for m in monomer_dict:
                s = str(monomer_dict[m]) + str(DICT_VECTOR_TO_INT[m]).zfill(3) + s
            return int(s)
        else:
            raise Exception("Dictionary can not be converted to int, missing bonds")
    except AssertionError as error:
        print(error)
        raise Exception("Dictionary canot be converted to int, missing the monomer type")


# checks if two lin segments intersect
# the segments intersect in 3d if they intersect on all planes xy, yz, zx
def segments_intersect(segment1, segment2):
    segment1_point1, segment1_point2 = segment1
    segment2_point1, segment2_point2 = segment2
    # intersects on plane xy
    line1 = sg.LineString([(segment1_point1[0], segment1_point1[1]), (segment1_point2[0], segment1_point2[1])])
    line2 = sg.LineString([(segment2_point1[0], segment2_point1[1]), (segment2_point2[0], segment2_point2[1])])

    # intersects on plane yz
    line3 = sg.LineString([(segment1_point1[1], segment1_point1[2]), (segment1_point2[1], segment1_point2[2])])
    line4 = sg.LineString([(segment2_point1[1], segment2_point1[2]), (segment2_point2[1], segment2_point2[2])])
    # intersects on plane xz
    line5 = sg.LineString([(segment1_point1[0], segment1_point1[2]), (segment1_point2[0], segment1_point2[2])])
    line6 = sg.LineString([(segment2_point1[0], segment2_point1[2]), (segment2_point2[0], segment2_point2[2])])

    if line1.intersects(line2) and line3.intersects(line4) and line5.intersects(line6):
        return True
    else:
        return False


# checks if a bonds crosses nearby bonds
# point1 exists, point2 is the one whose creation is considered
def bond_cross_bond(bond_matrix, point1, point2):
    # looks for bonds around point2

    list_bonds_to_examinate = []
    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                i, j, k = a + point2[0], b + point2[1], c + point2[2]
                if 0 <= i < len(bond_matrix) and 0 <= j < len(bond_matrix[0]) and 0 <= k < len(bond_matrix[0][0]):
                    dict_monomer = split_info(bond_matrix[i, j, k], i, j, k)
                    del dict_monomer["type"]
                    for d_c in dict_monomer:
                        # do not consider bonds that include point1 or point2
                        if (i, j, k) != point1 and (i + d_c[0], j + d_c[1], k + d_c[2]) != point1:
                            list_bonds_to_examinate.append(((i, j, k), (i + d_c[0], j + d_c[1], k + d_c[2])))

    for l_e in list_bonds_to_examinate:
        if segments_intersect((point1, point2), l_e):
        # if segment_intersects_3D_plane((point1, point2), l_e):
            return True
    return False


def update_positions_bonds(bond_matrix, i, j, k):
    class Found(Exception):
        pass

    # only positions that have monomers can be modified
    if bond_matrix[i, j, k] != 0:
        # see how many positions are allowed
        move_to = chooses_allowed_movement(bond_matrix, i, j, k)

        # if there are positions to move to
        if move_to:
            # check all the atoms that are connected to the bonds and modify them too

            # change the atom first
            old_mon_dict = split_info(bond_matrix[i, j, k], i, j, k)
            new_mon = str(old_mon_dict["type"])
            old_mon_dict.pop("type")
            for m in old_mon_dict:
                # redo the bond in the monomer
                new_bond_int = DICT_VECTOR_TO_INT[(m[0] - move_to[0], m[1] - move_to[1], m[2] - move_to[2])]
                new_mon = str(old_mon_dict[m]) + str(new_bond_int).zfill(3) + new_mon
                # redo the bonds in the neighbors
                # neighbor bonds
                # print(bond_matrix[i + m[0], j + m[1], k + m[2]])
                neighbor_dict = split_info(bond_matrix[i + m[0], j + m[1], k + m[2]], i + m[0], j + m[1], k + m[2])
                # print("neighbour_dict", neighbor_dict)
                to_be_changed = (-m[0], -m[1], -m[2])
                # change also the bonds of the neighbors
                neighbor_dict[(-m[0] + move_to[0], -m[1] + move_to[1], -m[2] + move_to[2])] = old_mon_dict[m]
                del neighbor_dict[to_be_changed]
                # print(neighbor_dict)
                c_m = str(neighbor_dict["type"])
                del neighbor_dict["type"]
                # join the bonds again
                for n in neighbor_dict:
                    c_m = str(neighbor_dict[n]) + str(DICT_VECTOR_TO_INT[n]).zfill(3) + c_m
                # convert dict back to int
                bond_matrix[i + m[0], j + m[1], k + m[2]] = int(c_m)

            bond_matrix[i, j, k] = 0
            bond_matrix[i + move_to[0], j + move_to[1], k + move_to[2]] = new_mon

            # after moving the atom, check if it is close to any monomer which is a type 2
            # forming bond
            # if two type 2 are nearby, then the types 3 they are attached to are removed and they are connected,
            # remember to do the bond both ways
            new_i = i + move_to[0]
            new_j = j + move_to[1]
            new_k = k + move_to[2]
            current_mon = split_info(bond_matrix[new_i, new_j, new_k], new_i, new_j, new_k)

            # the monomer that is being moved is type 2
            if current_mon["type"] == 2:
                del current_mon["type"]
                # check if it is bonded to type 2 already
                # calculate here probability to break
                bonded_to_type_2 = False
                bonded_to_coordinates = (0, 0, 0)
                for b_c_m in current_mon:
                    bonded_to_i = new_i + b_c_m[0]
                    bonded_to_j = new_j + b_c_m[1]
                    bonded_to_k = new_k + b_c_m[2]
                    bonded_to = split_info(bond_matrix[bonded_to_i, bonded_to_j, bonded_to_k], bonded_to_i,
                                           bonded_to_j, bonded_to_k)
                    if bonded_to["type"] == 2:
                        bonded_to_type_2 = True
                        bonded_to_coordinates = (new_i + b_c_m[0], new_j + b_c_m[1], new_k + b_c_m[2])

                # can not bond again to type 2 if already bonded, but can break the bond
                if bonded_to_type_2 == True:
                    print("already bonded to type 2")
                    # calculate probability to break
                    # check if there are free spaces around the double bond to create the type 3 monomers
                    spaces_around_new = set()
                    spaces_around_bonded_to = set()
                    for v in VECTORS:
                        if 0 <= new_i + v[0] < len(bond_matrix) \
                                and 0 <= new_j + v[1] < len(bond_matrix[0]) \
                                and 0 <= new_k + v[2] < len(bond_matrix[0][0]):
                            if bond_matrix[new_i + v[0], new_j + v[1], new_k + v[2]] == 0:
                                spaces_around_new.add((new_i + v[0], new_j + v[1], new_k + v[2]))
                        if 0 <= bonded_to_coordinates[0] + v[0] < len(bond_matrix) \
                                and 0 <= bonded_to_coordinates[1] + v[1] < len(bond_matrix[0]) \
                                and 0 <= bonded_to_coordinates[2] + v[2] < len(bond_matrix[0][0]):
                            if bond_matrix[bonded_to_coordinates[0] + v[0], bonded_to_coordinates[1] + v[1],
                                           bonded_to_coordinates[2] + v[2]] == 0:
                                spaces_around_bonded_to.add((bonded_to_coordinates[0] + v[0],
                                                             bonded_to_coordinates[1] + v[1],
                                                             bonded_to_coordinates[2] + v[2]))
                    # check if in the two lists there are places for the type 3 monomers
                    # check values in common, if they are more than 2 there is space
                    # if there is less than 2, substract the common set and see if th
                    inter = spaces_around_new.intersection(spaces_around_bonded_to)
                    new_exclusive = spaces_around_new.difference(spaces_around_bonded_to)
                    bonded_to_exclusive = spaces_around_bonded_to.difference(spaces_around_new)
                    there_is_space_for_breaking_bond = False
                    if len(inter) >= 2:
                        # there is space, break bonds
                        there_is_space_for_breaking_bond = True
                    else:
                        if len(new_exclusive) >= 1 or len(bonded_to_exclusive) >= 1:
                            # there is space, break bonds
                            there_is_space_for_breaking_bond = True
                        else:
                            # no space
                            there_is_space_for_breaking_bond = False
                    if there_is_space_for_breaking_bond == True:
                        # calculate probability

                        #checks if cross other bonds, if crosses, gives up on the movement


                        # #if there are no options that can not create new bonds,
                        # type_3_for_new = random.choice(tuple(spaces_around_new))
                        # if type_3_for_new in spaces_around_bonded_to:
                        #     spaces_around_bonded_to.remove(type_3_for_new)
                        # type_3_for_bonded_to = random.choice(tuple(spaces_around_bonded_to))
                        #
                        # #does not cross bonds
                        # if (bond_cross_bond(bond_matrix, (new_i, new_j, new_k), type_3_for_new)== False) and \
                        #         (bond_cross_bond(bond_matrix, bonded_to_coordinates, type_3_for_bonded_to) == False):

                        #if there are no options that can not create new bonds,
                        there_is_non_crossing_option = False

                        for_new_list = list(spaces_around_new)
                        random.shuffle(for_new_list)
                        for_bonded_list = list(spaces_around_bonded_to)
                        random.shuffle(for_bonded_list)

                        type_3_for_new = (0,0,0)
                        type_3_for_bonded_to = (0,0,0)

                        try:
                            for n_l in for_new_list:
                                for b_l in for_bonded_list:
                                    if n_l != b_l:
                                        if bond_cross_bond(bond_matrix, (i, j, k), n_l) == False:
                                            if bond_cross_bond(bond_matrix, (new_i, new_j, new_k), b_l) == False:
                                                raise Found
                        except Found:
                            type_3_for_new = n_l
                            type_3_for_bonded_to = b_l

                        #does not cross bonds
                        if type_3_for_new != (0,0,0) and type_3_for_bonded_to != (0,0,0):

                            # create bond from 2 to 3 type monomer and break bond 2-2
                            # add the type of monomer
                            current_mon["type"] = 2
                            # deletes de 2-2 bond
                            del current_mon[(bonded_to_coordinates[0] - new_i,
                                             bonded_to_coordinates[1] - new_j,
                                             bonded_to_coordinates[2] - new_k)]
                            # creates the 2-3 bond
                            current_mon[
                                (type_3_for_new[0] - new_i, type_3_for_new[1] - new_j, type_3_for_new[2] - new_k)] = 3

                            bond_matrix[new_i, new_j, new_k] = dict_to_int(current_mon)

                            bonded_to_dict = split_info(bond_matrix[bonded_to_coordinates], bonded_to_coordinates[0],
                                                        bonded_to_coordinates[1], bonded_to_coordinates[2])
                            # deletes the 2-2 bond
                            del bonded_to_dict[(-bonded_to_coordinates[0] + new_i, -bonded_to_coordinates[1] + new_j,
                                                -bonded_to_coordinates[2] + new_k)]
                            # print("check difference", type_3_for_bonded_to, bonded_to_coordinates)
                            # creates the 2-3 bond

                            # why this had to be deleted??
                            bonded_to_dict[(type_3_for_bonded_to[0] - bonded_to_coordinates[0],
                                            type_3_for_bonded_to[1] - bonded_to_coordinates[1],
                                            type_3_for_bonded_to[2] - bonded_to_coordinates[2])] = 3
                            bond_matrix[bonded_to_coordinates] = dict_to_int(bonded_to_dict)
                            print("finished unbonding")
                            # remove the  bond 2-2

                            # create bond from 3 to 2 type monomer
                            # print( int(str(3)  + str(DICT_VECTOR_TO_INT[(new_i - type_3_for_new[0],
                            #                                           new_j - type_3_for_new[1],
                            #                                           new_k - type_3_for_new[2])]).zfill(3) + str(3)))
                            bond_matrix[type_3_for_new] = \
                                int(str(3) + str(DICT_VECTOR_TO_INT[(new_i - type_3_for_new[0],
                                                                     new_j - type_3_for_new[1],
                                                                     new_k - type_3_for_new[2])]).zfill(3) + str(3))
                            # print(int(str(3)  + str(DICT_VECTOR_TO_INT[(bonded_to_coordinates[0] - type_3_for_bonded_to[0],
                            #                                           bonded_to_coordinates[1] - type_3_for_bonded_to[1],
                            #                                           bonded_to_coordinates[2] - type_3_for_bonded_to[2])]).zfill(3) + str(3)))
                            bond_matrix[type_3_for_bonded_to] = \
                                int(str(3) + str(DICT_VECTOR_TO_INT[(bonded_to_coordinates[0] - type_3_for_bonded_to[0],
                                                                     bonded_to_coordinates[1] - type_3_for_bonded_to[1],
                                                                     bonded_to_coordinates[2] - type_3_for_bonded_to[
                                                                         2])]).zfill(3) + str(3))
                        else:
                            print("All type 3 options crosses bonds")


                # not bonded to type 2 yet
                # need to check if the candidate is not bonded already
                else:
                    # calculate probability to bond
                    # check is there is a type 2 nearby
                    neighbors_type2 = []
                    for v in VECTORS:
                        neighbor_i = new_i + v[0]
                        neighbor_j = new_j + v[1]
                        neighbor_k = new_k + v[2]
                        if 0 <= neighbor_i < len(bond_matrix) and 0 <= neighbor_j < len(bond_matrix[0]) and \
                                0 <= neighbor_k < len(bond_matrix[0][0]):
                            neighbor_mon = split_info(bond_matrix[neighbor_i, neighbor_j, neighbor_k], neighbor_i,
                                                      neighbor_j, neighbor_k)
                            if neighbor_mon["type"] == 2:
                                del neighbor_mon["type"]
                                bonded_to_type_2 = False
                                for b_n_m in neighbor_mon:
                                    if str(bond_matrix[
                                               neighbor_i + b_n_m[0], neighbor_j + b_n_m[1],
                                               neighbor_k + b_n_m[2]])[-1] == "2":
                                        bonded_to_type_2 = True
                                if bonded_to_type_2 == False:
                                    neighbors_type2.append((neighbor_i, neighbor_j, neighbor_k))
                                    print("there is available neighbor type 2")
                            # if neighbor_mon["type"] == 3:
                            #     neighbors_type3.append((neighbor_i, neighbor_j, neighbor_k))
                            #     print("there is neighbor type 3")
                            # make a bond
                            # maybe there is more than one neighbor type 2, in this case one is randomly chosen
                    # there are neigbors type 2 that can be bonded to
                    if neighbors_type2:
                        chosen_type2 = random.choice(neighbors_type2)
                        # remove the monomers type 3 connected to

                        print(DICT_VECTOR_TO_INT)

                        # save from dict to the matrix
                        bonded_current = str(3) + str(DICT_VECTOR_TO_INT[(chosen_type2[0] - new_i,
                                                                          chosen_type2[1] - new_j,
                                                                          chosen_type2[2] - new_k)]) + str(2)
                        for c_m in current_mon:
                            if str(bond_matrix[new_i + c_m[0], new_j + c_m[1], new_k + c_m[2]])[-1] == "3":
                                bond_matrix[new_i + c_m[0], new_j + c_m[1], new_k + c_m[2]] = 0
                            else:
                                bonded_current = str(current_mon[c_m]) + str(DICT_VECTOR_TO_INT[c_m]).zfill(3) \
                                                 + bonded_current
                                # remove also from the perspective of type 3

                        # print(bonded_current)
                        bond_matrix[new_i, new_j, new_k] = int(bonded_current)
                        bond_to_dict = split_info(bond_matrix[chosen_type2[0], chosen_type2[1],
                                                              chosen_type2[2]], chosen_type2[0],
                                                  chosen_type2[1], chosen_type2[2])
                        del bond_to_dict["type"]
                        bonded_neighbor = str(3) + str(DICT_VECTOR_TO_INT[(-chosen_type2[0] + new_i,
                                                                           -chosen_type2[1] + new_j,
                                                                           -chosen_type2[2] + new_k)]) + str(2)
                        for d in bond_to_dict:
                            # print(chosen_type2, d)
                            if str(bond_matrix[chosen_type2[0] + d[0], chosen_type2[1] + d[1],
                                               chosen_type2[2] + d[2]])[-1] == "3":
                                bond_matrix[chosen_type2[0] + d[0], chosen_type2[1] + d[1], chosen_type2[2] + d[2]] = 0
                            else:
                                bonded_neighbor = str(bond_to_dict[d]) \
                                                  + str(DICT_VECTOR_TO_INT[d]).zfill(3) + bonded_neighbor
                        # print(bonded_neighbor)
                        # save dict to the matrix
                        bond_matrix[chosen_type2] = int(bonded_neighbor)
                        # delete the other monomers

                        # break bonds to type 3 of the types 2 that were bonded

                # breaking bond
                # if a bond is broken between 2 type 2, check if there is space for creating type 3 monomers that c
                # an be connected to them
                # is yes, then break the bond, create the new monomers type 3 and bind them to the respective type 3,
                # remember to do the bond both ways

            return True

        # there are no positions to move to
        else:
            return False

    # in case the chosen site has no monomer
    else:
        return False


# updates system for moving monomers but also bond formation
def update_system_with_bonds(bond_matrix):
    filled_cells = np.nonzero(bond_matrix)
    if filled_cells:
        chosen = random.choice(np.transpose(filled_cells))
        update_positions_bonds(bond_matrix, chosen[0], chosen[1], chosen[2])
    return bond_matrix


# https://stackoverflow.com/questions/20887555/dead-simple-example-of-using-multiprocessing-queue-pool-and-locking
# multiprocessing

# checks the possible movements and returns the movement that is allowed
def chooses_allowed_movement(bond_matrix, i, j, k):
    available_positions = []
    for a in range(-1, 2):
        for b in range(-1, 2):
            for c in range(-1, 2):
                # the new position must be within boundaries
                if 0 <= i + a < len(bond_matrix) and 0 <= j + b < len(bond_matrix[0]) and 0 <= k + c < len(
                        bond_matrix[0][0]):
                    # position must be free
                    if bond_matrix[i + a, j + b, k + c] == 0:
                        invades_exclusion = False
                        # overstretches = False
                        for d in range(-1, 2):
                            for e in range(-1, 2):
                                for f in range(-1, 2):
                                    if 0 <= i + a + d < len(bond_matrix) and 0 <= j + b + e < len(
                                            bond_matrix[0]) and 0 <= k + c + f < len(bond_matrix[0][0]):
                                        # not consider the monomer being moved
                                        if (a + d + i, b + e + j, c + f + k) != (i, j, k):
                                            # if the neighbor regions are occupied
                                            if bond_matrix[a + d + i, b + e + j, c + f + k] != 0:
                                                invades_exclusion = True
                        if invades_exclusion == False:
                            # if overstretches == False:
                            # needs to check if it overstretches the bond
                            # needs to check if there is bond crossing
                            # appends the relative coordinates
                            available_positions.append((a, b, c))
    # remove the possibility that any bond is overstretched
    old_mon_dict = split_info(bond_matrix[i, j, k], i, j, k)
    old_mon_dict.pop("type")
    to_remove = []
    for a_p in available_positions:
        for m in old_mon_dict:
            if (m[0] - a_p[0], m[1] - a_p[1], m[2] - a_p[2]) not in DICT_VECTOR_TO_INT:
                if a_p in available_positions:
                    to_remove.append(a_p)
    corrected = []
    for a in available_positions:
        if a not in to_remove:
            corrected.append(a)
    if corrected:
        # chosen_movement = random.choice(corrected)
        # makes the list order random
        random.shuffle(corrected)
        for c in corrected:
            # if the chosen monomer does not crosses bonds
            if crosses_bond(bond_matrix, (i, j, k), c) == False:
                return c
        # make a list of bonds to check if they cross when the monomer moves
        # but do the bond crossing can cause artificial entanglement if it is in 3D?
        list_of_bonds_to_check = []
        # check if this movement does not cross any bonds

        return False
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

            # change the atom first
            old_mon_dict = split_info(bond_matrix[i, j, k], i, j, k)
            new_mon = str(old_mon_dict["type"])
            old_mon_dict.pop("type")
            for m in old_mon_dict:
                # redo the bond in the monomer
                new_bond_int = DICT_VECTOR_TO_INT[(m[0] - move_to[0], m[1] - move_to[1], m[2] - move_to[2])]
                new_mon = str(old_mon_dict[m]) + str(new_bond_int).zfill(3) + new_mon
                # redo the bonds in the neighbors
                # neighbor bonds
                # print(bond_matrix[i + m[0], j + m[1], k + m[2]])
                neighbor_dict = split_info(bond_matrix[i + m[0], j + m[1], k + m[2]], i + m[0], j + m[1], k + m[2])
                # print("neighbour_dict", neighbor_dict)
                to_be_changed = (-m[0], -m[1], -m[2])
                # print(DICT_VECTOR_TO_INT[to_be_changed])
                # print(DICT_VECTOR_TO_INT[-m[0]+move_to[0], -m[1]+move_to[1], -m[2]+move_to[2]])
                neighbor_dict[(-m[0] + move_to[0], -m[1] + move_to[1], -m[2] + move_to[2])] = old_mon_dict[m]
                # print(m)
                # print(neighbor_dict)
                del neighbor_dict[to_be_changed]
                # print(neighbor_dict)
                c = str(neighbor_dict["type"])
                del neighbor_dict["type"]
                # join the bonds again
                for n in neighbor_dict:
                    c = str(neighbor_dict[n]) + str(DICT_VECTOR_TO_INT[n]).zfill(3) + c
                # convert dict back to int
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


# function that takes the absoulute coordinate of the bond being moved, the relative motion and says if it crosses bonds
def crosses_bond(bond_matrix, monomer, new_monomer):
    # make a list with all bond in the vicinity
    i, j, k = monomer
    new_i, new_j, new_k = new_monomer
    # list of absolute bond values
    list_of_bonds_to_examinate = []
    list_of_bond_monomer_triangle = []

    for a in range(-4, 5):
        for b in range(-4, 5):
            for c in range(-4, 5):
                if bond_matrix[a, b, c] != 0:
                    dict_mon = split_info(bond_matrix[a, b, c], a, b, c)
                    del dict_mon["type"]
                    for d_m in dict_mon:
                        # bond = ((a, b, c), (a + d_m[0], b + d_m[1], c + d_m[2]))
                        if (a, b, c) == (i, j, k):
                            # all the bonds of the monomer being observed
                            list_of_bond_monomer_triangle.append(
                                ((i, j, k), (new_i, new_j, new_k), (a + d_m[0], b + d_m[1], c + d_m[2])))
                        else:
                            # do not conside the monomer being observed
                            if (a + d_m[0], b + d_m[1], c + d_m[2]) != (i, j, k):
                                bond1 = ((a, b, c), (a + d_m[0], b + d_m[1], c + d_m[2]))
                                bond2 = ((a + d_m[0], b + d_m[1], c + d_m[2]), (a, b, c))
                                # avoid counting same bond twice
                                if (bond1 not in list_of_bonds_to_examinate) and (
                                        bond2 not in list_of_bonds_to_examinate):
                                    list_of_bonds_to_examinate.append(bond1)
    for l_e in list_of_bonds_to_examinate:
        for l_t in list_of_bond_monomer_triangle:
            if segment_intersects_3D_plane(l_e, l_t) == True:
                # found one bond that crosses the triangle represented by the bond movement
                return True
    # movement causes no bond crossing
    return False


def invades_exclusion_space(bond_matrix, i, j, k):
    def test_position(w, u, v):
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c in range(-1, 2):
                    # there is an exclusion space around the borders
                    # print(w + a, u + b, c + v)
                    # print(len(bond_matrix), len(bond_matrix[0][0]), len(bond_matrix[0][0]))
                    if not ((0 <= w + a < len(bond_matrix)) and (0 <= u + b < len(bond_matrix[0])) and (
                            0 <= v + c < len(bond_matrix[0][0]))):
                        return False
                    if bond_matrix[w + a, u + b, v + c] != 0:
                        return False
        return True

    if not test_position(i, j, k):
        return False
    if not test_position(i - 2, j + 2, k - 1):
        return False
    if not test_position(i - 4, j + 4, k - 2):
        return False
    if not test_position(i - 6, j + 6, k - 3):
        return False
    if not test_position(i - 2, j - 2, k + 1):
        return False
    if not test_position(i - 4, j - 4, k + 2):
        return False
    if not test_position(i - 6, j - 6, k + 3):
        return False
    if not test_position(i, j + 2, k):
        return False
    if not test_position(i, j + 4, k):
        return False
    if not test_position(i, j + 6, k):
        return False
    return True


def plot_shapely(segment, plane):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    u, v, w = np.array(segment).T
    ax.plot(u, v, w, color='b')

    x, y, z = np.array(plane).T
    ax.plot(x, y, z, color='k')
    # ax.plot(dates, zaxisvalues2, upperLimits, color='r')

    # for i, j, k, h in zip(dates, zaxisvalues0, lows, highs):
    #     ax.plot([i, i], [j, j], [k, h], color='g')
    #
    # ax.scatter(dates, zaxisvalues0, highs, color='g', marker="o")
    # ax.scatter(dates, zaxisvalues0, lows, color='y', marker="^")

    plt.show()


# http://geomalgorithms.com/a05-_intersect-1.html
# https://stackoverflow.com/questions/47359985/shapely-intersection-point-between-line-and-polygon-in-3d
# https://stackoverflow.com/questions/53962225/how-to-know-if-a-line-segment-intersects-a-triangle-in-3d-space
# https://stackoverflow.com/questions/53962225/how-to-know-if-a-line-segment-intersects-a-triangle-in-3d-space
# segment and plane are composed of lists of 3D points
def segment_intersects_3D_plane(segment, plane):
    plane_a = np.append(np.array(plane[0]), 1.)
    plane_b = np.append(np.array(plane[1]), 1.)
    plane_c = np.append(np.array(plane[2]), 1.)
    s_d = np.append(np.array(segment[0]), 1.)
    s_e = np.append(np.array(segment[1]), 1.)
    vol_Td = 1. / 6. * np.linalg.det((plane_a, plane_b, plane_c, s_d))
    vol_Te = 1. / 6. * np.linalg.det((plane_a, plane_b, plane_c, s_e))

    # plot_shapely(segment, plane)
    if vol_Td == 0. or vol_Te == 0.:
        # considers passing by a monomer is also forbidden
        return True
    else:
        if vol_Td * vol_Te > 0:
            # both points over or under plane
            return False
        else:
            vol1 = 1 / 6 * np.linalg.det((plane_a, plane_b, s_d, s_e))
            vol2 = 1 / 6 * np.linalg.det((plane_b, plane_c, s_d, s_e))
            vol3 = 1 / 6 * np.linalg.det((plane_c, plane_a, s_d, s_e))
            # check if all have the same sign
            if (vol1 > 0 and vol2 > 0 and vol3 > 0) or (vol1 < 0 and vol2 < 0 and vol3 < 0):
                # segment intersects triangle
                return True
            else:
                return False
