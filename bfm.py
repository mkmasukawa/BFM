import numpy as np
import itertools
from mayavi import mlab
import random
import gc

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
        while(i < max_steps):
            print(i)
            x, y , z, s, c = particle_coordinates(bond_matrix)
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


#TODO: make simpified plot of the whole system if too heavy to be loaded by mayavi



# updates the system within a given spot
# because it looks only at a slice of the matrix, this is why the coordinates need to be relative
def update_system(bond_matrix):
    filled_cells = np.nonzero(bond_matrix)
    if filled_cells:
        chosen = random.choice(np.transpose(filled_cells))
        update_position(bond_matrix, chosen[0], chosen[1], chosen[2])
    return bond_matrix


def update_positions_bonds(bond_matrix, i, j, k):
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

            #after moving the atom, check if it is close to any monomer which is a type 2
            #forming bond
            #if two type 2 are nearby, then the types 3 they are attached to are removed and they are connected,
            # remember to do the bond both ways
            new_i = i + move_to[0]
            new_j = j + move_to[1]
            new_k = k + move_to[2]
            current_mon = split_info(bond_matrix[new_i, new_j, new_k],  new_i, new_j, new_k )
            if current_mon["type"] == 2:
                del current_mon["type"]
                #check if it is bonded to type 2 already
                #calculate here probability to break
                bonded_to_type_2 = False
                for b_c_m in current_mon:
                    bonded_to_i = new_i + b_c_m[0]
                    bonded_to_j = new_j + b_c_m[1]
                    bonded_to_k = new_k + b_c_m[2]
                    bonded_to = split_info(bond_matrix[bonded_to_i, bonded_to_j, bonded_to_k], bonded_to_i,
                                           bonded_to_j, bonded_to_k)
                    if bonded_to["type"] == 2:
                        bonded_to_type_2 = True

                #can not bond again to type 2 if already bonded
                if bonded_to_type_2 == True:
                    print("already bonded to type 2")
                    # calculate probability to break
                    #check if there are free spaces around the double bond to create the type 3 monomers
                    # spaces_around_new =
                    #break the bond

                    #create the new type 3 monomers


                #not bonded to type 2 yet
                #need to check if the candidate is not bonded already
                else:
                    #calculate probability to bond
                    #check is there is a type 2 nearby
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
                                    if str(bond_matrix[neighbor_i + b_n_m[0], neighbor_j + b_n_m[1], neighbor_k + b_n_m[2]])[-1] == "2":
                                        bonded_to_type_2 = True
                                if bonded_to_type_2 == False:
                                    neighbors_type2.append((neighbor_i, neighbor_j, neighbor_k))
                                    print("there is available neighbor type 2")
                            # if neighbor_mon["type"] == 3:
                            #     neighbors_type3.append((neighbor_i, neighbor_j, neighbor_k))
                            #     print("there is neighbor type 3")
                            #make a bond
                            #maybe there is more than one neighbor type 2, in this case one is randomly chosen
                    #there are neigbors type 2 that can be bonded to
                    if neighbors_type2:
                        chosen_type2 = random.choice(neighbors_type2)
                        #remove the monomers type 3 connected to

                        print(DICT_VECTOR_TO_INT)

                        #save from dict to the matrix
                        bonded_current = str(3) + str(DICT_VECTOR_TO_INT[(chosen_type2[0] - new_i,
                                                                          chosen_type2[1] - new_j,
                                                                          chosen_type2[2] - new_k)] )+ str(2)
                        for c_m in current_mon:
                            if str(bond_matrix[new_i+c_m[0], new_j+c_m[1], new_k+c_m[2]])[-1] == "3":
                                bond_matrix[new_i + c_m[0], new_j + c_m[1], new_k + c_m[2]] = 0
                            else:
                                bonded_current = str(current_mon[c_m]) + str(DICT_VECTOR_TO_INT[c_m]).zfill(3) \
                                                 + bonded_current
                                #remove also from the perspective of type 3

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
                            if str(bond_matrix[chosen_type2[0]+d[0], chosen_type2[1]+d[1],
                                               chosen_type2[2]+d[2]])[-1] == "3":
                                bond_matrix[chosen_type2[0] + d[0], chosen_type2[1] + d[1], chosen_type2[2] + d[2]] = 0
                            else:
                                bonded_neighbor = str(bond_to_dict[d]) \
                                                  + str(DICT_VECTOR_TO_INT[d]).zfill(3) + bonded_neighbor
                        # print(bonded_neighbor)
                        #save dict to the matrix
                        bond_matrix[chosen_type2] = int(bonded_neighbor)
                        #delete the other monomers


                                #break bonds to type 3 of the types 2 that were bonded

                #breaking bond
                #if a bond is broken between 2 type 2, check if there is space for creating type 3 monomers that c
                # an be connected to them
                #is yes, then break the bond, create the new monomers type 3 and bind them to the respective type 3,
                # remember to do the bond both ways


            return True
        else:
            return False
    else:
        return False

#updates system for moving monomers but also bond formation
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
        chosen_movement = random.choice(corrected)
        # check if this movement does not overstretch the bonds
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


def creates_box(m, n):
    test_m = np.zeros((m, m, m), dtype=int)

    # avoid superposition of monomers during creation
    def particle(i, j, k):
        # it should be within limits of the matrix
        if 6 < i < m - 6 and 6 < j < m - 6 and 6 < k < m - 6:
            # it should not invade the exclusion space of other polymers
            # if invades_exclusion_space
            if invades_exclusion_space(test_m, i, j, k):
                test_m[i, j, k] = 1002106110601
                test_m[i - 2, j + 2, k - 1] = 105710601
                test_m[i - 4, j + 4, k - 2] = 105740602
                test_m[i - 6, j + 6, k - 3] = 40573
                test_m[i - 2, j - 2, k + 1] = 105610611
                test_m[i - 4, j - 4, k + 2] = 105640612
                test_m[i - 6, j - 6, k + 3] = 40563
                test_m[i, j + 2, k] = 100310021
                test_m[i, j + 4, k] = 100340022
                test_m[i, j + 6, k] = 40033
                return True
            else:
                return False
        else:
            return False
        # else:
        #     raise Exception("Values outside matrix")

    x = list(range(len(test_m)))
    y = list(range(len(test_m[0])))
    z = list(range(len(test_m[0][0])))

    # guarantee that the number of particles is specified but warns that the density of particle is too high
    particle_count = 0
    trial_max_times = 100
    trial_count = 0
    while particle_count < n:
        trial_max_times += 1
        x_random = random.choice(x)
        y_random = random.choice(y)
        z_random = random.choice(z)
        if test_m[x_random, y_random, z_random] == 0:
            if particle(x_random, y_random, z_random) == True:
                particle_count += 1
        if trial_count > trial_max_times:
            raise Exception("Could not generate required number of particles, try again or decrease density")
    # print(DICT_VECTOR_TO_INT)
    # print(DICT_VECTOR_TO_INT[(0, -2, 0)])
    #
    # test_m[3, 3, 3] = 12001
    # test_m[5, 3, 3] = 1
    return test_m
