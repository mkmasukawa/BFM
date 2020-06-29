import numpy as np
import bfm

#do test during the cooling of the soluton after athermal equilibration and test different cooling rates

#TODO: distribution of bond length

#TODO: distribution of angle length

#TODO: radius of gyration distribution

#TODO: end-to-end distance of polymer distribution

#returns the number of bonds
def recursive_find_branches(bond_list, point):


    # def in_list(p):
    #     for i in bond_list:
    #         if point in i:
    #             return True
    #     return False

    all_bonds =[]

    for b in bond_list:
        if point in b:
            if point == b[0]:
                all_bonds.append(b[1])
            if point == b[1]:
                all_bonds.append(b[0])
            bond_list.remove((b[1], b[0]))
            bond_list.remove((b[0], b[1]))

    if all_bonds:
            length = 0
            for a_b in all_bonds:
                length += recursive_find_branches(bond_list, a_b) + 1
            return length
    else:
        return 0


#TODO: size distribution of polymers
#considers all bonds are inside the matrix, if not, returns the
def size_distribution(bond_matrix):
    #get collection of bonds and remove the ones that are connected until all bonds are considered
    non_zero_indexes = np.transpose(np.nonzero(bond_matrix))
    list_of_bonds = []
    # list_of_bonds.append()
    for n_z in non_zero_indexes:
        dict_mon = bfm.split_info(bond_matrix[n_z[0], n_z[1], n_z[2]], n_z[0], n_z[1], n_z[2])
        del dict_mon["type"]
        for d_c in dict_mon:
            point1 = tuple(n_z)
            point2 = tuple(n_z + d_c)
            # if ((point1, point2) not in list_of_bonds) and ((point2, point1) not in list_of_bonds):
            if ((point1, point2) not in list_of_bonds):
                #all the bonds are duplicate, makes it easier to search for bonds
                list_of_bonds.append((point1, point2))
    #eliminate all the ones that have size 1, then 2 and so on
    # for bond in list_of_bonds:
    print(len(list_of_bonds))
    #collects the size of polymers
    size_of_polymer = []
    while len(list_of_bonds) != 0:
        monomer1 = list_of_bonds[0][0]
        size_of_polymer.append(recursive_find_branches(list_of_bonds, monomer1))

    print(size_of_polymer)




#TODO: Kratky plot

#TODO: measure diffusion coefficients g1, g2, g3

#TODO: measure entropy by checking how many available movements there are and calculate the entropy density

#TODO: scattering function

#TODO:  spin autocorrelation function (check which monomers can move) with temeprature change

#TODO: local viscosity based on the self diffusion coefficient


