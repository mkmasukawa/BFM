import numpy as np


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

vector_list_no_duplicate = np.array(vector_list_no_duplicate)