import numpy as np
from mayavi import mlab
import itertools

# generates all possible bonds


# View it.
# from mayavi import mlab
# s = mlab.mesh(x, y, z)
# mlab.show()

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

def test_quiver3d(vec):
    n = len(vec)
    print(n)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    u, v, w = vec.T

    obj = mlab.quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
    return obj


if __name__ == '__main__':
    test_quiver3d(vector_list_no_duplicate)
    mlab.show()
    print(len(vector_list_no_duplicate))
