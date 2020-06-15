import numpy as np
import itertools
from mayavi import mlab

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
                    connections.append([(i, j, k), (i+bond[1], j+bond[2], k+bond[3])])
                    # x.append((i, j, k))
                    # x.append((i+bond[1], j+bond[2], k+bond[3]))
                    # x.append((i+bond[1], j+bond[2], k+bond[3]))
                    # connections.append([i*x_len + x_len + j*y_len + k*z_len, (i+bond[1])*x_len + (j+bond[2])*y_len + (k+bond[3])*z_len])


    connection_index = []
    for c in connections:
        connection_index.append([x.index(c[0]), x.index(c[1])])

    u, v, w = np.array(x).T
    scalars = np.array(scalars)

    # connections = [(1, 2), (3, 4)]

    # for start, stop in edges:
    #     connections.append((labels[start], labels[stop]))
    # mlab.figure(1, bgcolor=(0, 0, 0))
    mlab.clf()

    pts = mlab.points3d(u*10, v*10, w*10, 1.5 * scalars.max() - scalars,
                        scale_factor=10, resolution=7)
    pts.mlab_source.dataset.lines = np.array(connection_index)

    # mlab.plot3d([1,2,3, 1], [2,3,4, 2], [5,6,7, 1], [1, 2, 1, 3],
    #             tube_radius=0.1, colormap='Reds', tube_sides = 6)

    #Use a tube fiter to plot tubes on the link, varying the radius with the
   # scalar value

    tube = mlab.pipeline.tube(pts, tube_radius=0.0001, tube_sides=4)
    tube.filter.radius_factor = 1
    tube.filter.vary_radius = 'vary_radius_by_scalar'
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0))

    # Visualize the local atomic density
    # mlab.pipeline.volume(mlab.pipeline.gaussian_splatter(pts))


    mlab.show()

    #obj = mlab.quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
    #mlab.show()
    # return obj




def creates_box(M, N):
    test_m = np.zeros((M, M, M), dtype=int)
    test_m[3, 3, 3] = 12001
    test_m[5, 3, 3] = 1
    return test_m

    # for i in range(M):
    #     test_m[i]