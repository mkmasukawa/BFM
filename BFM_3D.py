import numpy as np
import bfm
import bfm_characterization as bc
import h5py
import time
from mayavi import mlab
import time
# generates all possible bonds

#matrix monomer 0 == empty, 1 == A, 2 == T, 3 == C, 4 == G
#covalent bond == 1, hydrogen == 2

run = 3

# name_file_initial_matrix = "/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/" + time.strftime(
#     "%Y%m%d-%H%M%S") + "small_initial_config"
# name_file_equilibrated_matrix = "/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/" + time.strftime(
#     "%Y%m%d-%H%M%S") + "small_equilibrated_config"

name_file_initial_matrix = "/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/" + str(run) + "initial_config"
name_file_equilibrated_matrix = "/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/" + str(run) + "equilibrated_config"


def generates_equilibrated_config():



    box_dimension = 128
    particle_number = 150
    equilibration_steps = 1000


    #creates the box and places the particles in ordered position
    molecules_matrix = bfm.creates_box(box_dimension, particle_number)
    with h5py.File(name_file_initial_matrix, 'w') as hf:
        hf.create_dataset("initial_config", data=molecules_matrix)
    #saves the system state




    # bfm.animates(molecules_matrix, equilibration_steps)

    for i in range(equilibration_steps):
        print(i)
        bfm.update_system(molecules_matrix)
    with h5py.File(name_file_equilibrated_matrix, 'w') as hf:
        hf.create_dataset("equilibrated_config", data=molecules_matrix)
    bfm.plots_bonds(molecules_matrix)




if __name__ == '__main__':

    MCS_steps = 100000

    # generates_equilibrated_config()

    #small one (3 molecules)
    # f = h5py.File("/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/1small_equilibrated_config", "r")
    # molecules_matrix = f.get("equilibrated_config").value

    #large configuration, actual concentration
    f = h5py.File("/Users/marcosmasukawa/Documents/BFM_Simulation/empirical/2_pipeline/3equilibrated_config", "r")
    molecules_matrix = f.get("equilibrated_config").value
    # molecules_matrix = f.get["equilibrated_value"]

    # bfm.plots_bonds(molecules_matrix)
    # bc.size_distribution(molecules_matrix)


    size = []

    start = time.time()
    with open("size_dist.txt", "a") as f:
        #when only does simulation and shoes results
        for i in range(MCS_steps):
            print(i)
            molecules_matrix = bfm.update_system_with_bonds(molecules_matrix)
            #appends every thousand interactions
            if i % 1000 == 0:
                f.write(str(bc.size_distribution(molecules_matrix)))
                bfm.plots_bonds(molecules_matrix)
                # size.append(bc.size_distribution(molecules_matrix))


    print("time took", time.time() - start)



    bfm.plots_bonds(molecules_matrix)

    # with open("size_dist.txt", "r") as f:
    #     for line in f:
    #         print(line)
            # score.append(int(line.strip()))

    # when show animation in real time
    # bfm.animates(molecules_matrix, MCS_steps)


    #TODO: make video of the time evolution and then stitch it together so it can be plotted
    #depending on the monomer density, simplify each cell into only 4 x 4 x 4 units

