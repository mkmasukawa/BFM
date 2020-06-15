import numpy as np
import bfm
# generates all possible bonds

#matrix monomer 0 == empty, 1 == A, 2 == T, 3 == C, 4 == G
#covalent bond == 1, hydrogen == 2

if __name__ == '__main__':

    box_dimension = 10
    particle_number = 1
    #creates the box and places the particles in ordered position
    molecules_matrix = bfm.creates_box(box_dimension, particle_number)

    bfm.plots_bonds(molecules_matrix)
