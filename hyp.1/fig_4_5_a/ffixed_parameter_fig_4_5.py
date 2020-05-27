#This file contains the value of the parameters that I don't varry for making figure 4 and 5

from fixed_parameter import *
import numpy as np

rho = 0 #recombination rate
mig = 0.1 # migration rate
c_r = 0.1 #relative cost of choosiness
Pmut = 0.01 #frequency of introduction of the mutant
delta_c = 0 #

N_array = [N1, N2]
P = np.ones((3, 2)) / 3
Pmate = [1]
nb_mate = len(Pmate)

f = np.zeros((3, 3, nb_mate, nb_mate, 2))
for i in range(3):
    for j in range(3):
        for k in range(nb_mate):
            for l in range(nb_mate):
                for m in range(2):
                    f[i][j][k][l][m] = P[i][m] * P[j][m] * Pmate[k] * Pmate[l]

N0 = np.zeros((3, 3, nb_mate, nb_mate, 2))
for i in range(3):
    for j in range(3):
        for k in range(nb_mate):
            for l in range(nb_mate):
                for m in range(2):
                    N0[i][j][k][l][m] = (N_array[m] * f[i][j][k][l][m])

step_time = 1
