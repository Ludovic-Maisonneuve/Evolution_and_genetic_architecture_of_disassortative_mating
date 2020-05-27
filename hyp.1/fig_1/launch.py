from fixed_parameter_fig_1 import *
from functions import *
import os

Results = 'mig \t mating choice \t PA 1 \t PB 1 \t PC 1 \t PA 2 \t PB 2 \t PC 2'
#
L_mig = np.linspace(0, 0.5, 6)
#
for mig in L_mig:
    print('mig = ', mig,flush=True)

    N = np.copy(N0)
    var = 1

    while var > 0.000001:
        dN = f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)
        var = norm(dN) / np.sum(N)
        N += step_time * dN

    Pa, Pb, Pc, Pa1, Pb1, Pc1, Pa2, Pb2, Pc2, PA, PB, PC, PA1, PB1, PC1, PA2, PB2, PC2 = get_Pa_Pb_Pc_Pa_1_Pb_1_Pc_1_Pa_2_Pb_2_Pc_2_PA_PB_PC_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N)
    Results += '\n' + str(mig) + '\t Random + \t' + str(PA1) + '\t' + str(PB1) + '\t' + str(PC1) + '\t' + str(PA2) + '\t' + str(PB2) + '\t' + str(PC2)

    N_ = np.copy(N0)
    N = np.zeros((3, 3, 2, 2, 2))
    P = [0, 1]
    for i in range(3):
        for j in range(3):
            for p in range(2):
                A = N_[i][j][0][0][p]
                for k in range(2):
                    for l in range(2):
                        N[i][j][k][l][p] = A * P[k] * P[l]

    var = 1
    while var > 0.00001:
        dN = f_with_mutant(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta,
                           'disassortative', c_r)
        var = norm(dN) / np.sum(N)
        N += step_time * dN

    Pa, Pb, Pc, Pa1, Pb1, Pc1, Pa2, Pb2, Pc2, PA, PB, PC, PA1, PB1, PC1, PA2, PB2, PC2 = get_Pa_Pb_Pc_Pa_1_Pb_1_Pc_1_Pa_2_Pb_2_Pc_2_PA_PB_PC_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N)
    Results += '\n' + str(mig) + '\t Disassortative + \t' + str(PA1) + '\t' + str(PB1) + '\t' + str(PC1) + '\t' + str(
        PA2) + '\t' + str(PB2) + '\t' + str(PC2)

    N_ = np.copy(N0)
    N = np.zeros((3, 3, 2, 2, 2))
    P = [0, 1]
    for i in range(3):
        for j in range(3):
            for p in range(2):
                A = N_[i][j][0][0][p]
                for k in range(2):
                    for l in range(2):
                        N[i][j][k][l][p] = A * P[k] * P[l]

    var = 1
    while var > 0.00001:
        dN = f_with_mutant(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta,
                           'assortative', c_r)
        var = norm(dN) / np.sum(N)
        N += step_time * dN

    Pa, Pb, Pc, Pa1, Pb1, Pc1, Pa2, Pb2, Pc2, PA, PB, PC, PA1, PB1, PC1, PA2, PB2, PC2 = get_Pa_Pb_Pc_Pa_1_Pb_1_Pc_1_Pa_2_Pb_2_Pc_2_PA_PB_PC_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(
        N)
    Results += '\n' + str(mig) + '\t assortative + \t' + str(PA1) + '\t' + str(PB1) + '\t' + str(PC1) + '\t' + str(
        PA2) + '\t' + str(PB2) + '\t' + str(PC2)

file = open('results_.txt','a')
file.write(Results)
file.close()
