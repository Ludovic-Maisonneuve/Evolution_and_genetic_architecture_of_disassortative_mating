from fixed_parameter_fig_S9 import *
from functions import *

c_r = 0.1  # relative cost of choosiness
L_delta = [0,0.25,0.5,1]

Results = 'delta_a = delta_b \t delta_c \t PA1 \t PB1 \t PC1 \t PA2 \t PB2 \t PC2'

for i_delta_a, delta_a in enumerate(L_delta):
    delta_b = delta_a
    for i_delta_c, delta_c in enumerate(L_delta):

        print('delta_a = delta_b', delta_a, 'delta_c', delta_c, flush=True)

        N = np.copy(N0)

        var = 1

        while var > 0.00001:
            dN = f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)
            var = norm(dN) / np.sum(N)
            N += step_time * dN

        N_eq = np.copy(N)

        PA1, PB1, PC1, PA2, PB2, PC2 = get_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N_eq)

        Results += '\n' + str(delta_a) + '\t' + str(delta_c) + '\t' + str(PA1) + '\t' + str(PB1) + '\t' + str(
            PC1) + '\t' + str(PA2) + '\t' + str(PB2) + '\t' + str(PC2)

file = open('results.txt','a')
file.write(Results)
file.close()
