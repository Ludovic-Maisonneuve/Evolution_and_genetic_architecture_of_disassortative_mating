from fixed_parameter_fig_2 import *
from functions import *
import os

c_r = 0  # relative cost of choosiness
L_delta = np.linspace(0, 1, 26)

Results = np.zeros((len(L_delta), len(L_delta)))

for i_delta_a, delta_a in enumerate(L_delta):
    delta_b = delta_a
    for i_delta_c, delta_c in enumerate(L_delta):

        if i_delta_a % 1 == 0 and i_delta_c % 1 == 0 and os.path.isfile('c_r_0/Nf_' + str(i_delta_a) + '_' + str(i_delta_c) + '.npy') == False:

            print('delta_a = delta_b', delta_a, 'delta_c', delta_c, flush=True)

            N = np.copy(N0)

            var = 1

            while var > 0.00001:
                dN = f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)
                var = norm(dN) / np.sum(N)
                N += step_time * dN

            N_eq = np.copy(N)

            N_inv = np.zeros((3, 3, 2, 2, 2))
            P = [1 - Pmut, Pmut]
            for i in range(3):
                for j in range(3):
                    for p in range(2):
                        A = N_eq[i][j][0][0][p]
                        for k in range(2):
                            for l in range(2):
                                N_inv[i][j][k][l][p] = A * P[k] * P[l]

            N_ = np.copy(N_inv)
            np.save('c_r_0/Neq_' + str(i_delta_a) + '_' + str(i_delta_c), N_)

            for t in range(int(ninv / step_time)):
                if (t + 1) % 500 == 0:
                    print((t + 1) * step_time, ' on ', ninv, flush=True)
                N_ += step_time * f_with_mutant(N_, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta, 'disassortative',c_r)

            np.save('c_r_0/Nf_' + str(i_delta_a) + '_' + str(i_delta_c), N_)
