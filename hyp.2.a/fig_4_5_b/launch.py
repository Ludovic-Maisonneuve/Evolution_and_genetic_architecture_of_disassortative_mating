from fig_4_5_b.fixed_parameter_fig_4_5 import *
from functions import *
import os

L_delta = np.linspace(0,1,11)

for i_delta_a, delta_a in enumerate(L_delta):
    delta_b = delta_a

    print('delta_a', delta_a, flush=True)

    if os.path.isfile('results/Nf_' + str(i_delta_a) + '.npy') == False:
        N = np.copy(N0)

        var = 1

        while var > 0.00001:
            dN = f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)
            var = norm(dN) / np.sum(N)
            N += step_time * dN

        N_eq = np.copy(N)

        N_inv = np.zeros((3, 3, 4, 4, 2))
        P = [Pmut / 3, Pmut / 3, Pmut / 3, 1 - Pmut]
        for i in range(3):
            for j in range(3):
                for p in range(2):
                    A = N_eq[i][j][0][0][p]
                    for k in range(4):
                        for l in range(4):
                            N_inv[i][j][k][l][p] = A * P[k] * P[l]

        N_ = np.copy(N_inv)
        np.save('results/Neq_' + str(i_delta_a), N_)

        var = 1
        t = 0
        while var > 0.00001:
            t += 1
            print('t =', t, flush=True)
            print(var, flush = True)
            dN = f_with_mutant(N_, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta, c_r)
            var = norm(dN) / np.sum(N)

            N_ += step_time * dN

        np.save('results/Nf_' + str(i_delta_a), N_)
