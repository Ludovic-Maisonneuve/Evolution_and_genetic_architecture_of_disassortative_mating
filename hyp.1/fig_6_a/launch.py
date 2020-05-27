from fixed_parameter_fig_6 import *
from functions import *
import os

L = np.linspace(0,0.5,6)

for i_rho, rho in enumerate(L):

    print('rho',rho,flush=True)

    if os.path.isfile('results/Nf_' + str(i_rho) + '.npy') == False:

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
        np.save('results/Neq_' + str(i_rho), N_)

        var = 1
        t = 0
        while var > 0.00001:
            t += 1
            if t % 100 == 0:
                print('t =', t,flush=True)
                print(var,flush=True)
            dN = f_with_mutant(N_, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta,
                               'disassortative', c_r)
            var = norm(dN) / np.sum(N)

            N_ += step_time * dN

        np.save('results/Nf_' + str(i_rho), N_)
