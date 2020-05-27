from fixed_parameter_fig_S10 import *
from functions import *
import os

L_delta = np.linspace(0,1,11)
for i_delta_a, delta_a in enumerate(L_delta):
    delta_b = delta_a

    print('delta_a', delta_a, flush=True)

    if os.path.isfile('results/Nf_' + str(i_delta_a) + '.npy') == False:

        N = np.copy(N0)

        for t in range(int(neq / step_time)):
            if (t + 1) % 2500 == 0:
                print((t + 1) * step_time, ' on ', neq, flush=True)
            N += step_time * f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)

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

        for t in range(int(ninv / step_time)):
            if (t + 1) % 1 == 0:
                print((t+ 1) * step_time, ' on ', ninv, flush=True)
            N_ += step_time * f_with_mutant(N_, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta, c_r)
            print(get_N_A_B_C_(N_))

        np.save('results/Nf_' + str(i_delta_a), N_)
