from functions import *
from fixed_parameter import *
from fixed_parameter_fig_S8 import *

delta_a = 0.5 #genetic burden associated with allele a
delta_b = 0.5 #genetic burden associated with allele b
delta_c = 0 #genetic burden associated with allele c

step_time = 0.1

L_1 = []

print('Reaching equilibrium')

N = np.copy(N0)

for t in range(int(neq/step_time)):
    if (t + 1) % 2500 == 0:
        print((t + 1)*step_time, ' on ', neq, flush=True)
    N = N + step_time * f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta)

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

N = np.copy(N_inv)

for t in range(int(ninv/step_time)):
    if (t + 1) % 25 == 0:
        print(t + 1,' on ',ninv, flush=True)
    N += step_time * f_with_mutant(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta, 'disassortative', c_r)
    L_1.append(get_P0_P1(N)[1])

np.save('results/a_'+str(step_time),L_1)
