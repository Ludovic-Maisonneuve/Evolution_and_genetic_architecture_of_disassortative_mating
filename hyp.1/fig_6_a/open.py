import os
from functions import *

L = [0] * 6
L_rho = np.linspace(0, 1, 6)

thisdir = 'results'
for r_, d_, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        f = file.split("_")
        f[-1] = f[-1][:-4]
        #if f[0] == 'Neq':
        #    N = np.load(thisdir + '/' + file)
        #    R[int(f[1]), int(f[2])] += - get_P0_P1(N)[1]
        if f[0] == 'Nf':
            N = np.load(thisdir + '/' + file)
            L[int(f[1])] += N

Result_haplo = 'rho \t a-r \t b-r \t c-r \t a-dis \t b-dis \t c-dis'
Result_behavior = 'rho = delta_b \t self-a \t self-r'

for i, N in enumerate(L):

    if str(type(N)) == "<class 'numpy.ndarray'>":

        Pa01, Pb01, Pc01, Pa11, Pb11, Pc11, Pa02, Pb02, Pc02, Pa12, Pb12, Pc12, Pa0, Pb0, Pc0, Pa1, Pb1, Pc1 = get_Pa0_1_Pb0_1_Pc0_1_Pa1_1_Pb1_1_Pc1_1_Pa0_2_Pb0_2_Pc0_2_Pa1_2_Pb1_2_Pc1_2_Pa0_Pb0_Pc0_Pa1_Pb1_Pc1(N)
        Result_haplo += '\n' + str(L_rho[i]) + '\t' + str(Pa01) + '\t' + str(Pb01) + '\t' + str(Pc01) + '\t' + str(Pa11) + '\t' + str(Pb11) + '\t' + str(Pc11)
        Result_haplo += '\n' + str(L_rho[i]) + '\t' + str(Pa02) + '\t' + str(Pb02) + '\t' + str(Pc02) + '\t' + str(Pa12) + '\t' + str(Pb12) + '\t' + str(Pc12) + '\n'

        P_s_r = get_self_av(N)
        Result_behavior += '\n' + str(L_rho[i]) + '\t' + str(1 - P_s_r) + '\t' + str(P_s_r)

file = open('haplo.txt','a')
file.write(Result_haplo)
file.close()

file = open('behavior.txt','a')
file.write(Result_behavior)
file.close()
