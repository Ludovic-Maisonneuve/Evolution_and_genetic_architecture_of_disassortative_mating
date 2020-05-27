import os
from functions import *


L = [0] * 11
L_eq = [0] * 11
L_delta = np.linspace(0,1,11)


thisdir = 'results'
for r_, d_, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        f = file.split("_")
        f[-1] = f[-1][:-4]
        if f[0] == 'Neq':
            N = np.load(thisdir + '/' + file)
            L_eq[int(f[1])] += N
        if f[0] == 'Nf':
            N = np.load(thisdir + '/' + file)
            L[int(f[1])] += N

Result_haplo = 'delta_a = delta_b \t a-Ma	a-Mb	a-Mc	a-Mr	b-Ma	b-Mb	b-Mc	b-Mr	c-Ma	c-Mb	c-Mc	c-Mr'
Result_behavior = 'delta_a = delta_b \t self-a \t self-r'

for i, N in enumerate(L):
    PaA_1, PaB_1, PaC_1, Pa0_1, PbA_1, PbB_1, PbC_1, Pb0_1, PcA_1, PcB_1, PcC_1, Pc0_1, PaA_2, PaB_2, PaC_2, Pa0_2, PbA_2, PbB_2, PbC_2, Pb0_2, PcA_2, PcB_2, PcC_2, Pc0_2, PaA, PbA, PcA, PaB, PbB, PcB, PaC, PbC, PcC, Pa0, Pb0, Pc0  = get_haplotype(N)
    Result_haplo += '\n' + str(L_delta[i]) + '\t' + str(PaA_1) + '\t' + str(PaB_1)+ '\t' + str(PaC_1)+ '\t' + str(Pa0_1)+ '\t' + str(PbA_1)+ '\t' + str(PbB_1)+ '\t' + str(PbC_1)+ '\t' + str(Pb0_1)+ '\t' + str(PcA_1)+ '\t' + str(PcB_1)+ '\t' + str(PcC_1)+ '\t' + str(Pc0_1)
    Result_haplo += '\n' + str(L_delta[i]) + '\t' + str(PaA_2) + '\t' + str(PaB_2)+ '\t' + str(PaC_2)+ '\t' + str(Pa0_2)+ '\t' + str(PbA_2)+ '\t' + str(PbB_2)+ '\t' + str(PbC_2)+ '\t' + str(Pb0_2)+ '\t' + str(PcA_2)+ '\t' + str(PcB_2)+ '\t' + str(PcC_2)+ '\t' + str(Pc0_2) + '\n'
    P_a, P_s = get_behavior(N)
    Result_behavior += '\n' + str(L_delta[i]) + '\t' + str(P_a) + '\t' + str(P_s)

file = open('haplo.txt','a')
file.write(Result_haplo)
file.close()

file = open('behavior.txt','a')
file.write(Result_behavior)
file.close()
