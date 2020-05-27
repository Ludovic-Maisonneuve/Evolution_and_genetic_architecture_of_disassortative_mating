import numpy as np

def indic(a, b):  # indicatrice function of a = b
    if a == b:
        return 1.0
    else:
        return 0.0


def coef_haplotype(P_offspring, M_offspring, P_mother_parent, P_father_parent, M_mother_parent, M_father_parent, rho):
    return (1 - rho) / 2 * indic(P_offspring, P_mother_parent) * indic(M_offspring, M_mother_parent) \
           + (1 - rho) / 2 * indic(P_offspring, P_father_parent) * indic(M_offspring,M_father_parent) \
           + rho / 2 * indic(P_offspring, P_mother_parent) * indic(M_offspring, M_father_parent) \
           + rho / 2 * indic(P_offspring, P_father_parent) * indic(M_offspring, M_mother_parent)


def coef(P_mother, P_father, M_mother, M_father, P_mother_female, P_father_female, M_mother_female, M_father_female,
         P_mother_male, P_father_male, M_mother_male, M_father_male, rho):
    return coef_haplotype(P_mother, M_mother, P_mother_female, P_father_female, M_mother_female, M_father_female,
                          rho) * coef_haplotype(P_father, M_father, P_mother_male, P_father_male, M_mother_male,
                                                M_father_male, rho)


def genetic_burden(a1, a2, delta_a, delta_b, delta_c, delta):
    if a1 == 0 and a2 == 0:
        #return (1 - delta1) * (1 - deltah)
        return 1 - delta_a - delta
    elif a1 == 1 and a2 == 1:
        #return (1 - delta2) * (1 - deltah)
        return 1 - delta_b - delta
    elif a1 == 2 and a2 == 2:
        #return (1 - delta3) * (1 - deltah)
        return 1 - delta_c - delta
    else:
        return (1 - delta)


def Reproduction_random(N, r, K, rho):
    nb_mate = np.shape(N)[2]
    N1, N2 = get_N1_N2(N)

    Res = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i == j:
                Res[i][j] = 1

    f = np.zeros((3, 3, nb_mate, nb_mate, 2))
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    f[i][j][k][l][0] = N[i][j][k][l][0] / N1
                    f[i][j][k][l][1] = N[i][j][k][l][1] / N2

    Rep = np.zeros((3, 3, nb_mate, nb_mate, 2))

    f_ = np.zeros((3, 3, nb_mate, nb_mate, 2))
    f_1 = 0
    f_2 = 0
    for a in range(3):
        for b in range(3):
            for c in range(nb_mate):
                for d in range(nb_mate):
                    for e in range(2):
                        for i in range(3):
                            for j in range(3):
                                for k in range(nb_mate):
                                    for l in range(nb_mate):
                                        for m in range(3):
                                            for n in range(3):
                                                for o in range(nb_mate):
                                                    for p in range(nb_mate):
                                                        f_[a][b][c][d][e] += \
                                                            coef(a, b, c, d, i, j, k, l, m, n, o, p,
                                                                 rho) * \
                                                            f[i][j][k][l][e] * f[m][n][o][p][e]
                        if e == 0:
                            f_1 += f_[a][b][c][d][e]
                        else:
                            f_2 += f_[a][b][c][d][e]

    for a in range(3):
        for b in range(3):
            for c in range(nb_mate):
                for d in range(nb_mate):
                    f_[a][b][c][d][0] = f_[a][b][c][d][0] / f_1
                    f_[a][b][c][d][1] = f_[a][b][c][d][1] / f_2

    for m in range(2):
        Ntot = 0
        for i_ in range(3):
            for j_ in range(3):
                for k_ in range(nb_mate):
                    for l_ in range(nb_mate):
                        Ntot = Ntot + N[i_][j_][k_][l_][m]
        for i in range(3):
            for j in range(3):
                for k in range(nb_mate):
                    for l in range(nb_mate):
                        Rep[i][j][k][l][m] = r * (1 - Ntot / K) * Ntot * f_[i][j][k][l][m]

    return Rep


def Reproduction_with_mutant(N, r, K, rho, mut, c_r):
    nb_mate = np.shape(N)[2]
    N1, N2 = get_N1_N2(N)
    PA1, PB1, PC1, PA2, PB2, PC2 = get_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N)
    Res = np.eye(3)

    f = np.zeros((3, 3, nb_mate, nb_mate, 2))
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    f[i][j][k][l][0] = N[i][j][k][l][0] / N1
                    f[i][j][k][l][1] = N[i][j][k][l][1] / N2

    Rep = np.zeros((3, 3, nb_mate, nb_mate, 2))

    Pref = np.zeros((3, 3, 2, 2, 3))
    if mut == 'disassortative':
        for p1 in range(3):
            for p2 in range(3):
                for m1 in range(2):
                    for m2 in range(2):
                        for morph in range(3):
                            if m1 == m2 and m1 == 0:
                                Pref[p1][p2][m1][m2][morph] = 1
                            else:
                                Pref[p1][p2][m1][m2][morph] = (1 - Res[min(p1, p2)][morph])

    elif mut == 'assortative':
        for p1 in range(3):
            for p2 in range(3):
                for m1 in range(2):
                    for m2 in range(2):
                        for morph in range(3):
                            if m1 == m2 and m1 == 0:
                                Pref[p1][p2][m1][m2][morph] = 1
                            else:
                                Pref[p1][p2][m1][m2][morph] = Res[min(p1, p2)][morph]

    f_ = np.zeros((3, 3, nb_mate, nb_mate, 2))
    f_1 = 0
    f_2 = 0
    for a in range(3):
        for b in range(3):
            for c in range(nb_mate):
                for d in range(nb_mate):
                    for e in range(2):
                        for i in range(3):
                            for j in range(3):
                                for k in range(nb_mate):
                                    for l in range(nb_mate):
                                        if e == 0:
                                            F = Pref[i][j][k][l][0] * PA1 + Pref[i][j][k][l][1] * PB1 + \
                                                Pref[i][j][k][l][2] * PC1
                                        else:
                                            F = Pref[i][j][k][l][0] * PA2 + Pref[i][j][k][l][1] * PB2 + \
                                                Pref[i][j][k][l][2] * PC2
                                        if F != 0:
                                            M = (1 - c_r + c_r * F) / F
                                        else:
                                            M = 0
                                        for m in range(3):
                                            for n in range(3):
                                                for o in range(nb_mate):
                                                    for p in range(nb_mate):
                                                        f_[a][b][c][d][e] += M * (Pref[i][j][k][l][min(m, n)]) * coef(a,
                                                                                                                      b,
                                                                                                                      c,
                                                                                                                      d,
                                                                                                                      i,
                                                                                                                      j,
                                                                                                                      k,
                                                                                                                      l,
                                                                                                                      m,
                                                                                                                      n,
                                                                                                                      o,
                                                                                                                      p,
                                                                                                                      rho) * \
                                                                             f[i][j][k][l][e] * f[m][n][o][p][e]
                        if e == 0:
                            f_1 += f_[a][b][c][d][e]
                        else:
                            f_2 += f_[a][b][c][d][e]

    for a in range(3):
        for b in range(3):
            for c in range(nb_mate):
                for d in range(nb_mate):
                    f_[a][b][c][d][0] = f_[a][b][c][d][0] / f_1
                    f_[a][b][c][d][1] = f_[a][b][c][d][1] / f_2
    for m in range(2):
        Ntot = 0
        for i_ in range(3):
            for j_ in range(3):
                for k_ in range(nb_mate):
                    for l_ in range(nb_mate):
                        Ntot = Ntot + N[i_][j_][k_][l_][m]
        for i in range(3):
            for j in range(3):
                for k in range(nb_mate):
                    for l in range(nb_mate):
                        Rep[i][j][k][l][m] = r * (1 - Ntot / K) * Ntot * f_[i][j][k][l][m]

    return Rep


def Predation(N, d_m, d_nm, lmbda):
    nb_mate = np.shape(N)[2]

    Res = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            if i == j:
                Res[i][j] = 1

    Pred = np.zeros((3, 3, nb_mate, nb_mate,
                     2))  # Create the matrix of the predation at the same dimention as N, each element is the dP on the N element corresponding
    for i in range(3):  # every possible value of the first allele of the first loci
        for j in range(3):  # every possible value of the second allele of the first loci
            for k in range(nb_mate):  # every possible value of the first allele of the second loci
                for l in range(nb_mate):  # every possible value of the second allele of the second loci
                    for m in range(2):  # the number of the population
                        sum = 0
                        for i_ in range(3):
                            for j_ in range(3):
                                for k_ in range(nb_mate):
                                    for l_ in range(nb_mate):
                                        sum = sum + lmbda * Res[min(i, j)][min(i_, j_)] * N[i_][j_][k_][l_][m]
                        Pred[i][j][k][l][m] = - (
                                 d_nm * (1-Res[min(i, j)][m]) + d_m * Res[min(i, j)][
                            m]) * N[i][j][k][l][m] / (1 + sum)
    return Pred


def Migration(N, mig):
    nb_mate = np.shape(N)[2]

    Mig = np.zeros((3, 3, nb_mate, nb_mate, 2))
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    for m in range(2):
                        Mig[i][j][k][l][m] = mig * (N[i][j][k][l][-m + 1] - N[i][j][k][l][m])
        return Mig


def Survival(N, delta_a, delta_b, delta_c, delta):
    nb_mate = np.shape(N)[2]

    Suv = np.zeros((3, 3, nb_mate, nb_mate, 2))
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    for m in range(2):
                        Suv[i][j][k][l][m] = - (1 - genetic_burden(i, j, delta_a, delta_b, delta_c, delta)) * \
                                             N[i][j][k][l][
                                                 m]
    return Suv


def f_random(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta):
    dN = Reproduction_random(N, r, K, rho) + Migration(N, mig) + Predation(N, d_m, d_nm, lmbda) + Survival(N, delta_a,
                                                                                                           delta_b,
                                                                                                           delta_c, delta)

    return dN


def f_with_mutant(N, r, K, rho, d_m, d_nm, lmbda, mig, delta_a, delta_b, delta_c, delta, mut, c_r):
    R = Reproduction_with_mutant(N, r, K, rho, mut, c_r)
    M = Migration(N, mig)
    P = Predation(N, d_m, d_nm, lmbda)
    S = Survival(N, delta_a, delta_b, delta_c, delta)
    dN = R + M + P + S
    return dN


def get_N1_N2(N):
    nb_mate = np.shape(N)[2]

    N1 = 0
    N2 = 0
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    N1 += N[i][j][k][l][0]
                    N2 += N[i][j][k][l][1]
    return N1, N2


def get_Pa_Pb_Pc_Pa_1_Pb_1_Pc_1_Pa_2_Pb_2_Pc_2_PA_PB_PC_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N):
    nb_mate = np.shape(N)[2]
    N1, N2 = get_N1_N2(N)
    Ntot = np.sum(N)

    Na1 = 0
    Nb1 = 0
    Nc1 = 0

    Na2 = 0
    Nb2 = 0
    Nc2 = 0

    NA1 = 0
    NB1 = 0
    NC1 = 0

    NA2 = 0
    NB2 = 0
    NC2 = 0

    for i in range(nb_mate):
        for j in range(nb_mate):
            Na1 += N[0][0][i][j][0] + 0.5 * (
                    N[0][1][i][j][0] + N[1][0][i][j][0] + N[0][2][i][j][0] + N[2][0][i][j][0])
            Nb1 += N[1][1][i][j][0] + 0.5 * (
                    N[0][1][i][j][0] + N[1][0][i][j][0] + N[1][2][i][j][0] + N[2][1][i][j][0])
            Nc1 += N[2][2][i][j][0] + 0.5 * (
                    N[0][2][i][j][0] + N[2][0][i][j][0] + N[1][2][i][j][0] + N[2][1][i][j][0])

            Na2 += N[0][0][i][j][1] + 0.5 * (
                    N[0][1][i][j][1] + N[1][0][i][j][1] + N[0][2][i][j][1] + N[2][0][i][j][1])
            Nb2 += N[1][1][i][j][1] + 0.5 * (
                    N[0][1][i][j][1] + N[1][0][i][j][1] + N[1][2][i][j][1] + N[2][1][i][j][1])
            Nc2 += N[2][2][i][j][1] + 0.5 * (
                    N[0][2][i][j][1] + N[2][0][i][j][1] + N[1][2][i][j][1] + N[2][1][i][j][1])

            NA1 += N[0][0][i][j][0] + N[0][1][i][j][0] + N[1][0][i][j][0] + N[0][2][i][j][0] + N[2][0][i][j][0]
            NB1 += N[1][1][i][j][0] + N[1][2][i][j][0] + N[2][1][i][j][0]
            NC1 += N[2][2][i][j][0]

            NA2 += N[0][0][i][j][1] + N[0][1][i][j][1] + N[1][0][i][j][1] + N[0][2][i][j][1] + N[2][0][i][j][1]
            NB2 += N[1][1][i][j][1] + N[1][2][i][j][1] + N[2][1][i][j][1]
            NC2 += N[2][2][i][j][1]

    Pa = (Na1 + Na2) / Ntot
    Pb = (Nb1 + Nb2) / Ntot
    Pc = (Nc1 + Nc2) / Ntot

    Pa1 = Na1 / N1
    Pb1 = Nb1 / N1
    Pc1 = Nc1 / N1

    Pa2 = Na2 / N2
    Pb2 = Nb2 / N2
    Pc2 = Nc2 / N2

    PA = (NA1 + NA2) / Ntot
    PB = (NB1 + NB2) / Ntot
    PC = (NC1 + NC2) / Ntot

    PA1 = NA1 / N1
    PB1 = NB1 / N1
    PC1 = NC1 / N1

    PA2 = NA2 / N2
    PB2 = NB2 / N2
    PC2 = NC2 / N2

    return Pa, Pb, Pc, Pa1, Pb1, Pc1, Pa2, Pb2, Pc2, PA, PB, PC, PA1, PB1, PC1, PA2, PB2, PC2


def get_P0_P1(N):
    nb_mate = np.shape(N)[2]
    Ntot = np.sum(N)

    if nb_mate == 1:
        return 1, 0

    N0 = 0
    N1 = 0

    for i in range(3):
        for j in range(3):
            for p in range(2):
                N0 += N[i][j][0][0][p] + 0.5 * (N[i][j][0][1][p] + N[i][j][1][0][p])
                N1 += N[i][j][1][1][p] + 0.5 * (N[i][j][0][1][p] + N[i][j][1][0][p])

    P0 = N0 / Ntot
    P1 = N1 / Ntot

    return P0, P1


def get_Pa0_1_Pb0_1_Pc0_1_Pa1_1_Pb1_1_Pc1_1_Pa0_2_Pb0_2_Pc0_2_Pa1_2_Pb1_2_Pc1_2_Pa0_Pb0_Pc0_Pa1_Pb1_Pc1(N):
    nb_mate = np.shape(N)[2]
    N1, N2 = get_N1_N2(N)
    Ntot = np.sum(N)

    Na01 = 0
    Nb01 = 0
    Nc01 = 0
    Na11 = 0
    Nb11 = 0
    Nc11 = 0

    Na02 = 0
    Nb02 = 0
    Nc02 = 0
    Na12 = 0
    Nb12 = 0
    Nc12 = 0

    for i in range(3):
        for j in range(nb_mate):
            Na01 += N[i][0][j][0][0] + N[0][i][0][j][0]
            Nb01 += N[i][1][j][0][0] + N[1][i][0][j][0]
            Nc01 += N[i][2][j][0][0] + N[2][i][0][j][0]
            Na11 += N[i][0][j][1][0] + N[0][i][1][j][0]
            Nb11 += N[i][1][j][1][0] + N[1][i][1][j][0]
            Nc11 += N[i][2][j][1][0] + N[2][i][1][j][0]

            Na02 += N[i][0][j][0][1] + N[0][i][0][j][1]
            Nb02 += N[i][1][j][0][1] + N[1][i][0][j][1]
            Nc02 += N[i][2][j][0][1] + N[2][i][0][j][1]
            Na12 += N[i][0][j][1][1] + N[0][i][1][j][1]
            Nb12 += N[i][1][j][1][1] + N[1][i][1][j][1]
            Nc12 += N[i][2][j][1][1] + N[2][i][1][j][1]

    Pa01 = Na01 / (2 * N1)
    Pb01 = Nb01 / (2 * N1)
    Pc01 = Nc01 / (2 * N1)
    Pa11 = Na11 / (2 * N1)
    Pb11 = Nb11 / (2 * N1)
    Pc11 = Nc11 / (2 * N1)

    Pa02 = Na02 / (2 * N2)
    Pb02 = Nb02 / (2 * N2)
    Pc02 = Nc02 / (2 * N2)
    Pa12 = Na12 / (2 * N2)
    Pb12 = Nb12 / (2 * N2)
    Pc12 = Nc12 / (2 * N2)

    Pa0 = (Na01 + Na02) / (2 * Ntot)
    Pb0 = (Nb01 + Nb02) / (2 * Ntot)
    Pc0 = (Nc01 + Nc02) / (2 * Ntot)
    Pa1 = (Na11 + Na12) / (2 * Ntot)
    Pb1 = (Nb11 + Nb12) / (2 * Ntot)
    Pc1 = (Nc11 + Nc12) / (2 * Ntot)

    return Pa01, Pb01, Pc01, Pa11, Pb11, Pc11, Pa02, Pb02, Pc02, Pa12, Pb12, Pc12, Pa0, Pb0, Pc0, Pa1, Pb1, Pc1


def get_heterozygote_rate_homozygote_rate(N):
    nb_mate = np.shape(N)[2]
    Ntot = np.sum(N)

    Nhetero = 0
    Nhomo = 0
    Nhomo_i = [0, 0, 0]
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    for p in range(2):
                        if i == j:
                            Nhomo += N[i][j][k][l][p]
                            Nhomo_i[i] += N[i][j][k][l][p]

                        else:
                            Nhetero += N[i][j][k][l][p]

    Phetero = Nhetero / Ntot
    Phomo = Nhomo / Ntot
    PhomoA = Nhomo_i[0] / Ntot
    PhomoB = Nhomo_i[1] / Ntot
    PhomoC = Nhomo_i[2] / Ntot

    return Phetero, Phomo, PhomoA, PhomoB, PhomoC


def get_PA_1_PB_1_PC_1_PA_2_PB_2_PC_2(N):
    nb_mate = np.shape(N)[2]
    N1, N2 = get_N1_N2(N)

    NA1 = 0
    NB1 = 0
    NC1 = 0

    NA2 = 0
    NB2 = 0
    NC2 = 0

    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    if min(i, j) == 0:
                        NA1 += N[i][j][k][l][0]
                        NA2 += N[i][j][k][l][1]
                    elif min(i, j) == 1:
                        NB1 += N[i][j][k][l][0]
                        NB2 += N[i][j][k][l][1]
                    else:
                        NC1 += N[i][j][k][l][0]
                        NC2 += N[i][j][k][l][1]

    PA1 = NA1 / N1
    PB1 = NB1 / N1
    PC1 = NC1 / N1

    PA2 = NA2 / N2
    PB2 = NB2 / N2
    PC2 = NC2 / N2

    return PA1, PB1, PC1, PA2, PB2, PC2

def norm(N):
    nb_mate = np.shape(N)[2]
    Norm = 0
    for i in range(3):
        for j in range(3):
            for k in range(nb_mate):
                for l in range(nb_mate):
                    for p in range(2):
                        Norm += N[i][j][k][l][p]**2
    return np.sqrt(Norm)

def get_self_av(N):
    nb_mate = np.shape(N)[2]
    Ntot = np.sum(N)

    if nb_mate == 1:
        0

    N1 = 0

    for i in range(3):
        for j in range(3):
            for p in range(2):
                N1 += N[i][j][1][1][p] + N[i][j][0][1][p] + N[i][j][1][0][p]

    P1 = N1 / Ntot

    return P1
