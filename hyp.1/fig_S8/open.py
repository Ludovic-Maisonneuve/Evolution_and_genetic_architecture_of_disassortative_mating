import numpy as np
import matplotlib.pyplot as plt

a_1 = np.load('results/a_1.npy')
a_0_1 = np.load('results/a_0.1.npy')
a_0_01 = np.load('results/a_0.01.npy')
b_1 = np.load('results/b_1.npy')
b_0_1 = np.load('results/b_0.1.npy')
b_0_01 = np.load('results/b_0.01.npy')

plt.figure()
plt.plot(np.linspace(0,2500,2500), a_1, label=r'$\Delta t = 1$', alpha=0.7)
plt.plot(np.linspace(0,2500,25000), a_0_1, label=r'$\Delta t = 0.1$', alpha=0.7)
plt.plot(np.linspace(0,2500,250000), a_0_01, label=r'$\Delta t = 0.01$', alpha=0.7)

plt.legend()
plt.xlabel('t')
plt.ylabel(r'Frequency of allele $dis$')
plt.show()

plt.figure()
plt.plot(np.linspace(0,2500,2500), b_1, label=r'$\Delta t = 1$', alpha=0.7)
plt.plot(np.linspace(0,2500,25000), b_0_1, label=r'$\Delta t = 0.1$', alpha=0.7)
plt.plot(np.linspace(0,2500,250000), b_0_01, label=r'$\Delta t = 0.01$', alpha=0.7)

plt.legend()
plt.xlabel('t')
plt.ylabel(r'Frequency of allele $dis$')
plt.show()
