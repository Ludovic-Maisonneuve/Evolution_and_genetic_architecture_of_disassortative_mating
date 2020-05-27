import numpy as np
import os
import matplotlib.pyplot as plt
from functions import *
import matplotlib.colors as colors

SMALL_SIZE = 8
MEDIUM_SIZE = 12
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

L_delta = np.linspace(0,1,26)


R_0_1 = np.zeros((16, 16))
thisdir = 'c_r_0_1_'
for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        f = file.split("_")
        f[-1] = f[-1][:-4]
        if f[0] == 'Nf':
            N = np.load(thisdir + '/' + file)
            R_0_1[int(f[1]), int(f[2])] += get_P0_P1(N)[1]


P_max = np.max(R_0_1)

plt.figure(figsize=(6.2,5))
plt.imshow(R_0_1.T,extent=[0, 1, 0, 1],aspect=1,origin='lower',vmin=0, vmax=P_max,interpolation='None',cmap='RdBu',norm=MidpointNormalize(midpoint=0.01,vmin=0, vmax=P_max))
plt.colorbar(label=r'frequency of allele $\it{dis}$')
plt.xlabel(r'$\delta_a=\delta_b$')
plt.ylabel(r'$\delta_c$')
plt.title(r'$c_r = 0.1$')

plt.show()
