import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb

# get histogram of densities from isochrones:
spec_rho = np.genfromtxt('../data/rhostarsamples_spec.txt')

# get histogram of densities from TTVs:
ttv_rho = np.genfromtxt('../data/rhostarsamples_lightcurve.txt', invalid_raise=False).flatten()

# plot:
c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
fig = plt.figure()
ax = plt.gca()
ax.hist(spec_rho, 20, normed=True, color=c3, label='Spectroscopic')
ax.hist(ttv_rho, 20, normed=True, color=c2, label='Lightcurve')


ax.set_xlabel(r'Stellar Density (g cm$^{-3}$)',size=28)
ax.legend(loc='upper left',prop={'size':24})


plt.savefig('density.pdf')