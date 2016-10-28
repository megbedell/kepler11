import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
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
ax.hist(spec_rho, normed=True, color=c3, alpha=0.5, histtype='stepfilled', bins=np.arange(0.8,1.8,0.05), label='Spectroscopic')
ax.hist(spec_rho, normed=True, color=c3, lw=3, histtype='step', bins=np.arange(0.8,1.8,0.05))
#plt.setp(patches, 'facecolor', c3, 'alpha', 0.5)
ax.hist(ttv_rho, normed=True, color=c2, alpha=0.5, histtype='stepfilled', bins=np.arange(0.8,1.8,0.05), label='Lightcurve')
ax.hist(ttv_rho, normed=True, color=c2, lw=3, histtype='step', bins=np.arange(0.8,1.8,0.05))


ax.set_xlabel(r'Stellar Density (g cm$^{-3}$)',size=28)

#spec_patch = mpatches.Patch(color=c3, alpha=0.5, label='Spectroscopic')
#ttv_patch = mpatches.Patch(color=c2, alpha=0.5, label='Lightcurve')
#ax.legend(handles=[spec_patch, ttv_patch],loc='upper left',prop={'size':24})

#plt.text(1.122,0.9,'L13 result',ha='center',size=24)
plt.errorbar([1.122],[4.0], xerr=[[0.060],[0.049]], lw=3, color=c1, capsize=5, capthick=3, fmt='o', markersize=9, label='Lissauer+2013')

ax.legend(loc='upper left',prop={'size':24}, numpoints=1)
plt.xlim([0.7,1.8])

plt.savefig('density.pdf')