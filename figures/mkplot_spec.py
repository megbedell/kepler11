import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import matplotlib
import numpy as np
from astropy.io import fits

#import some data
sp = fits.open('/Users/mbedell/Documents/Research/Kepler11/Spectra/K11_2d.fits')
header = sp[0].header
naxis1 = header['NAXIS1']
crval1 = header['CRVAL1']
cdelt1 = header['CDELT1']
K11_wave = crval1 + np.arange(naxis1)*cdelt1
K11_flux = sp[0].data
sp = fits.open('/Users/mbedell/Documents/Research/Kepler11/Spectra/Sun_2d.fits')
header = sp[0].header
naxis1 = header['NAXIS1']
crval1 = header['CRVAL1']
cdelt1 = header['CDELT1']
Sun_wave = crval1 + np.arange(naxis1)*cdelt1
Sun_flux = sp[0].data
sp = fits.open('/Users/mbedell/Documents/Research/Kepler11/Spectra/HD1178_2d.fits')
header = sp[0].header
naxis1 = header['NAXIS1']
crval1 = header['CRVAL1']
cdelt1 = header['CDELT1']
HD1178_wave = crval1 + np.arange(naxis1)*cdelt1
HD1178_flux = sp[0].data
"""sp = fits.open('/Users/mbedell/Documents/Research/Kepler11/Spectra/HD10145_2d.fits')
header = sp[0].header
naxis1 = header['NAXIS1']
crval1 = header['CRVAL1']
cdelt1 = header['CDELT1']
HD10145_wave = crval1 + np.arange(naxis1)*cdelt1
HD10145_flux = sp[0].data"""

#set up the figure
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20

fig = plt.figure()
gs = gridspec.GridSpec(2,1,height_ratios=[4,1],hspace=0.1)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax2.ticklabel_format(useOffset=False)
majorLocator = MultipleLocator(1)
minorLocator = MultipleLocator(0.1)
ax1.xaxis.set_minor_locator(minorLocator)
ax1.xaxis.set_major_locator(majorLocator)
ax2.xaxis.set_minor_locator(minorLocator)
ax2.xaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(0.05)
ax2.yaxis.set_major_locator(majorLocator)
plt.setp(ax1.get_xticklabels(), visible=False)


c1 = 'black'
c2 = '#003399'
c3 = '#CC0033'

#plot spectra
ax1.plot(K11_wave,K11_flux+0.01,label=r'Kepler-11',color=c1)
ax1.plot(Sun_wave,Sun_flux+0.01,label=r'Sun',color=c2)
ax1.plot(HD1178_wave,HD1178_flux+0.01,label=r'HD1178',color=c3)
ax1.text(6247.55,0.6,'Fe II',ha='center',size=24)
ax1.text(6246.3,0.35,'Fe I',ha='center',size=24)
ax1.text(6245.6,0.67,'Sc II',ha='center',size=24)
ax1.text(6244.45,0.65,'Si I',ha='center',size=24)

ax1.set_xlim([6244.0,6248.0])
ax1.set_ylim([0.31,1.05])
ax1.set_ylabel(r'Normalized Flux',size=28)
ax1.legend(loc='lower left',prop={'size':24})

Sun_resids = Sun_flux - np.interp(Sun_wave, K11_wave, K11_flux)
HD1178_resids = HD1178_flux - np.interp(HD1178_wave, K11_wave, K11_flux)

ax2.plot(K11_wave,np.zeros_like(K11_wave),color=c1)
ax2.plot(Sun_wave,Sun_resids,color=c2)
ax2.plot(HD1178_wave,HD1178_resids,color=c3)
ax2.set_ylim([-0.05,0.05])
ax2.set_ylabel(r'Diff.',size=28)
ax2.set_xlabel(r'Wavelength ($\AA$)',size=28)

plt.savefig('spec.pdf')
#plt.show()