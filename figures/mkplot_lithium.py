import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb

spec_file = '../data/K11_2d.dat'
synth_file = '../data/K11.li'
sun_file = '../data/Sun_li.dat'
sun_synth_file = '../data/solar.li'
spec_wave, spec = np.loadtxt(spec_file, unpack=True)
sun_wave, sun_spec = np.loadtxt(sun_file, unpack=True)
synth_wave, synth = np.loadtxt(synth_file, unpack=True, skiprows=2)
synth_wave, sun_synth = np.loadtxt(sun_synth_file, unpack=True, skiprows=2)
wave_factor = 0.00 # wavelength shift = second-to-last number on second-to-last line of synth pars
c_factor = 1.0205 # continuum factor = last number on second-to-last line of synth pars

c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
fig = plt.figure()
ax = plt.gca()
ax.set_xlim([6706.8,6708.4])
ax.set_ylim([0.94,1.02])
ax.ticklabel_format(useOffset=False)
ax.tick_params(which='major', length=10, width=1.5, labelsize=22)
ax.tick_params(which='minor', length=6, width=1.5)
majorLocator = MultipleLocator(0.05)
minorLocator = MultipleLocator(0.01)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(1.0)
minorLocator = MultipleLocator(0.1)
ax.xaxis.set_minor_locator(minorLocator)
ax.xaxis.set_major_locator(majorLocator)
plt.scatter(spec_wave+wave_factor, spec*c_factor, color=c1, label=r'Kepler-11')
plt.plot(synth_wave, sun_synth, color=c2, lw=2, label=r'A(Li)$_{\odot}$ = +1.03')
plt.plot(synth_wave, synth, color=c3, lw=2, label=r'A(Li) = +1.28')
ax.set_ylabel(r'Normalized Flux',size=28)
ax.set_xlabel(r'Wavelength ($\AA$)',size=28)
ax.legend(loc='lower right',prop={'size':24})

plt.text(6707.75,1.0,'Li I',ha='center',size=28)
plt.plot(np.zeros(5)+6707.75, np.linspace(0.994,0.997,num=5), color=c1, lw=3)


plt.savefig('lithium.pdf')