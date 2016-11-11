import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb

data = q2.Data("../data/K11_solution.csv","../data/K11_lines.csv")
sun = q2.Star("Sun")
sun.get_data_from(data)
K11 = q2.Star("K11")
K11.get_data_from(data)


# set up figure:
c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
fig = plt.figure()
ax = plt.gca()
ax.set_xlim([6000,5500])
ax.set_ylim([4.50,4.30])
ax.set_ylabel(r'$\log g$',size=30)
ax.set_xlabel(r'$T_{\rm eff}$',size=30)
#ax.tick_params(which='major', length=10, width=1.5, labelsize=20)
#ax.tick_params(which='minor', length=6, width=1.5, labelsize=20)
ax.yaxis.set_ticks(np.arange(4.3, 4.5, 0.05))

# plot isochrones:
feh_offset = 0.0
#ips = q2.isopars.get_isochrone_points(K11, feh_offset, key_parameter_known='logg')
ages = range(1, 10)
#text_x = iter(np.linspace(5840,5090,num=10))
#text_y = iter(np.linspace(4.485,4.425,num=10))
for age in ages:
    iso = q2.isopars.get_isochrone(age, round(K11.feh+feh_offset,2), 'yy02.sql3')
    plt.plot(10.0**iso['logt'], iso['logg'], color=c1, alpha=0.8)
    #plt.text(next(text_x), next(text_y), str(age)+' Gyr', size=22, color=c1, alpha=0.8)

plt.text(5970, 4.48, '1 Gyr', size=22, color=c1)
plt.text(5630, 4.47, '5 Gyr', size=22, color=c1)
plt.text(5600, 4.34, '9 Gyr', size=22, color=c1)
# plot data:
plt.errorbar(K11.teff, K11.logg, xerr=K11.err_teff, yerr=K11.err_logg, color=c3, ecolor=c3, fmt='o', markersize=12, label='this work')
plt.errorbar(5663, 4.366, xerr=60, yerr=0.016, color=c2, ecolor=c2, fmt='o', markersize=12, label='L13')
#plt.scatter(sun.teff, sun.logg, color=c2)
plt.legend(loc='upper left')
plt.savefig('isochrones.pdf')
    
    