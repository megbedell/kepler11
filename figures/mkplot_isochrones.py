import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb

data = q2.Data("../data/K11_stars.csv","../data/K11_lines.csv")
sun = q2.Star("Sun")
sun.get_data_from(data)
K11 = q2.Star("K11")
K11.get_data_from(data)
K11.teff = 5836.0
K11.err_teff = 7.0
K11.logg = 4.44
K11.err_logg = 0.023
K11.feh = 0.062
K11.err_feh = 0.007

# set up figure:
c1 = 'black'
c2 = '#003399' # blue
c3 = '#CC0033' # red
fig = plt.figure()
ax = plt.gca()
ax.set_xlim([5500,6000])
ax.set_ylim([4.40,4.50])
ax.set_ylabel(r'$\log g$',size=30)
ax.set_xlabel(r'$T_{\rm eff}$',size=30)
#ax.tick_params(which='major', length=10, width=1.5, labelsize=20)
#ax.tick_params(which='minor', length=6, width=1.5, labelsize=20)

# plot isochrones:
feh_offset = 0.0
#ips = q2.isopars.get_isochrone_points(K11, feh_offset, key_parameter_known='logg')
ages = range(1, 7)
text_x = iter(np.arange(5840,5490,-60))
text_y = iter(np.arange(4.485,4.455,-0.0055))
for age in ages:
    iso = q2.isopars.get_isochrone(age, round(K11.feh+feh_offset,2), 'yy02.sql3')
    plt.plot(10.0**iso['logt'], iso['logg'], color=c1, alpha=0.8)
    plt.text(next(text_x), next(text_y), str(age)+' Gyr', size=22, color=c1, alpha=0.8)

# plot data:
plt.errorbar(K11.teff, K11.logg, xerr=K11.err_teff, yerr=K11.err_logg, color=c3, ecolor=c3, fmt='o', markersize=12)
#plt.scatter(sun.teff, sun.logg, color=c2)
plt.savefig('isochrones.pdf')
    
    