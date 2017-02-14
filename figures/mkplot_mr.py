import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from bokeh.plotting import *
from bokeh.models import HoverTool, PanTool, BoxZoomTool, ResetTool
from scipy.io.idl import readsav



#data = np.recfromcsv('exoplanets_org_mr.csv', filling_values=np.nan)

import xml.etree.ElementTree as ET, urllib, gzip, io
url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urllib.urlopen(url).read())))

# Output mass and radius of all planets 
mass = []
radius = []
name = []
mass_errm = []
mass_errp = []
radius_errm = []
radius_errp = []
for planet in oec.findall(".//planet[mass]"):
    try: # check that all error bars exist
        mem = planet.find("mass").attrib['errorminus']
        mep = planet.find("mass").attrib['errorplus']
        rem = planet.find("radius").attrib['errorminus']
        rep = planet.find("radius").attrib['errorplus']
    except: # if not, skip this planet
        continue
    if (planet.findtext("name") == 'Kepler-11 c'):
        # skip this, mass errors are bad & will re-plot anyway
        continue
    # if yes, save its relevant stats
    mass =  np.append(mass, float(planet.findtext("mass")))
    radius = np.append(radius, float(planet.findtext("radius")))
    name = np.append(name, planet.findtext("name"))
    mass_errm = np.append(mass_errm, float(mem))
    mass_errp = np.append(mass_errp, float(mep))
    radius_errm = np.append(radius_errm, float(rem))
    radius_errp = np.append(radius_errp, float(rep))

'''    
save_cat =  np.transpose(np.asarray([name, mass, mass_errm, mass_errp, radius, \
            radius_errm, radius_errp]))
np.savetxt('mr_from_oec.csv', save_cat, \
            delimiter=',', fmt='%s', header='name, mass, mass_errm, mass_errp, radius, radius_errm, radius_errp')
'''      

models = readsav('structure_models.dat')

err_scale = np.sqrt((mass_errp + mass_errm)**2/mass**2 + (radius_errp + radius_errm)**2/radius**2) 

#err_range = np.nanmax(err_scale) - np.nanmin(err_scale)
alphas = np.exp(-err_scale*1.5)
colors = np.asarray([(0,0,0, alpha) for alpha in alphas])


plt.scatter(mass*317.83, radius*11.209, c=colors, edgecolors=colors)
for ind in np.where(np.isfinite(err_scale))[0]:
    xerr = np.max([mass_errm[ind] * 317.83,mass_errp[ind] * 317.8])
    yerr = np.max([radius_errm[ind] * 11.209, radius_errp[ind] * 11.209])
    plt.errorbar(mass[ind]*317.83,radius[ind]*11.209, xerr=xerr, yerr=yerr, lw=2, color=colors[ind], capsize=5, capthick=2, fmt='o')


plt.plot(models.hydrogen_mass, models.hydrogen_radius,ls='--', color='#7A68A6') # purple
plt.text(4,5.1,r'H/He',color='#7A68A6',size=22)
plt.plot(models.water_mass, models.water_radius,ls='--', color='#348ABD') # blue
plt.text(20.5,3.0,r'water world',color='#348ABD',size=22)
plt.plot(models.silicate_mass, models.silicate_radius,ls='--', color='#188487') # turquoise
plt.text(20.5,2.3,r'Earth-like',color='#188487',size=22)
plt.plot(models.iron_mass, models.iron_radius,ls='--', color='#467821') # green
plt.text(20.5,1.6,r'Mercury-like',color='#467821',size=22)


K11_name = np.asarray(['K-11 b','K-11 c','K-11 d','K-11 e','K-11 f'])
#K11_mass_ttv = np.asarray([1.9,2.9,7.3,8.0,2.0])
#K11_massupper =np.asarray([1.4,2.9,0.8,1.5,0.8])
#K11_masslower = np.asarray([1.0,1.6,1.5,2.1,0.9])
#K11_r_ttv = np.asarray([1.8,2.87,3.12,4.19,2.49])
#K11_rupper = np.asarray([.03,.05,.06,.07,.04])
#K11_rlower = np.asarray([.05,.06,.07,.09,.07])
#K11_mass_spec = K11_mass_ttv * 1.04/0.961
#K11_r_spec = K11_r_ttv * 1.008/1.065
K11_mass_ttv = np.asarray([2.78,5.0,8.13,9.48,2.53])
K11_massupper =np.asarray([0.64,1.3,0.67,0.86,0.49])
K11_masslower = np.asarray([0.66,1.35,0.66,0.88,0.45])
K11_r_ttv = np.asarray([1.83,2.89,3.21,4.26,2.54]) * 1.008
K11_rupper = np.asarray([.07,.12,.12,.16,.10])
K11_rlower = np.asarray([.04,.04,.04,.07,.04])
K11_mass_spec = K11_mass_ttv * 1.04/0.961
K11_r_spec = K11_r_ttv * 1.008/1.065
#K11_mass_spec = np.asarray([2.83,5.05,7.52,8.37,1.59])
#K11_r_spec = np.asarray([1.71,2.694,3.00,3.947,2.356])

c1 = '#003399'
c2 = '#CC0033'

plt.errorbar(K11_mass_ttv,K11_r_ttv, xerr=[K11_masslower,K11_massupper], yerr=[K11_rlower,K11_rupper], lw=2.5, color=c1, capsize=5, capthick=2.5, fmt='o', markersize=7, label='TTV')
plt.errorbar(K11_mass_spec,K11_r_spec, xerr=[K11_masslower,K11_massupper], yerr=[K11_rlower,K11_rupper], lw=2.5, color=c2, capsize=5, capthick=2.5, fmt='o', markersize=7, label='spectroscopic')


plt.xlabel(r'Mass ($M_{\oplus}$)', size=24) 
plt.ylabel(r'Radius ($R_{\oplus}$)', size=24)

plt.xscale('log')
plt.xlim(0.5,20.0) 
plt.ylim(0.0,5.0) 
#x0, x1, y0, y1 = plt.axis()
#plot_margin = 0.4
#plt.axis((x0 - plot_margin, x1 + plot_margin,y0,y1))
ax = plt.gca()
majorFormatter = FormatStrFormatter('%d')
ax.xaxis.set_major_formatter(majorFormatter)

plt.legend(numpoints=1,loc='upper left',prop={'size':24})

plt.savefig('K11_massradius.pdf',bbox_inches='tight')
