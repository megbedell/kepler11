import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20

data = np.recfromcsv('../data/exoplanets_org_mr.csv', filling_values=np.nan) 
err_scale = np.sqrt((data.massupper + data.masslower)**2/data.mass**2 + (data.rupper + data.rlower)**2/data.r**2) 

err_range = np.nanmax(err_scale) - np.nanmin(err_scale)
alphas = np.exp(-err_scale)
colors = np.asarray([(0,0,0, alpha) for alpha in alphas])

plt.scatter(data.mass*317.83, data.r*11.209, c=colors, edgecolors=colors)
for ind in np.where(np.isfinite(err_scale))[0]:
    xerr = np.max([data.masslower[ind] * 317.83,data.massupper[ind] * 317.8])
    yerr = np.max([data.rlower[ind] * 11.209, data.rupper[ind] * 11.209])
    plt.errorbar(data.mass[ind]*317.83,data.r[ind]*11.209, xerr=xerr, yerr=yerr, lw=2, color=colors[ind], capsize=5, capthick=2, fmt='o')
    
m_all = np.arange(1000)/10.0  # in grams
r_p05 = (m_all * 5.972e27/(4.0/3.0*np.pi*0.5))**(1.0/3.0) / 6.371e8 # radius curve for rho = 0.5 g/cm3
r_p1 = (m_all * 5.972e27/(4.0/3.0*np.pi*1.0))**(1.0/3.0) / 6.371e8 # radius curve for rho = 1.0 g/cm3
r_p2 = (m_all * 5.972e27/(4.0/3.0*np.pi*2.0))**(1.0/3.0) / 6.371e8 # radius curve for rho = 2.0 g/cm3
r_p4 = (m_all * 5.972e27/(4.0/3.0*np.pi*4.0))**(1.0/3.0) / 6.371e8 # radius curve for rho = 4.0 g/cm3
r_p8 = (m_all * 5.972e27/(4.0/3.0*np.pi*8.0))**(1.0/3.0) / 6.371e8 # radius curve for rho = 8.0 g/cm3
r_p16 = (m_all * 5.972e27/(4.0/3.0*np.pi*16.0))**(1.0/3.0) / 6.371e8 # radius curve for rho = 16.0 g/cm3


plt.plot(m_all,r_p05,ls='--',color=(0,0,0,0.3))
plt.plot(m_all,r_p1,ls='--',color=(0,0,0,0.3))
plt.plot(m_all,r_p2,ls='--',color=(0,0,0,0.3))
plt.plot(m_all,r_p4,ls='--',color=(0,0,0,0.3))
plt.plot(m_all,r_p8,ls='--',color=(0,0,0,0.3))
plt.plot(m_all,r_p16,ls='--',color=(0,0,0,0.3))
    
K11_name = np.asarray(['K-11 b','K-11 c','K-11 d','K-11 e','K-11 f'])
K11_mass_ttv = np.asarray([1.9,2.9,7.3,8.0,2.0])
K11_massupper =np.asarray([1.4,2.9,0.8,1.5,0.8])
K11_masslower = np.asarray([1.0,1.6,1.5,2.1,0.9])
K11_r_ttv = np.asarray([1.8,2.87,3.12,4.19,2.49])
K11_rupper = np.asarray([.03,.05,.06,.07,.04])
K11_rlower = np.asarray([.05,.06,.07,.09,.07])
K11_mass_spec = K11_mass_ttv * 1.04/0.961
K11_r_spec = K11_r_ttv * 1.00/1.065

c1 = '#003399'
c2 = '#CC0033'

plt.errorbar(K11_mass_ttv,K11_r_ttv, xerr=[K11_masslower,K11_massupper], yerr=[K11_rlower,K11_rupper], lw=2.5, color=c1, mec=c1, capsize=5, capthick=2.5, fmt='o', markersize=7, label='TTV')
plt.errorbar(K11_mass_spec,K11_r_spec, xerr=[K11_masslower,K11_massupper], yerr=[K11_rlower,K11_rupper], lw=2.5, color=c2, mec=c2, capsize=5, capthick=2.5, fmt='o', markersize=7, label='spectroscopic')


plt.xlim(0.5,20.0) 
plt.ylim(0.0,5.0) 

plt.xlabel(r'Mass ($M_{\oplus}$)', size=24) 
plt.ylabel(r'Radius ($R_{\oplus}$)', size=24)

plt.text(9,5.1,r'0.5 g cm$^{-3}$',color=(0,0,0,0.5),size=20)
plt.text(20.5,4.7,r'1 g cm$^{-3}$',color=(0,0,0,0.5),size=20)
plt.text(20.5,3.8,r'2 g cm$^{-3}$',color=(0,0,0,0.5),size=20)
plt.text(20.5,2.9,r'4 g cm$^{-3}$',color=(0,0,0,0.5),size=20)
plt.text(20.5,2.4,r'8 g cm$^{-3}$',color=(0,0,0,0.5),size=20)
plt.text(16.5,1.2,'16 \n'r'g cm$^{-3}$',color=(0,0,0,0.5),size=20,ha='center')


plt.xscale('log')
ax = plt.gca()
majorFormatter = FormatStrFormatter('%d')
ax.xaxis.set_major_formatter(majorFormatter)
 
plt.legend(loc='upper left',prop={'size':24})

plt.savefig('K11_massradius.pdf')
