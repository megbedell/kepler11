import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb
from scipy.optimize import curve_fit

def linear(x, m, b):
     model = m*x + b
     return model

def linbreak(p, fjac=None, x=None, y=None, err=None):
     m,offs,bk=(p[0],p[1],p[2])
     model = np.piecewise(x, [x < bk, x >= bk], [offs, lambda x: (offs + m*(x - bk))])
     status = 0
     return([status, (y-model)/err])

if __name__ == "__main__":
    
    data = q2.Data("../data/solution_4step.csv","../data/K11_lines.csv")
    sun = q2.Star("Sun")
    sun.get_data_from(data)
    K11 = q2.Star("K11")
    K11.get_data_from(data)
    K11.get_model_atmosphere('odfnew')
    sun.get_model_atmosphere('odfnew')

    #sp = q2.specpars.SolvePars()
    #sp.step_teff = 4
    #sp.step_logg = 0.04
    #sp.step_vt = 0.04
    #sp.niter = 100
    #sp.grid = 'odfnew'
    #sp.errors = True
    #q2.specpars.solve_one(K11, sp, Ref=sun)
    q2.abundances.one(K11, Ref=sun, silent=True, errors=False)

    abund_linear, err, Tc = q2.gce.correct(K11, age=2.7, method='linear', silent=False)
    # exclude K:
    abund_linear = np.delete(abund_linear, 7)
    err = np.delete(err, 7)
    Tc = np.delete(Tc, 7)
    

    matplotlib.rcParams['xtick.labelsize'] = 20
    matplotlib.rcParams['ytick.labelsize'] = 20

    fig = plt.figure()
    c1 = '#003399'  #blue
    c2 = '#CC0033'  #red
    ax = plt.gca()
    majorLocator = MultipleLocator(0.05)
    minorLocator = MultipleLocator(0.01)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.set_ylim([-0.01,0.15])
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    ax.set_xlim([0,1800])
    x_all = np.arange(1800)*1.0
    
    ages, slopes, intercepts = np.loadtxt('../data/tcfit_msboot.txt', unpack=True)
    slope_lo = np.percentile(slopes,16)
    slope_med = np.median(slopes)
    slope_hi = np.percentile(slopes,84)
    int_lo = np.percentile(intercepts,16)
    int_med = np.median(intercepts)
    int_hi = np.percentile(intercepts,84)
    
    plt.fill_between(x_all,linear(x_all, slope_hi, int_lo), linear(x_all, slope_lo, int_hi),\
        color=c1,alpha=0.2) # ASSUMES PERFECT ANTI-CORRELATION

    plt.errorbar(Tc,abund_linear,yerr=err,color=c1,mec=c1,fmt='o',markersize=10)
    plt.plot(x_all, linear(x_all, slope_med, int_med), color=c1, linewidth=2)
    
    #popt, pcov = curve_fit(linear, Tc, abund_linear, sigma=err, absolute_sigma=True)
    #plt.plot(x_all, linear(x_all, popt[0], popt[1]), color=c2, linewidth=2)
    plt.xlabel(r'$T_c$ (K)', fontsize=26)
    plt.ylabel('[X/H]', fontsize=26)
    
    fig.savefig('K11_Tc_linear.pdf')
    