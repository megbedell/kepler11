import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb
from scipy.optimize import curve_fit
from numpy import genfromtxt


def linear(x, m, b):
     model = m*x + b
     return model

def linbreak(p, fjac=None, x=None, y=None, err=None):
     m,offs,bk=(p[0],p[1],p[2])
     model = np.piecewise(x, [x < bk, x >= bk], [offs, lambda x: (offs + m*(x - bk))])
     status = 0
     return([status, (y-model)/err])

def xh(el_name, data):
    xh = a[el_name+'_1'][1:]
    return xh
    
def xh_err(el_name, data):
    err = a['e_'+el_name+'_1'][1:]
    return err
    
def onesigma_region(x_all, slopes, intercepts):
    y_lo, y_med, y_hi = np.empty(len(x_all)), np.empty(len(x_all)), np.empty(len(x_all))
    for i,x in enumerate(x_all):
        sample = linear(x, slopes, intercepts)
        y_lo[i] = np.percentile(sample,16)
        y_med[i] = np.percentile(sample,50)
        y_hi[i] = np.percentile(sample,84)
    return y_lo, y_med, y_hi
    
def get_full_errors(Star, Ref):
    q2.abundances.one(Star, Ref=Ref, silent=True, errors=True) 
    
    err = []
    abund = []
    Tc = []
    species_codes = sorted(set(Star.linelist['species']))
    species_ids = q2.abundances.getsp_ids(species_codes)
    for species_id in species_ids:
        # get abundances:
        Tc = np.append(Tc, q2.abundances.gettc(species_id))
        species_difab = np.array(getattr(Star, species_id)['difab'])
        species_difab = species_difab[species_difab != np.array(None)]  #remove Nones from where ref star was unavailable
        abund = np.append(abund, np.mean(species_difab))
        err = np.append(err, getattr(Star, species_id)['err_difab'])
        
    for t in set(Tc):
        # eliminate duplicate measurements of the same element
        ind = np.where(Tc == t)[0]
        if len(ind) == 2:  # won't take care of 3+ states of the same thing
            (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
            abund = np.delete(abund, ind[1])
            err = np.delete(err, ind[1])
            Tc = np.delete(Tc, ind[1])
            species_ids = np.delete(species_ids, ind[1])
            
    return err
    
    

if __name__ == "__main__":
    
    data = q2.Data("../data/K11_solution.csv","../data/K11_lines.csv")
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
    q2.abundances.one(K11, Ref=sun, silent=True, errors=True) 
    
    K11.OI['difab'] = np.array([ 0.053,  0.048,  0.028])  # manually apply NLTE corrections
    
    matplotlib.rcParams['xtick.labelsize'] = 20
    matplotlib.rcParams['ytick.labelsize'] = 20

    fig = plt.figure()
    c1 = '#003399'  #blue
    c2 = '#CC0033'  #red
    ax = fig.add_subplot(1, 2, 1)
    majorLocator = MultipleLocator(0.05)
    minorLocator = MultipleLocator(0.01)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.set_ylim([-0.01,0.15])
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    ax.set_xlim([-100,1900])
    
    # without GCE:
    x_all = np.arange(2000)*1.0 - 100.0
    

    abund_linear, err, Tc = q2.gce.correct(K11, age=4.6, method='linear', silent=True, errors=True)
    #abund_linear = np.delete(abund_linear, 7)
    #err = np.delete(err, 7)
    #Tc = np.delete(Tc, 7)
    #abund_linear[7] -= 0.019  # LTE correction
    err[7] *= np.sqrt(2)
    #ax.errorbar(Tc,abund_linear,yerr=err,color=c1,mec=c1,fmt='o',markersize=10)
    #print "Scatter-based errors pre-GCE:", err
    
    #err_full = get_full_errors(K11, sun)
    print "Full errors pre-GCE:", err
    eb2 = ax.errorbar(Tc,abund_linear,yerr=err,color=c1,mec=c1,fmt='o',markersize=10)
    #eb2[-1][0].set_linestyle('--')
    
    
    popt, pcov = curve_fit(linear, Tc, abund_linear, sigma=err, absolute_sigma=True)
    plt.plot(x_all, linear(x_all, popt[0], popt[1]), color=c1, linewidth=2)
    
    ax.set_xlabel(r'$T_c$ (K)', fontsize=26)
    ax.set_ylabel('[X/H]', fontsize=26)
    
    
    # with GCE:
    abund_agelo, err, Tc = q2.gce.correct(K11, age=2.3, method='linear', silent=True, errors=True)
    abund_agehi, err, Tc = q2.gce.correct(K11, age=4.1, method='linear', silent=True, errors=True)
    abund_linear, err, Tc = q2.gce.correct(K11, age=3.2, method='linear', silent=True, errors=True)
    age_err = np.max([np.abs(abund_linear - abund_agelo), np.abs(abund_linear - abund_agehi)])
    # exclude K:
    #abund_linear = np.delete(abund_linear, 7)
    #err = np.delete(err, 7)
    #Tc = np.delete(Tc, 7)
    #abund_linear[7] -= 0.019  # LTE correction
    err[7] *= np.sqrt(2)
    
    # get error bars with propagated age uncertainties:
    _, _, err2, _, _ = np.loadtxt('../data/tc_medianabund.txt', unpack=True)
    err = np.insert(err2, 7, err[7]) # insert KI
    print "Errors post-GCE without age uncertainty:", err
    err = np.sqrt(err**2 + age_err**2)
    print "Errors post-GCE without age uncertainty:", err
    

    
    
    ax = fig.add_subplot(1, 2, 2)
    majorLocator = MultipleLocator(0.05)
    minorLocator = MultipleLocator(0.01)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.set_ylim([-0.01,0.15])
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    ax.set_xlim([-100,1900])
    ax.set_yticklabels( () )
    
    
    ages, slopes, intercepts = np.loadtxt('../data/tcfit_msboot.txt', unpack=True)
    y_lo, y_med, y_hi = onesigma_region(x_all, slopes, intercepts)    
    plt.fill_between(x_all, y_lo, y_hi, color=c1,alpha=0.2)
    
    '''''   
    rand = np.random.randint(0,high=len(slopes),size=500)
    rand = rand[(slopes[rand] >= slope_lo) & (slopes[rand] <= slope_hi)]
    for i in rand:
        ax.plot(x_all,linear(x_all,slopes[i],intercepts[i]),color=c1,alpha=0.2)
    '''
    

    ax.errorbar(Tc,abund_linear,yerr=err,color=c1,mec=c1,fmt='o',markersize=10)
    ax.plot(x_all, y_med, color=c1, linewidth=2)
    
    #popt, pcov = curve_fit(linear, Tc, abund_linear, sigma=err, absolute_sigma=True)
    #plt.plot(x_all, linear(x_all, popt[0], popt[1]), color=c2, linewidth=2)
    ax.set_xlabel(r'$T_c$ (K)', fontsize=26)
    fig.subplots_adjust(wspace=0.05,bottom=0.5)

    
    fig.savefig('K11_Tc_linear.pdf')
    plt.clf()
    