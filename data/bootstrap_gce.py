import q2
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pdb

def linear(x, m, b):
     model = m*x + b
     return model

if __name__ == "__main__":
    data = q2.Data("K11_solution.csv","K11_lines.csv")
    sun = q2.Star("Sun")
    sun.get_data_from(data)
    K11 = q2.Star("K11")
    K11.get_data_from(data)
    K11.get_model_atmosphere('odfnew')
    sun.get_model_atmosphere('odfnew')
    
    sp = q2.specpars.SolvePars()
    sp.step_teff = 4
    sp.step_logg = 0.04
    sp.step_vt = 0.04
    sp.niter = 100
    sp.grid = 'odfnew'
    sp.errors = True
    q2.specpars.solve_one(K11, sp, Ref=sun)

    q2.abundances.one(K11, Ref=sun, silent=True, errors=False)
    sp = q2.isopars.SolvePars(key_parameter_known='logg', db='yy02.sql3', feh_offset = 0)
    pp = q2.isopars.PlotPars(make_figures=False)
    q2.isopars.solve_one(K11, sp, PlotPars=pp)
    
    
    # get the UNCORRECTED abundances:
    species_codes = sorted(set(K11.linelist['species']))
    species_ids = q2.abundances.getsp_ids(species_codes)
    Tc = []
    abund = []
    err = []
    for species_id in species_ids:
        Tc = np.append(Tc, q2.abundances.gettc(species_id))
        species_difab = np.array(getattr(K11, species_id)['difab'])
        species_difab = species_difab[species_difab != np.array(None)]  #remove Nones from where ref star was unavailable
        abund = np.append(abund, np.mean(species_difab))
        err = np.append(err, np.std(species_difab) / np.sqrt(len(species_difab) - 1) )
    # TO DO: include stellar parameter uncertainty in errors?
    for t in set(Tc):
        ind = np.where(Tc == t)[0]
        if len(ind) > 1:
            (abund[ind[0]], err[ind[0]]) = np.average(abund[ind], weights=err[ind], returned=True)
            abund = np.delete(abund, ind[1])  #assumes there's only one other state
            err = np.delete(err, ind[1])  #assumes there's only one other state
            Tc = np.delete(Tc, ind[1])  #assumes there's only one other state
            species_ids = np.delete(species_ids, ind[1])

    # remove K:
    abund = np.delete(abund, 7)
    err = np.delete(err, 7)
    Tc = np.delete(Tc, 7)
    species_ids = np.delete(species_ids, 7)
    
    # bootstrap GCE:
    n_iter = 1000
    ages = np.random.choice(K11.pdf_age['x'], size=n_iter, p=K11.pdf_age['y']/np.sum(K11.pdf_age['y']))
    slope = np.zeros(n_iter, dtype='f16')
    intercept = np.zeros(n_iter, dtype='f16')
    # fetch some GCE quantities
    bs = np.asarray([ q2.gce.getb_linear(species_id)[0] for species_id in species_ids ])
    b_errs = np.asarray([ q2.gce.getb_linear(species_id)[1] for species_id in species_ids ])
    Ref_age = 4.6 #Sun's age in Gyr
    for i in range(n_iter):
        b_i = np.zeros_like(bs)
        for j in range(len(b_i)):
            if b_errs[j] == 0:
                continue
            else:
                b_i[j] = np.random.normal(bs[j], b_errs[j])
        abund_gce = abund - (ages[i] - Ref_age)*b_i
        popt, pcov = curve_fit(linear, Tc, abund_gce, sigma=err, absolute_sigma=True)
        slope[i] = popt[0]
        intercept[i] = popt[1]
        if (i % 100 == 0):
            print "{0}/{1} iterations done".format(i, n_iter)
            print "latest fit: m = {0}, b = {1}".format(slope[i], intercept[i])
    print "GCE-corrected Tc slope = {0} +/- {1}".format(np.median(slope),np.std(slope))
    plt.clf()
    plt.scatter(ages, slope*1.0e5)
    plt.ylim([-6,3])
    plt.xlabel('Age (Gyr)')
    plt.ylabel(r'$T_c$ slope after GCE correction x 10^5')
    plt.savefig('age-Tc.png')
    plt.clf()
    n, bins, patches = plt.hist(ages, 50, normed=1)
    plt.plot(K11.pdf_age['x'], K11.pdf_age['y'])
    plt.xlabel('Age (Gyr)')
    plt.savefig('ages.png')
    plt.clf()
    n, bins, patches = plt.hist(slope*1.0e5, 10, normed=1)
    plt.xlabel(r'$T_c$ slope after GCE correction x 10^5')
    plt.savefig('tc.png')
    plt.clf()
    np.savetxt('Tc_fits.txt',np.transpose([ages, slope, intercept]))
        
