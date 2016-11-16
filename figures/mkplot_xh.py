import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
from numpy import genfromtxt

def xh(el_name, data):
    xh = a[el_name+'_1'][1:]
    return xh
    
def xh_err(el_name, data):
    err = a['e_'+el_name+'_1'][1:]
    return err

if __name__ == "__main__":

    a = genfromtxt('../data/K11_abundances_err.csv', delimiter=',', dtype=None, names=True)

    a['e_KI_1'] *= np.sqrt(2)  # adjust for the duplicate potassium lines  

    #(a['C_1'], a['e_C_1']) = np.average([xh('CI',a), xh('CH',a)], weights=[xh_err('CI',a), xh_err('CH',a)], returned=True, axis=0)
    #a_elements = np.asarray(['CI','OI','NaI','MgI','AlI','SiI','SI','KI','CaI','ScI','ScII','TiI','TiII','VI','CrI','CrII','MnI','FeI','FeII','CoI','NiI','CuI','ZnI','YII','CH'])
    elements = np.asarray(['C','O','Na','Mg','Al','Si','S','K','Ca','Sc','Ti','V','Cr','Mn','Co','Ni','Cu','Zn','Y'])
    fehs, feh_errs = np.average([xh('FeI',a), xh('FeII',a)], weights=[xh_err('FeI',a), xh_err('FeII',a)], returned=True, axis=0)
    xhs = np.zeros((len(elements),len(fehs)))
    xh_errs = np.zeros_like(xhs)
    for i,element in enumerate(elements):
        #print element
        if element == 'C':
            xhs[i,:], xh_errs[i,:] = np.average([xh('CI',a), xh('CH',a)], weights=[xh_err('CI',a), xh_err('CH',a)], returned=True, axis=0)
        elif element == 'Sc':
            xhs[i,:], xh_errs[i,:] = np.average([xh('ScI',a), xh('ScII',a)], weights=[xh_err('ScI',a), xh_err('ScII',a)], returned=True, axis=0)
        elif element == 'Ti':
            xhs[i,:], xh_errs[i,:] = np.average([xh('TiI',a), xh('TiII',a)], weights=[xh_err('TiI',a), xh_err('TiII',a)], returned=True, axis=0)
        elif element == 'Cr':
            xhs[i,:], xh_errs[i,:] = np.average([xh('CrI',a), xh('CrII',a)], weights=[xh_err('CrI',a), xh_err('CrII',a)], returned=True, axis=0)
        elif element == 'Y':
            xhs[i,:] = xh(element+'II', a)    
            xh_errs[i,:] = xh_err(element+'II', a)
        else:
            xhs[i,:] = xh(element+'I', a)    
            xh_errs[i,:] = xh_err(element+'I', a)    
            

    fig = plt.figure()
    c1 = '#003399'  #blue
    c2 = '#CC0033'  #red
    matplotlib.rc('xtick', labelsize=12) 
    matplotlib.rc('ytick', labelsize=12)
    thin = [1,2,3,5,6,7,8]  # index of thin-disk stars in xhs etc
    thick = [4,9]  # index of thick-disk stars in xhs etc
    loc = 0
    for i in range(19):
        ax = fig.add_subplot(5, 4, loc + 1)
        xfe_errs = np.sqrt(feh_errs**2 + xh_errs[i,:]**2)
        ax.errorbar(fehs[thin],xhs[i,thin]-fehs[thin],xerr=feh_errs[thin],yerr=xfe_errs[thin],color=c2,mec=c2,fmt='o',markersize=6, label='thin disk')
        ax.errorbar(fehs[thick],xhs[i,thick]-fehs[thick],xerr=feh_errs[thick],yerr=xfe_errs[thick],color=c1,mec=c1,fmt='o',markersize=6, label='thick disk')
        ax.errorbar(fehs[0],xhs[i,0]-fehs[0],xerr=feh_errs[0],yerr=xfe_errs[0],color='black',mec='black',markeredgewidth=1.,fmt='*',markersize=12, label='Kepler-11')
        ax.set_ylabel('[{0}/Fe]'.format(elements[i]),fontsize=20)
        majorLocator = MultipleLocator(0.2)
        minorLocator = MultipleLocator(0.05)
        if elements[i] in ['C', 'O']:
            majorLocator = MultipleLocator(0.4)
            minorLocator = MultipleLocator(0.1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        
        majorLocator = MultipleLocator(0.4)
        minorLocator = MultipleLocator(0.1)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.tick_params('both', length=5, width=2, which='major', pad=5)
        ax.tick_params('both', length=3, width=1, which='minor')
        
        if i < 15:
            ax.set_xticklabels( () )
        else:
            ax.set_xlabel('[Fe/H]',fontsize=20)
        if i == 2:
            loc += 1  # skip this subplot
        loc += 1
        
    fig.subplots_adjust(left=0.2, wspace=0.7, hspace=0.1)
    #fig.tight_layout()
    #fig.text(0.55, 0.04, '[Fe/H]', ha='center', va='center', fontsize=28)
    Artist1 = plt.Line2D((0,1),(0,0), color=c1, mec=c1, marker='o',markersize=10, linestyle='')
    Artist2 = plt.Line2D((0,1),(0,0), color=c2, mec=c2, marker='o',markersize=10, linestyle='')
    Artist3 = plt.Line2D((0,1),(0,0), color='black', mec='black',markeredgewidth=1., marker='*',markersize=14, linestyle='')
    fig.legend([Artist1,Artist2,Artist3],['thick disk','thin disk','Kepler-11'],loc='upper right',borderpad=0.2,borderaxespad=1.0, fontsize=22)
    plt.savefig('xh.pdf')

'''
print r"Element & $N_lines$ & {s[1]} & {s[2]} & {s[3]} & {s[4]} & {s[5]} & {s[6]} & {s[7]} & {s[8]} & {s[9]} & {s[10]} \\".format(s=a['id'])
print r"\hline"
for el in elements:
     print r"{el} & {n} & {s[1]:.3f} $\pm$ {e[1]:.3f} & {s[2]:.3f} $\pm$ {e[2]:.3f} & {s[3]:.3f} $\pm$ {e[3]:.3f} & {s[4]:.3f} \
     $\pm$ {e[4]:.3f} & {s[5]:.3f} $\pm$ {e[5]:.3f} & {s[6]:.3f} $\pm$ {e[6]:.3f} & {s[7]:.3f} $\pm$ {e[7]:.3f} & {s[8]:.3f} \
     $\pm$ {e[8]:.3f} & {s[9]:.3f} $\pm$ {e[9]:.3f} & {s[10]:.3f} $\pm$ {e[10]:.3f} \\".format( \
     el=el, n=a['n_'+el][0], s=a[el+'_1'], e=a['e_'+el+'_1'])
'''

