import q2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.gridspec as gridspec
import numpy as np
import pdb

def linear(p, fjac=None, x=None, y=None, err=None):
     m,b = (p[0],p[1])
     model = m*x + b
     status = 0
     return([status, (y-model)/err])

def linbreak(p, fjac=None, x=None, y=None, err=None):
     m,offs,bk=(p[0],p[1],p[2])
     model = np.piecewise(x, [x < bk, x >= bk], [offs, lambda x: (offs + m*(x - bk))])
     status = 0
     return([status, (y-model)/err])

if __name__ == "__main__":
    
    Tc = np.array([   40.    ,   180.    ,   958.    ,  1336.    ,  1653.    , \
        1310.    ,   664.    ,  1006.    ,  1517.    ,  1659.    , \
        1582.    ,  1429.    ,  1296.    ,  1158.    ,  1334.    , \
        1352.    ,  1353.    ,  1037.    ,   726.    ,  1659.0001])
    abund = np.array([ 0.02574641,  0.05433333,  0.04682   ,  0.03035   ,  0.0618    , \
        0.04720214,  0.034655  ,  0.0674    ,  0.06088909,  0.09034457, \
        0.06480416,  0.08120111,  0.04741255,  0.0596175 ,  0.05746099, \
        0.06522167,  0.0639525 ,  0.078205  ,  0.02857   ,  0.065288  ])
    err = np.array([ 0.04336853,  0.00856997,  0.00930914,  0.02530275,  0.0078    , \
            0.00501722,  0.0229398 ,  0.05      ,  0.00523373,  0.02910326, \
            0.01357965,  0.00599431,  0.01702205,  0.00846271,  0.00835767, \
            0.00362001,  0.00471549,  0.0071002 ,  0.02918747,  0.01093004])
    abund_gce = np.array([-0.01242189,  0.02701096, -0.01038549,  0.00763424,  0.02946838, \
        0.03026747,  0.01826281,  0.0674    ,  0.04935534,  0.05814537, \
        0.05325056,  0.06992715,  0.04025283,  0.0392135 ,  0.05746099, \
        0.03102024,  0.03283167,  0.03001643, -0.00644105,  0.065288  ])
    err_gce = np.array([ 0.04690626,  0.01496574,  0.04262738,  0.0280323 ,  0.01866978, \
        0.01164087,  0.02515044,  0.05      ,  0.01007805,  0.03571811, \
        0.0146558 ,  0.01118995,  0.01755926,  0.01834229,  0.00835767, \
        0.02410788,  0.02344593,  0.02736322,  0.0356984 ,  0.01093004])
    
    Tc = np.delete(Tc, 7)
    abund = np.delete(abund, 7)
    err = np.delete(err, 7)
    abund_gce = np.delete(abund_gce, 7)
    err_gce = np.delete(err_gce, 7)
    
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
    ax.set_ylim([-0.07,0.12])
    majorLocator = MultipleLocator(500)
    minorLocator = MultipleLocator(100)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    ax.set_xlim([0,1800])
    plt.errorbar(Tc,abund,yerr=err,color=c1,mec=c1,fmt='o',markersize=10)
    plt.errorbar(Tc,abund_gce,yerr=err_gce,color=c2,mec=c2,fmt='o',markersize=10)
    x_all = np.arange(1800)*1.0
    slope,offs,bk= (9.48695845e-06,   4.97951530e-02,   0.00000000e+00)
    plt.plot(x_all, np.piecewise(x_all,[x_all<bk,x_all>=bk],[offs,lambda x_all: (offs + slope*(x_all - bk))]), color=c1, linewidth=2)
    slope,offs,bk=(4.29898056e-05,   1.95984843e-02,   7.26000003e+02)
    plt.plot(x_all, np.piecewise(x_all,[x_all<bk,x_all>=bk],[offs,lambda x_all: (offs + slope*(x_all - bk))]), color=c2, linewidth=2)
    plt.xlabel(r'$T_c$ (K)', fontsize=26)
    plt.ylabel('[X/H]', fontsize=26)
    
    fig.savefig('K11_Tc_break.pdf')
    