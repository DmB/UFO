#!/usr/bin/python

"""
This simulation repeats the algorithm from Kokolakis thesis.
It estimates dependence of Percentage of registered blazars vs average ISP in the field.
"""

from scipy import stats
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math, random
import distr

N_field_stars = 40

def ReadPFData():
    Ps = []
    fop = open('./PF/table2.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        P = float(sl[4])
        Perr = float(sl[5])
        if P/Perr > 0.1:
            Ps.append(P)
    fop.close()
    print 'Using ', len(Ps), ' P/sP > 0.1 stars'
    return Ps

def log_norm(x, a, x0, sigma):
    return a/(x)*np.exp(-(np.log(x)-x0)**2/(2*sigma**2))

def norm(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def FitDist(Ps):
    '''
    Fits the distribution of polarization degree observed in Polaris Flare
    '''
    weights = np.ones_like(Ps)/len(Ps)
    n, bins, patches   = plt.hist(Ps,bins=33,normed=True,histtype='bar',color=['red'],label='mean '+str(np.mean(Ps)),alpha=0.6)

    x = [0.5*(a+b) for a,b in zip(bins,bins[1:])]
    popt, pcov = curve_fit(norm, x, n, p0 = [1, 1, 0.24], maxfev=1000)
    plt.plot(x,norm(x,*popt))


    # plot dist with three times larger mean testing purposes
    gendata = stats.norm.rvs(loc=popt[1]*3, scale=popt[2], size=200)
    n, bins, patches = plt.hist(gendata,bins=33,normed=True,histtype='bar',color=['green'],label='mean '+str(np.mean(gendata)),alpha=0.6)
    
    # plot dist with 5 times larger mean for testing purposes
    gendata2 = stats.norm.rvs(loc=popt[1]*5, scale=popt[2], size=200)
    n, bins, patches = plt.hist(gendata2,bins=33,normed=True,histtype='bar',color=['yellow'],label='mean '+str(np.mean(gendata2)),alpha=0.6)

    plt.grid()
    plt.legend()
    plt.xlabel('P [%]')
    plt.ylabel('Fraction')
    plt.savefig('PF_pol_dist.png')
    plt.close()
    plt.cla()
    plt.clf()
    return popt

def GenISP(Plevels, mu, sigma):
    """
    Generates N_field_stars for every level of average polarization and finds the mean and std for them
    This function almost makes no sence, just introduces more noise (but we are repeating what has been done before).
    """
    fields_params = []
    base_mean = np.mean(Ps)
    for p in Plevels:
        gen_pols = stats.norm.rvs(loc=p, scale=sigma, size=N_field_stars)
        fields_params.append([np.mean(gen_pols),np.std(gen_pols)])
    return fields_params

def ReadBlazarCat():
    """
    Reads catalog from Angelakis et al. MNRAS, 463, 3365 (2016)
    and return list of tuples ( average_polar, modulation_index )
    """
    Pblaz = []
    fop = open('RoboPol_mean-polarization.txt')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        if sl[5] != 'NA':
            p0 = float(sl[5])
            m  = float(sl[9])
            Pblaz.append((p0,m))
    fop.close()
    return Pblaz

def GetRandBlazars(N,Pblaz):
    '''
    Looks into catalog by Angelakis et al. MNRAS, 463, 3365 (2016) 
    And generates a list of N polarizations according to it

    N:     is the number of generated polarizations needed
    Pblaz: is a list of tuples ( average_polar, modulation_index ) from Angelakis et al.
    '''
    Ps = []
    for i in range(N):
        p0,m0 = random.choice(Pblaz)
        alpha = ((1.0-p0)/(p0*m0**2) - 1) * p0
        beta  = ((1.0-p0)/(p0*m0**2) - 1) * (1 - p0)
        P = np.random.beta(alpha, beta, size=1)[0]
        Ps.append(P)
    return Ps

def Plot(Plevels,detect_fract):
    fig, ax = plt.subplots()
    fig.set_size_inches(5,3.333)
    font = {'size'   : 22, 'family' : 'sans'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=False)
    ax.set_xlabel('<P> (%)')
    ax.set_ylabel('Detectability fraction')
    plt.plot(Plevels,detect_fract,marker='o',linestyle = 'None')
    plt.savefig('Detect_Av_Pol.eps',bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    Ps = ReadPFData()
    fact, mu, sigma = FitDist(Ps)
    Plevels = np.arange(0.1,8,0.2)
    fields_params = GenISP(Plevels, mu, sigma)
    Pblaz = ReadBlazarCat()
    Ps = GetRandBlazars(len(Plevels)*500,Pblaz)
    
    detect_fract = np.zeros(len(Plevels))
    for i in range(1,501):
        for k in range(len(Plevels)):
            mean, sigma = fields_params[k]
            if 100 * Ps[i*k] > mean + 5 * sigma:
                detect_fract[k] += 1
    Plot(Plevels,map(lambda x: x/500.,detect_fract))
        




