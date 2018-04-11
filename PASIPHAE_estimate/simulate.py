#!/usr/bin/python

import random
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import distr

Nsim = 100

def ReadData():
    PDs = []
    fop = open('unidentPASIPHAEextinct.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        PD = float(sl[3])
        PDs.append(PD)
    fop.close()
    return PDs

#Gaussian function
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def log_norm(x, a, x0, sigma):
    return a/(x)*np.exp(-(np.log(x)-x0)**2/(2*sigma**2))

def log_norm_swapped(x, a, x0, sigma):
    try:
        some_object_iterator = iter(x)
        return a/(map(lambda xx: -1*(xx-25.34), x))*np.exp(-(np.log(map(lambda xx: -1*(xx-25.34), x))-x0)**2/(2*sigma**2))
    except:
        return a/(-1*(x-25.34))*np.exp(-(np.log(-1*(x-25.34))-x0)**2/(2*sigma**2))

def log_norm_swapped_def(x):
    a, x0, sigma = [ 0.83856933,  2.02583607,  0.226895  ]
    return a/(-1*(x-25.34))*np.exp(-(np.log(-1*(x-25.34))-x0)**2/(2*sigma**2))


def FitIdentMags():
    Rs = []
    usno = np.genfromtxt('USNO_ident.csv', delimiter=',',skip_header=1)
    for i in range(len(usno)):
        R = usno[i,40]
        if not np.isnan(R):
            Rs.append(-R + 25.34)

    plt.figure()
    weights = np.ones_like(Rs)/float(len(Rs))
    n, bins = np.histogram(Rs, 36, weights=weights, normed=False)
    
    x = [0.5*(y+x) for x,y in zip(bins,bins[1:])]
    mean = 8
    sigma = 2
    popt, pcov = curve_fit(log_norm, x, n, p0 = [1, mean, sigma])
    x = np.arange(0,20,0.08)
    n, bins, patches = plt.hist(map(lambda x: -1*(x-25.34), Rs), 36, weights=weights, normed=False, histtype='bar',color='green', alpha=0.75)
    plt.plot(x, log_norm_swapped(x, *popt))
    plt.xlabel("R [mag]")
    plt.ylabel("N")
    plt.savefig('dist_mags_ident.png',bbox_inches='tight')
    return popt

def plot(Mags):
    plt.figure()
    weights = np.ones_like(Mags)/float(len(Mags))
    n, bins, patches = plt.hist(Mags, 13, weights=weights, normed=False, histtype='bar',color='green', alpha=0.75)
    plt.xlim([-1,27])
    plt.xlabel("R [mag]")
    plt.ylabel("N")
    plt.savefig('dist_mags_generated.png',bbox_inches='tight')

def getRandPol(N):
    
    return

def getIdent(Mags,PDs,Gen_PDs):
    return 

if __name__ == "__main__":
    PDs = ReadData()
    popt = FitIdentMags()
    x = np.arange(9, 18.0, .01)
    p = log_norm_swapped_def(x)
    for i in range(Nsim):
        a = distr.GeneralRandom( x, p, Nrl = 1000)
        Mags = a.random(len(PDs))[0]
        Gen_PDs = getRandPol(len(PDs))
        getIdent(Mags,PDs,Gen_PDs)
    plot(Mags)
    

