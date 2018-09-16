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

Nsim = 5000
N_field_stars = 40

def ReadPFData():
    Glats  = []
    Glongs = []
    Ps = []
    fop = open('./PF/table2.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        G_lon = float(sl[2])
        G_lat = float(sl[3])
        P = float(sl[4])
        Perr = float(sl[5])
        #if P/Perr > 3:
        Ps.append(P)
        Glongs.append(G_lon)
        Glats.append(G_lat)
    fop.close()
    print 'Using ', len(Ps), ' PF stars' #, ' P/sP > 3 stars'
    return Glongs,Glats,Ps

def log_norm(x, a, x0, sigma):
    return a/(x)*np.exp(-(np.log(x)-x0)**2/(2*sigma**2))

def norm(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def findSigmaPDist(Glongs,Glats,Ps):
    """
    Drops a circle with area equal to the average UFO 95% ellipse 0.0456 deg^2 on the PF region
    then calculates sigmaP of stars within the circle if there are more than 5 stars inside
    Finally it fits the distribution of sigmaPs and returns sigma of this distribution
    """
    G_Lat_Low = 24.6
    G_Lat_Upp = 28.0
    G_Lon_Low = 122.0
    G_Lon_Upp = 126.0
    Radius = math.sqrt(0.0456/math.pi)
    
    std_Ps = []  # list of standard deviations of polarization within the circle at different positions of the PF
    i = 0
    while i < 3000:
        G_Lat = random.uniform(G_Lat_Low, G_Lat_Upp)
        G_Lon = random.uniform(G_Lon_Low, G_Lon_Upp)
        curPs = []
        for m in range(len(Glongs)):
            if math.sqrt((Glongs[m]-G_Lon)**2 + (Glats[m]-G_Lat)**2) < Radius: # radius is small so we don't care about spherical geometry
                curPs.append(Ps[m])
        # now make sure that we have enough stars within the circle to find sigma
        if len(curPs) > 4:
            i = i + 1
            std_Ps.append(np.std(curPs))
    
    '''
    Fits the distribution of polarization degree observed in Polaris Flare
    '''
    weights = np.ones_like(std_Ps)/len(std_Ps)
    n, bins, patches   = plt.hist(std_Ps,bins=33,normed=True,histtype='bar',color=['red'],label='mean '+str(np.mean(std_Ps)),alpha=0.6)

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
    plt.savefig('PF_sigma_pol_dist.png')
    plt.close()
    plt.cla()
    plt.clf()

    return popt[2]

def GenISP(Plevels, sigma):
    """
    Generates N_field_stars for every level of average polarization and finds the mean and std for them
    This function almost makes no sence, just introduces more noise (but we are repeating what has been done before).
    """
    fields_params = []
    base_mean = np.mean(Ps)
    for p in Plevels:
        gen_pols = stats.norm.rvs(loc=p, scale=sigma, size=N_field_stars)
        positive_gen_pols = []
        for gp in gen_pols:
            if gp < 0:
                continue
            else:
                positive_gen_pols.append(gp)
        fields_params.append([np.mean(positive_gen_pols),np.std(positive_gen_pols)])
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

def Plot(Plevels,detect_fract_5s,detect_fract_3s):
    fig, ax = plt.subplots()
    fig.set_size_inches(4.5,3)
    font = {'size'   : 24, 'family' : 'sans'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=False)
    
    plt.plot(Plevels,detect_fract_5s,marker='o',color="#5d76cb",linestyle = 'None',label='$5\sigma$')
    plt.plot(Plevels,detect_fract_3s,marker='o',color="#cd5b45",linestyle = 'None',label='$3\sigma$')
    
    ax.set_xlabel('Field stars <P> (%)')
    ax.set_ylabel('Detectability fraction')
    ax.set_xticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlim(0,8)
    plt.legend(prop={'size': 14})
    plt.savefig('Detect_Av_Pol.eps',bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    Glongs,Glats,Ps = ReadPFData()
    sigma = findSigmaPDist(Glongs,Glats,Ps)
    Plevels = np.arange(0.1,8,0.2)
    fields_params = GenISP(Plevels, sigma)
    Pblaz = ReadBlazarCat()
    Pfakeblaz = GetRandBlazars(len(Plevels)*Nsim,Pblaz)
    
    detect_fract_5s = np.zeros(len(Plevels))
    detect_fract_3s = np.zeros(len(Plevels))
    for i in range(1,Nsim+1):
        for k in range(len(Plevels)):
            mean, sigma = fields_params[k]
            if 100 * Pfakeblaz[i*k] > (mean + 5 * sigma):
                detect_fract_5s[k] += 1
            if 100 * Pfakeblaz[i*k] > (mean + 3 * sigma):
                detect_fract_3s[k] += 1
    Plot(Plevels,map(lambda x: x/float(Nsim),detect_fract_5s),map(lambda x: x/float(Nsim),detect_fract_3s))
        




