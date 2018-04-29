#!/usr/bin/python

import pyfits
import healpy
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
import pickle

"""
Simulation of possibility do descriminate blazar from stars
variant 3

1. We found P(E(B-V)) dependence from 
...
3. We draw a random blazar from the flux limited unbiased sample of RoboPol
4. using appendix from Blinov et al., MNRAS, 457(2):2252-2262, 2016 we transform intrinsic p0 and m0 to random observed pr
5. if the blazar polarization is larger than the field polarization we consider blazar to be detectable
6. Repeat simulation multiple times
"""

Nsim = 10000

class Star():
    def __init__(self):
        self.name     = ""
        self.ra       = np.nan
        self.dec      = np.nan
        self.p        = np.nan
        self.perr     = np.nan
        self.pa       = np.nan
        self.paerr    = np.nan
        self.ebv      = np.nan
        self.l        = np.nan
        self.b        = np.nan
        self.p_mean   = np.nan
        self.p_std    = np.nan

    def calcGal(self):
        c = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree, frame='icrs')
        self.l = c.galactic.l.degree
        self.b = c.galactic.b.degree

class Bin():
    def __init__(self):
        self.Ps    = []
        self.Perrs = []

    def calcStat(self):
        if len(self.Ps) > 0:
            self.mean_P = np.mean(self.Ps)
            self.std_P  = np.std(self.Ps)
        else:
            self.mean_P = -0.5
            self.std_P  = -0.5

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

def loadPScat():
    if not os.path.exists('ps1-ebv-4.5kpc.fits'):
        response = urlopen('https://faun.rc.fas.harvard.edu/eschlafly/2dmap/ps1-ebv-4.5kpc.fits')
        CHUNK = 16 * 1024
        with open('ps1-ebv-4.5kpc.fits', 'wb') as f:
            while True:
                chunk = response.read(CHUNK)
                if not chunk:
                    break
                f.write(chunk)

def getHeiles():
    stars = []
    fop = open('Heiles.csv')
    for line in fop.readlines():
        if line.startswith('#'):
            continue
        sl=line.split(',')
        st = Star()
        st.ra    = float(sl[0])
        st.dec   = float(sl[1])
        st.p     = float(sl[5])
        st.perr  = float(sl[6])
        st.pa    = float(sl[7])
        st.paerr = float(sl[8])
        if st.perr > 0 and st.p/st.perr > 3:
            stars.append(st)
            #if len(stars) > 3000:
            #    break
    fop.close()
    return stars

def findEBV(stars):
    good_stars = []
    mp = pyfits.getdata('ps1-ebv-4.5kpc.fits')
    tmp = []
    for st in stars:
        st.calcGal()
        #st.ebv = healpy.get_interp_val(mp['ebv'], np.radians(90. - st.b), np.radians(st.l)) # example is like this
        st.ebv = healpy.get_interp_val(mp['ebv'], st.l, st.b, lonlat=True)
        if not np.isnan(st.ebv):
            tmp.append(st.ebv)
            good_stars.append(st)
    return good_stars

def binEBV(stars):
    bins = np.arange(0.0,3.0,0.03)
    st_binned = map(lambda x: Bin(), np.zeros(shape=len(bins)-1))
    for i in range(len(bins)-1):
        for st in stars:
            if bins[i] < st.ebv <= bins[i+1]:
                st_binned[i].Ps.append(st.p)
                st_binned[i].Perrs.append(st.perr)

    map(lambda x: x.calcStat(), st_binned)
    return st_binned, bins

def plot(st_binned, bins):
    figure = plt.figure(figsize=(12, 4), dpi=150)
    figure.subplots_adjust(hspace=0.1)
    fig = plt.subplot(111)
    plt.xlabel('E (B-V)')
    plt.ylabel('P (%)',fontsize=16)
    
    
    xmean = [0.5*(x+y) for x,y in zip(bins,bins[1:])]
    y = map(lambda x: x.mean_P, st_binned)
    y_p_std = map(lambda x: x.mean_P + x.std_P, st_binned)
    y_n_std = map(lambda x: x.mean_P - x.std_P, st_binned)
    
    plt.scatter(xmean, y,c='r',s=1)
    plt.scatter(xmean, y_p_std,c='b',s=1)
    plt.scatter(xmean, y_n_std,c='b',s=1)
    plt.savefig('EBV_P.png',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

def getUFOs():
    UFOs = []
    hdulist = pyfits.open('../PASIPHAE_estimate/gll_psc_v16.fit')
    header = hdulist[1].header
    tbdata = hdulist[1].data
    mp = pyfits.getdata('ps1-ebv-4.5kpc.fits')
    for td in tbdata:
        fgl3name=td[0]
        ra  = float(td[1])
        dec = float(td[2])
        glon= float(td[3])
        glat= float(td[4])
        CLASS1 = td[73]

        if CLASS1 == '': # and (glat < -40 or glat > 40):
            ufo = Star()
            ufo.name = fgl3name.replace(" ","")
            ufo.ebv = healpy.get_interp_val(mp['ebv'], glon, glat, lonlat=True)
            ufo.ra  = ra
            ufo.dec = dec
            ufo.l   = glon
            ufo.b   = glat
            UFOs.append(ufo)
    return UFOs

def fillNanEBV(UFOs):
    """
    Fill missing EBV by data from SF11
    """
    SF11 = {}
    fop = open('../PASIPHAE_estimate/unidentALLextinct.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        name = sl[0]
        AB   = float(sl[4])
        AV   = float(sl[5])
        SF11[name] = AB - AV
    fop.close()
    for ufo in UFOs:
        if np.isnan(ufo.ebv):
            ufo.ebv = SF11[ufo.name]
    return UFOs

def assignPolDistPar(UFOs, st_binned, bins):
    for ufo in UFOs:
        for i in range(len(bins)-1):
            if bins[i] < ufo.ebv <= bins[i+1]:
                if bins[i+1] < 1.91 and st_binned[i].mean_P > 0 and st_binned[i].std_P > 0:
                    ufo.p_mean = st_binned[i].mean_P
                    ufo.p_std  = st_binned[i].std_P
                else:
                    # assigning very large numbers so we definitely cannot detect it
                    ufo.p_mean = 50
                    ufo.p_std  = 1
    return UFOs

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
        P = 100*np.random.beta(alpha, beta, size=1)[0]
        Ps.append(P)
    return Ps

def simulate(UFOs,Pblaz):
    det_fr = []
    N = len(UFOs)
    FieldPols = [np.random.normal(x.p_mean, x.p_std, size=Nsim) for x in UFOs]
    for i in range(Nsim):
        detect = 0
        Ps = GetRandBlazars(N,Pblaz)
        for m in range(N):
            if FieldPols[m][i] < Ps[m]: #  WARNING: there are negative values here
                detect += 1
        det_fr.append(float(detect)/N)
    print np.mean(det_fr)
    print np.std(det_fr)
        

if __name__ == "__main__":
    Pblaz = ReadBlazarCat()
    '''
    loadPScat()
    stars = getHeiles()
    stars = findEBV(stars)

    output = open('data.pkl', 'wb')
    pickle.dump(stars, output)
    output.close()
    '''
    pkl_file = open('data.pkl', 'rb')
    stars = pickle.load(pkl_file)
    pkl_file.close()
    
    st_binned, bins = binEBV(stars)
    plot(st_binned, bins)
    UFOs = getUFOs()
    UFOs = fillNanEBV(UFOs)
    UFOs = assignPolDistPar(UFOs, st_binned, bins)
    simulate(UFOs,Pblaz)


