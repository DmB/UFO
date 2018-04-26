#!/usr/bin/python

import pyfits, healpy
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy import units as u
import pickle

class Star():
    def __init__(self):
        self.ra       = np.nan
        self.dec      = np.nan
        self.p        = np.nan
        self.perr     = np.nan
        self.pa       = np.nan
        self.paerr    = np.nan
        self.ebv      = np.nan

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

if __name__ == "__main__":
    
    mp = pyfits.getdata('ps1-ebv-4.5kpc.fits')
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
    '''
    st_binned, bins = binEBV(stars)
    plot(st_binned, bins)




