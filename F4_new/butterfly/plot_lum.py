#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
from uncertainties import ufloat
import numpy as np

cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
z = 0.06169
cm_in_Mpc = 3.08568e24
Dl = cosmo.luminosity_distance(z).value # in Mpc
h = 4.135667662e-15 # plank const in eV s

class Obs():
    def __init__(self):
        self.name  = None
        self.lmbda = None # Angstroms
        self.freq  = None
        self.SFD = None # spectral energy flux density erg/cm2/s/A
        self.SFDerr = None
        self.nuFnu = None # erg/cm2/s
        self.nuFnuerr = None
        self.Lumin = None # erg/s
        self.Luminerr = None
        self.energy = None  # GeV
        self.EnFl_up = None  # erg/cm2/s
        self.EnFl_lo = None  # erg/cm2/s

    def calcNuFnu(self):
        SFD = ufloat(self.SFD,self.SFDerr)
        nuFnu = self.lmbda * SFD
        self.nuFnu    = nuFnu.n
        self.nuFnuerr = nuFnu.s

    def calcLumin(self):
        Lumin = 4*np.pi*(cm_in_Mpc * Dl)**2 * ufloat(self.nuFnu,self.nuFnuerr)
        self.Lumin    = Lumin.n
        self.Luminerr = Lumin.s
        
    def calcLumin2(self):
        self.Lumin_G_up = 4*np.pi*(cm_in_Mpc * Dl)**2 * self.EnFl_up
        self.Lumin_G_lo = 4*np.pi*(cm_in_Mpc * Dl)**2 * self.EnFl_lo

    def calcFreq(self):
        self.freq  = const.c.cgs.value * 1e8 / self.lmbda # Hz

    def EnToFreq(self):
        self.freq  = self.energy * 1e9 / h
        

def ReadData():
    obs = []
    fn = '../UFO_1152-0027816.dat'
    fop = open(fn)
    for line in fop.readlines():
        if line.startswith("#"):
            continue

        sl = line.split()
        o = Obs()
        o.name   = sl[0]
        o.lmbda  = float(sl[1])
        o.calcFreq()
        o.SFD    = float(sl[2])
        o.SFDerr = float(sl[3])
        o.calcNuFnu()
        o.calcLumin()
        obs.append(o)
    fop.close()
    return obs

def ReadDataGamma():
    obs = []
    fn = "butterfly.dat"
    fop = open(fn)
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        o = Obs()
        o.energy = float(sl[0])
        o.EnToFreq()
        o.EnFl_up = float(sl[1])
        o.EnFl_lo = float(sl[2])
        o.calcLumin2()
        obs.append(o)
    fop.close()
    return obs

def Lum2Flux(Lum):
    return Lum / (4*np.pi*(cm_in_Mpc * Dl)**2)

def axis_callback(ax_lumin):
    """
    Update second axis according with first axis.
    """
    global ax_flux
    y1, y2 = ax_lumin.get_ylim()
    ax_flux.set_ylim(Lum2Flux(y1), Lum2Flux(y2))
    ax_flux.figure.canvas.draw()

def plot(UFO,Fermi):
    global ax_flux
    figure = plt.figure(figsize=(4, 3), dpi=200)
    figure.subplots_adjust(hspace=0.1)
    fig = plt.subplot(111)
    plt.xlabel(r'$\mathrm{\nu}$ $\mathrm{(Hz)}$')
    plt.ylabel(r'$\mathrm{\nu L_{\nu} (erg/s)}$',fontsize=11) #,position=(-0.07,0.5)

    ax_lumin = plt.gca()
    ax_flux = ax_lumin.twinx()
    ax_lumin.callbacks.connect("ylim_changed", axis_callback)

    fig.set_xscale('log')
    ax_lumin.set_yscale('log')
    ax_flux.set_yscale('log')

    for ufo in UFO:
        print ufo.freq,ufo.Lumin
        if ufo.name.startswith("GALEX"):
            color = '#8b00ff'
        elif ufo.name.startswith("GAIA"):
            color = '#00cc00'
        elif ufo.name.startswith("SLOAN"):
            color = '#ffa500'
        elif ufo.name.startswith("2MASS"):
            color = '#EC0ED6'
        elif ufo.name.startswith("WISE"):
            color = '#8b4513'
            
        ax_lumin.errorbar(ufo.freq,ufo.Lumin,yerr=ufo.Luminerr,markersize=2,capthick=0,fmt='o',color=color)
    
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='#8b4513', label="WISE")
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='#EC0ED6', label="2MASS")
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='#ffa500', label="SDSS")
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='#00cc00', label="Gaia")
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='#8b00ff', label="GALEX")
    ax_lumin.errorbar(0,0,yerr=0,markersize=2,capthick=0,fmt='o',color='r', label="Fermi")
    
    ax_lumin.text(1e13,5e43,'UFSC')
    ax_lumin.text(2e19,5e43,'3FGLJ0221.2+2518')

    #for ferm in Fermi:
    freqs      = map(lambda x: x.freq, Fermi)
    Lumin_G_up = map(lambda x: x.Lumin_G_up, Fermi)
    Lumin_G_lo = map(lambda x: x.Lumin_G_lo, Fermi)
    ax_lumin.plot(freqs,Lumin_G_up,markersize=0.5,color='r')
    ax_lumin.plot(freqs,Lumin_G_lo,markersize=0.5,color='r')
    
    ax_lumin.plot([Fermi[0].freq,Fermi[0].freq],[Fermi[0].Lumin_G_up,Fermi[0].Lumin_G_lo],markersize=0.5,color='r')
    ax_lumin.plot([Fermi[-1].freq,Fermi[-1].freq],[Fermi[-1].Lumin_G_up,Fermi[-1].Lumin_G_lo],markersize=0.5,color='r')
    
    plt.xlim([1e10,1e27])
    ax_lumin.set_ylim([1e41,1e44])
    ax_lumin.legend(loc=8,ncol=3,fontsize=9)

    
    ax_flux.set_ylabel(r'$\mathrm{\nu F_{\nu} (erg/s/cm^2)}$',fontsize=11)
    
    plt.savefig('SED.eps',bbox_inches='tight')

if __name__ == "__main__":
    UFO = ReadData()
    Fermi = ReadDataGamma()
    plot(UFO,Fermi)
    Lumin = 4*np.pi*(cm_in_Mpc * Dl)**2 * ufloat(3.78455e-12,8.77048e-13)
    print 'Lumin. in 100 Mev - 100 Gev range is ',Lumin.n,' +- ', Lumin.s

