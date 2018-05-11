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

class Obs():
    def __init__(self):
        self.lmbda = None # Angstroms
        self.freq  = None
        self.SFD = None # spectral energy flux density erg/cm2/s/A
        self.SFDerr = None
        self.nuFnu = None # erg/cm2/s
        self.nuFnuerr = None
        self.Lumin = None # erg/s
        self.Luminerr = None

    def calcNuFnu(self):
        SFD = ufloat(self.SFD,self.SFDerr)
        nuFnu = self.lmbda * SFD
        self.nuFnu    = nuFnu.n
        self.nuFnuerr = nuFnu.s

    def calcLumin(self):
        Lumin = 4*np.pi*(cm_in_Mpc * Dl)**2 * ufloat(self.nuFnu,self.nuFnuerr)
        self.Lumin    = Lumin.n
        self.Luminerr = Lumin.s

    def calcFreq(self):
        self.freq  = const.c.cgs.value * 1e8 / self.lmbda # Hz
        

def ReadData(src):
    obs = []
    if src == "UFO":
        fn = 'UFO_1152-0027816.dat'
    elif src == "3FGL":
        fn = "3FGLJ0221.2+2518.dat"
    fop = open(fn)
    for line in fop.readlines():
        if line.startswith("#"):
            continue

        sl = line.split()
        o = Obs()
        o.lmbda  = float(sl[1])
        o.SFD    = float(sl[2])
        o.SFDerr = float(sl[3])
        o.calcNuFnu()
        o.calcLumin()
        o.calcFreq()
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
    figure = plt.figure(figsize=(3, 2.56), dpi=200)
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
        ax_lumin.errorbar(ufo.freq,ufo.Lumin,yerr=ufo.Luminerr,markersize=2,capthick=0,fmt='o',color='k')
    for ferm in Fermi:
        ax_lumin.errorbar(ferm.freq,ferm.Lumin,yerr=ferm.Luminerr,markersize=2,capthick=0,fmt='o',color='r')
    
    plt.xlim([1e9,1e27])
    ax_lumin.set_ylim([1e41,1e44])
    
    ax_flux.set_ylabel(r'$\mathrm{\nu F_{\nu} (erg/s/cm^2)}$',fontsize=11)
    
    plt.savefig('SED.eps',bbox_inches='tight')

if __name__ == "__main__":
    UFO = ReadData("UFO")
    Fermi = ReadData("3FGL")
    plot(UFO,Fermi)
    Lumin = 4*np.pi*(cm_in_Mpc * Dl)**2 * ufloat(3.78455e-12,8.77048e-13)
    print 'Lumin. in 100 Mev - 100 Gev range is ',Lumin.n,' +- ', Lumin.s

