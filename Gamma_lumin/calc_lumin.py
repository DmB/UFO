#!/usr/bin/python

'''
calculates gamma-ray luminosities from photon fluxes and photon indexes
according to http://arxiv.org/pdf/1403.4961.pdf
'''

import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import numpy as np

ev2erg = 1.602176565e-12 # ergs in ev
flux     = 3.32160e-10 #photon/cm^2/MeV/s 1-100GeV from 3FGL
flux_err = 7.72456e-11 #photon/cm^2/MeV/s from 3FGL
Gamma    = 1.93567 #photon ind from 3FGL
z = 0.062 # our data

def calcLumin():
    """
    E - energy flux
    """
    if Gamma == 1:
        E = flux * (10**3 - 1.0)/(np.log(10)**3) * ev2erg * 1e6
        E_err = flux_err * (10**3 - 1.0)/(np.log(10)**3) * ev2erg * 1e6
    elif Gamma == 2:
        E = flux * (np.log(10)**3)/(1.0 - 1e-3) * ev2erg * 1e6
        E_err = flux_err * (np.log(10)**3)/(1.0 - 1e-3) * ev2erg * 1e6
    else:
        E = flux * 100 * (Gamma - 1.0) / (Gamma - 2.0) * (1 - 10**(3*(2.0-Gamma)))/(1 - 10**(3*(1.0-Gamma))) * ev2erg * 1e6
        E_err = flux_err * 100 * (Gamma - 1.0) / (Gamma - 2.0) * (1 - 10**(3*(2.0-Gamma)))/(1 - 10**(3*(1.0-Gamma))) * ev2erg * 1e6

    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308) # Planck Collab 2015, Paper XIII
    Lumin = 4 * np.pi * (cosmo.luminosity_distance(z).cgs.value)**2 * (1.0 + z)**(Gamma-2.0) * E
    L_err = 4 * np.pi * (cosmo.luminosity_distance(z).cgs.value)**2 * (1.0 + z)**(Gamma-2.0) * E_err
    return Lumin, L_err

if __name__ == "__main__":
    Lumin, L_err = calcLumin()
    print "Lumin = ",Lumin,"+-",L_err
