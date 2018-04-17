#!/usr/bin/python

"""
Simulation of possibility do descriminate blazar from stars
variant 2

1. We found A_B and A_V for real UFO fields see ../PASIPHAE_estimate/get_ext.py
2. Pmax from  Hiltner, 1956 Pmax <= 9 E(B - V)(%/mag)
   Also possible to repeat Fosalba et al. 2002  and find STD(P) within those bins in fig. 3
3. We draw a random blazar from the flux limited unbiased sample of RoboPol
4. using appendix from D Blinov et al., MNRAS, 457(2):2252-2262, 2016 we transform intrinsic p0 and m0 to random observed pr
5. if Pmax from step 2 and pr related as pr >= Pmax we consider blazar to be detectable
6. Repeat simulation multiple times
"""

import random
import numpy as np

Nsim  = 10000 # number of trials in MC
N_UFO = 1010 # number of UFOs in 3FGL
N_BL  = int(round(0.5659 * N_UFO)) # expected number of blazars among UFOs taking into account that there are 1717 identified blazars in 3FGL

def ReadExt():
    Pmax = []
    fop = open('unidentALLextinct.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        Pmax.append(float(sl[3]))
    fop.close()
    return Pmax

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

if __name__ == "__main__":
    Pmax = ReadExt()
    Pblaz = ReadBlazarCat()
    all_detections = []
    for i in range(Nsim):
        Ps = GetRandBlazars(N_BL,Pblaz)
        detections = [ps >= pmax for pmax, ps in zip(random.sample(Pmax,N_BL),Ps)].count(True)
        all_detections.append(detections)

    print "Mean number of detection is:", np.mean(all_detections)," out of ",str(N_BL), " expected blazars"
    print "std of detection is:", np.std(all_detections)



