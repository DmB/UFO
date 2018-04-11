#!/usr/bin/python

import random

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def FitIdentMags():
    Rs = []
    usno = np.genfromtxt('USNO_ident.csv', delimiter=',',skip_header=1)
    for i in range(len(usno)):
        R = usno[i,40]
        if not np.isnan(R):
            Rs.append(R)

    plt.figure()
    weights = np.ones_like(Rs)/float(len(Rs))
    n, bins, patches = plt.hist( Rs, 50, weights=weights, normed=False, histtype='bar',color='green', alpha=0.75)
    plt.xlabel("R [mag]")
    plt.ylabel("N")
    plt.savefig('dist_mags_ident.png',bbox_inches='tight')
    return

if __name__ == "__main__": 
    FitIdentMags()
