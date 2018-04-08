#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import scipy.stats as st
import sys, math, copy, random

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

mpl.rcParams['ps.fonttype'] = 42

c = 2.99792458E+18

WLs = []
Frs = []
Fs = []
Ferrs = []
nuFs = []
nuFerrs = []
fitflag = []
fop = open('UFO_1152-0027816.dat')
for line in fop.readlines():
    if line.startswith("#"):
        continue
    sl = line.split()
    WLs.append(float(sl[1]))
    Frs.append(c/float(sl[1]))
    Fs.append(float(sl[2]))
    nuFs.append(float(sl[2]) * float(sl[1]))
    Ferrs.append(float(sl[3]))
    nuFerrs.append(float(sl[3]) * float(sl[1]))
    fitflag.append(sl[4])
fop.close()



x = []
y = []
for i in range(len(WLs)):
    if fitflag[i] != "y":
        continue
    x.append(np.log10(Frs[i]))
    y.append(np.log10(nuFs[i]))

fit = np.polyfit(x, y, 2)
p = np.poly1d(fit)

X = np.arange(13.2,15.2,0.1)
X2 = []
Y = []
for x in X:
    y = p(x)
    X2.append(10**x)
    Y.append(10**y)

figure = plt.figure(figsize=(4, 3.9), dpi=150)
figure.subplots_adjust(hspace=0.1)
fig = plt.subplot(111)
#plt.xlabel(r'$\mathrm{\lambda$ $\mathrm{(A)}$')
plt.xlabel(r'$\mathrm{\nu}$ $\mathrm{(Hz)}$')
plt.ylabel(r'$\mathrm{\nu F_{\nu} (erg/cm^2/s)}$',fontsize=16) #,position=(-0.07,0.5)
#fig.get_yaxis().set_label_coords(-0.07,0.5)

plt.errorbar(Frs,nuFs,yerr=nuFerrs,capthick=0,fmt='o')
plt.plot(X2,Y)
plt.xlim([1e8,1e26])
plt.ylim([10**(-14.3),1e-9])
fig.set_yscale('log')
fig.set_xscale('log')
#plt.xticks([100,200,300], ['100','200','300'])
#plt.yticks([1,10,20,30], ['1','10','20','30'])

plt.savefig('SED.png',bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()

