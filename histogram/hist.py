#!/usr/bin/python

import numpy as np

from scipy import stats
import matplotlib.pyplot as plt

p = []
perr = []
fop = open('../F3/field.dat')
for line in fop.readlines():
    if line.startswith("#"):
        continue
    sl = line.split()
    p.append( float(sl[8]) * 100)
    perr.append( float(sl[9]) * 100 )
fop.close()

figure = plt.figure(figsize=(6, 6.0), dpi=150)
figure.subplots_adjust(hspace=0.1)
fig = plt.subplot(111)

n, bins, patches = plt.hist( p, 10, normed=False, histtype='bar',color='#b3b3b3')
plt.ylabel("N",fontsize=18)
plt.xlabel("Polarization degree (%)",fontsize=18)

ax = plt.gca()

ax.set_yticks([0,5,10,15,20])
#ax.set_yticklabels([0])
fig.tick_params(axis='both', which='major', labelsize=18)

#plt.arrow( 5, 2.8, 0, -1.2, fc="k", ec="k", width=0.1,length_includes_head=False, head_width=0.3, head_length=0.3)
#ax.arrow(4.6, 2.8, 0.4, -1.2, head_width=0.1, head_length=0.3, fc='k', ec='k')
#plt.annotate('UGSC', xy=(4.3, 3), xytext=(4.1, 3),fontsize=20)
plt.annotate('UGSC', xy=(5, 1), xytext=(4.1, 3), arrowprops=dict(arrowstyle='->'),fontsize=20)

plt.savefig('p_hist.eps',bbox_inches='tight')


