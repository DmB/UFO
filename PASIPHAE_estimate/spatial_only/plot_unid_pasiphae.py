#!/usr/bin/python

import pyfits
import sys
import os
from astropy import coordinates
from astropy import units
import matplotlib.pyplot as plt
import numpy as np


hdulist = pyfits.open('../gll_psc_v16.fit')

header = hdulist[1].header
tbdata = hdulist[1].data

glongs = []
glats  = []

for td in tbdata:
    fgl3name=td[0]
    ra  = float(td[1])
    dec = float(td[2])
    glon=float(td[3])
    if glon > 180:
        glon = glon - 360
    glat=float(td[4])
    sum_flux = td[35] + td[39] + td[43] + td[47] + td[51]
    l_err    = td[36][0] + td[40][0] + td[44][0] + td[48][0] + td[52][0]
    u_err    = td[36][1] + td[40][1] + td[44][1] + td[48][1] + td[52][1]
    power_law_index = td[30]
    Variability_Index = td[55]
    CLASS1 = td[73]
    ASSOC1 = td[74]
    ASSOC2 = td[75]
    if CLASS1 == '':# and (glat < -40 or glat > 40):
        glongs.append(glon)
        glats.append(glat)



#r = 90
#theta = np.arange(0,2*np.pi,0.1)
#x = np.array([0,np.pi/3,np.pi/3,0,0])
x = np.array([              0,          np.pi,     np.pi,                    0,                     -np.pi,             -np.pi,                 0])
y = np.array([np.radians(50.),np.radians(50.),    np.radians(89.999),  np.pi/2-0.000001,    np.radians(89.999),  np.radians(50.), np.radians(50.)])


x2 = np.array([              0,          np.pi,     np.pi,                    0,                     -np.pi,             -np.pi,                 0])
y2 = np.array([np.radians(-50.),np.radians(-50.),    np.radians(-89.999),  -np.pi/2+0.000001,    np.radians(-89.999),  np.radians(-50.), np.radians(-50.)])

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection="aitoff")#, axisbg ='LightCyan')
ax.grid(True)

ax.plot(x,y,c='#ffc0cb')

ax.fill_between(x,y,where=y>=0, color='#ffc0cb')

ax.plot(x2,y2,c='#ffc0cb')

ax.fill_between(x2,y2,where=y2<=0, color='#ffc0cb')

ax.scatter(np.radians(glongs),np.radians(glats),zorder=4)  # convert degrees to radians



plt.savefig('Pasiphae_cover.png')



