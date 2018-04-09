#!/usr/bin/python

import pyfits
import sys
import os
from astropy import coordinates
from astropy import units
try:
    from urllib.request import urlopen # Python 3
except ImportError:
    from urllib2 import urlopen # Python 2


if not os.path.exists('gll_psc_v16.fit'):
    response = urlopen('https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4yr_catalog/gll_psc_v16.fit')
    CHUNK = 16 * 1024
    with open('gll_psc_v16.fit', 'wb') as f:
        while True:
            chunk = response.read(CHUNK)
            if not chunk:
                break
            f.write(chunk)

hdulist = pyfits.open('gll_psc_v16.fit')

header = hdulist[1].header
tbdata = hdulist[1].data

fop = open('unidentPASIPHAE.dat','w')
fop.write('# 3FGL name   RA   DEC\n')
for td in tbdata:
    fgl3name=td[0]
    ra  = float(td[1])
    dec = float(td[2])
    glon=float(td[3])
    glat=float(td[4])
    sum_flux = td[35] + td[39] + td[43] + td[47] + td[51]
    l_err    = td[36][0] + td[40][0] + td[44][0] + td[48][0] + td[52][0]
    u_err    = td[36][1] + td[40][1] + td[44][1] + td[48][1] + td[52][1]
    power_law_index = td[30]
    Variability_Index = td[55]
    CLASS1 = td[73]
    ASSOC1 = td[74]
    ASSOC2 = td[75]
    if CLASS1 == '' and (glat < -40 or glat > 40):
        #print fgl3name.lstrip('3FGL '), ",", glat, ",", ASSOC1, ",", ASSOC2
        fop.write(fgl3name + '     ' + str(ra) + ' ' + str(dec) + '\n')

fop.close()
