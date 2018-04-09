#!/usr/bin/python

from urllib import urlretrieve
from astropy import coordinates
from astropy import units
import time
import random

prefix = "http://ned.ipac.caltech.edu/cgi-bin/calc?in_csys=Equatorial&in_equinox=J2000.0&obs_epoch=2016.0&lon="
suffix = "d&pa=0.0&out_csys=Equatorial&out_equinox=J2000.0"

names = []
ra = []
dec = []

fop = open("unidentPASIPHAE.dat")
fop.readline()
for line in fop.readlines():
   if line.startswith("#"):
       continue
   sl = line.split()
   names.append(sl[0]+sl[1])
   st_ra = float(sl[2])
   st_dec = float(sl[3])
   try:  # handle different versions of astropy
       c = coordinates.ICRSCoordinates(ra=st_ra,dec=st_dec,unit=(units.hour,units.degree))
   except:
       c = coordinates.SkyCoord(ra=st_ra*units.degree, dec=st_dec*units.degree, frame='icrs')
   ra.append(c.ra.degree)
   dec.append(c.dec.degree)
fop.close()

def parse():
   fop = open('html.html')
   for line in fop.readlines():
      if line.startswith("Landolt R"):
         A_R = line.split()[3]
      if line.startswith("Landolt B"):
         A_B = line.split()[3]
      if line.startswith("Landolt V"):
         A_V = line.split()[3]
   fop.close()
   return A_B,A_V,A_R

fop = open('unidentPASIPHAEextinct.dat','w')
fop.write("# 3FGL               RA             DEC        PD[frac]  AB   AV   AR\n")
for i in range(len(names)):
   url = prefix + str(round(ra[i],4)) + "d&lat=" + str(round(dec[i],4)) + suffix
   urlretrieve(url, 'html.html')
   A_B,A_V,A_R = parse()
   Pmax = 0.09 * (float(A_B) - float(A_V))
   fop.write(names[i] + ' ' + str(ra[i]) + ' ' + str(dec[i]) + ' ' + str(Pmax) + ' ' + str(A_B) + ' ' + str(A_V) + ' ' + str(A_R) + "\n")
   time.sleep(random.randint(3,7))
fop.close()

