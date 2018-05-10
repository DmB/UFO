#!/usr/bin/python

import math
import pyfits

class Obj():
    def __init__(self):
        self.name       = None
        self.ra         = None
        self.dec        = None
        self.flux       = None
        self.var        = None
        self.semi_major = None
        self.semi_minor = None

    def calcRank(self):
        self.rank = self.var * 1e6 * self.flux / (math.pi * self.semi_major * self.semi_minor )

    def __lt__(self, other):
        return self.rank < other.rank


ufos = []

hdulist = pyfits.open('../PASIPHAE_estimate/gll_psc_v16.fit')
header = hdulist[1].header
tbdata = hdulist[1].data

for td in tbdata:
    CLASS1 = td[73]
    if CLASS1 == '':
        o = Obj()
        o.name = td[0]
        o.ra   = td[1]
        o.dec  = td[2]
        o.semi_major = td[8]
        o.semi_minor = td[9]
        o.flux = td[35] + td[39] + td[43] + td[47] + td[51]
        o.var  = td[55]
        o.calcRank()
        ufos.append(o)

ufos.sort(reverse=True)

fop2 = open('ranking.csv','w')
fop2.write(str('#Name ranking RA DEC Semi_major Semi_minor Flux Var'))
fop2.write('\n')
order=0
for ufo in ufos:
    fop2.write(str(ufo.name))
    fop2.write(',')
    fop2.write(str(ufo.rank))
    fop2.write(',')
    fop2.write(str(ufo.ra))
    fop2.write(',')
    fop2.write(str(ufo.dec))
    fop2.write(',')
    fop2.write(str(ufo.semi_major))
    fop2.write(',')
    fop2.write(str(ufo.semi_minor))
    fop2.write(',')
    fop2.write(str(ufo.flux))
    fop2.write(',')
    fop2.write(str(ufo.var))
    fop2.write('\n')
    order+=1
    if ufo.name == '3FGL J0221.2+2518':
        print 'F4 J0221.2+2518 parameter=', str(ufo.rank), ' order=',order
    elif ufo.name == '3FGL J1848.6+3232':
        print 'F1 J1848.6+3232 parameter=', str(ufo.rank), ' order=',order
    elif ufo.name == '3FGL J0336.1+7500':
        print 'F2 J0336.1+7500 parameter=', str(ufo.rank), ' order=',order
    elif ufo.name == '3FGL J0419.1+6636':
        print 'F3 J0419.1+6636 parameter=', str(ufo.rank), ' order=',order

fop2.close()
