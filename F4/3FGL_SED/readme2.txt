

following https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/python_tutorial.html

we select 1 - 523 weeks of data

import gt_apps as my_apps
my_apps.filter['evclass'] = 128
my_apps.filter['evtype'] = 3
my_apps.filter['ra'] = 35.3203
my_apps.filter['dec'] = 25.3056
my_apps.filter['rad'] = 15
my_apps.filter['emin'] = 100
my_apps.filter['emax'] = 100000
my_apps.filter['zmax'] = 90
my_apps.filter['tmin'] = 303667469
my_apps.filter['tmax'] = 312941069
my_apps.filter['infile'] = '@events.list'
my_apps.filter['outfile'] = 'UFO_filtered.fits'
my_apps.filter.run()


my_apps.maketime['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
my_apps.maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
my_apps.maketime['roicut'] = 'no'
my_apps.maketime['evfile'] = 'UFO_filtered.fits'
my_apps.maketime['outfile'] = 'UFO_filtered_gti.fits'
my_apps.maketime.run()

my_apps.expCube['evfile'] = 'UFO_filtered_gti.fits'
my_apps.expCube['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
my_apps.expCube['outfile'] = 'UFO_ltCube.fits'
my_apps.expCube['zmax'] = 90
my_apps.expCube['dcostheta'] = 0.025
my_apps.expCube['binsz'] = 1
my_apps.expCube.run()


my_apps.expMap['evfile'] = 'UFO_filtered_gti.fits'
my_apps.expMap['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
my_apps.expMap['expcube'] = 'UFO_ltCube.fits'
my_apps.expMap['outfile'] = 'UFO_expMap.fits'
my_apps.expMap['irfs'] = 'CALDB'
my_apps.expMap['srcrad'] = 25
my_apps.expMap['nlong'] = 120
my_apps.expMap['nlat'] = 120
my_apps.expMap['nenergies'] = 10
my_apps.expMap.run()


my_apps.diffResps['evfile'] = 'UFO_filtered_gti.fits'
my_apps.diffResps['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
my_apps.diffResps['srcmdl'] = '/home/blinov/Fermi/data/models/RBPLJ0221+2518_model.xml'
my_apps.diffResps['irfs'] = 'CALDB'
my_apps.diffResps.run()

import pyLikelihood
from UnbinnedAnalysis import *
obs = UnbinnedObs('UFO_filtered_gti.fits','/home/blinov/ssd/data/lat_spacecraft_merged.fits',expMap='UFO_expMap.fits',expCube='UFO_ltCube.fits',irfs='CALDB')
like = UnbinnedAnalysis(obs,'/home/blinov/Fermi/data/models/RBPLJ0221+2518_model.xml',optimizer='NewMinuit')

like.tol = 0.0001
likeobj = pyLike.NewMinuit(like.logLike)
like.fit(verbosity=0,covar=True,optObject=likeobj)

print likeobj.getRetCode()


#>>> for sn in sourceDetails.keys():
#...    if sourceDetails[sn] < 0:
#...       print sn, sourceDetails[sn]

sourceDetails = {}
for source in like.sourceNames():
    sourceDetails[source] = like.Ts(source)

print sourceDetails

for source,TS in sourceDetails.iteritems():
    print source, TS
    if (TS < 0):
       print "Deleting...", source
       like.deleteSource(source)


like.logLike.writeXml('UFO_fit.xml')

from likeSED import *

inputs = likeInput(like,'RBPLJ0221+2518',nbins=4, model='UFO_fit.xml')
#low_edges = [200.,   427.69,  914.61, 1955.87, 4182.56,  8944.27, 19127.05,  40902.61]
#high_edges = [427.69,914.61, 1955.87, 4182.56, 8944.27, 19127.05, 40902.61, 187049.69]
#low_edges = [100.,   914.61,   4182.56,   19127.05]
#high_edges = [914.61, 4182.56, 19127.05, 100000.00]

low_edges = [100.,   1800,   10000.0]
high_edges = [1800, 10000.0, 100000.00]

inputs.customBins(low_edges,high_edges)


inputs.srcMaxE() #to find most energetic photon
inputs.fullFit(CoVar=True)
sed = likeSED(inputs)
sed.getECent()
sed.fitBands(ftol=1e-3,opt='Minuit')#default ftol=1e-3,opt='NewMinuit'

sed.Plot()
