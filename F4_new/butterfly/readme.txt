

following https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/python_tutorial.html

we select 1 - 523 weeks of data

>>> import gt_apps as my_apps
>>> my_apps.filter['evclass'] = 128
>>> my_apps.filter['evtype'] = 3
>>> my_apps.filter['ra'] = 35.3203
>>> my_apps.filter['dec'] = 25.3056
>>> my_apps.filter['rad'] = 15
>>> my_apps.filter['emin'] = 100
>>> my_apps.filter['emax'] = 100000
>>> my_apps.filter['zmax'] = 90
>>> my_apps.filter['tmin'] = 234316801
>>> my_apps.filter['tmax'] = 550022405
>>> my_apps.filter['infile'] = '@events.list'
>>> my_apps.filter['outfile'] = 'UFO_filtered.fits'
>>> my_apps.filter.run()


>>> my_apps.maketime['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
>>> my_apps.maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
>>> my_apps.maketime['roicut'] = 'no'
>>> my_apps.maketime['evfile'] = 'UFO_filtered.fits'
>>> my_apps.maketime['outfile'] = 'UFO_filtered_gti.fits'
>>> my_apps.maketime.run()

>>> my_apps.expCube['evfile'] = 'UFO_filtered_gti.fits'
>>> my_apps.expCube['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
>>> my_apps.expCube['outfile'] = 'UFO_ltCube.fits'
>>> my_apps.expCube['zmax'] = 90
>>> my_apps.expCube['dcostheta'] = 0.025
>>> my_apps.expCube['binsz'] = 1
>>> my_apps.expCube.run()


>>> my_apps.expMap['evfile'] = 'UFO_filtered_gti.fits'
>>> my_apps.expMap['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
>>> my_apps.expMap['expcube'] = 'UFO_ltCube.fits'
>>> my_apps.expMap['outfile'] = 'UFO_expMap.fits'
>>> my_apps.expMap['irfs'] = 'CALDB'
>>> my_apps.expMap['srcrad'] = 25
>>> my_apps.expMap['nlong'] = 120
>>> my_apps.expMap['nlat'] = 120
>>> my_apps.expMap['nenergies'] = 10
>>> my_apps.expMap.run()


>>> my_apps.diffResps['evfile'] = 'UFO_filtered_gti.fits'
>>> my_apps.diffResps['scfile'] = '/home/blinov/ssd/data/lat_spacecraft_merged.fits'
>>> my_apps.diffResps['srcmdl'] = '/home/blinov/Fermi/data/models/RBPLJ0221+2518_model.xml'
>>> my_apps.diffResps['irfs'] = 'CALDB'
>>> my_apps.diffResps.run()

>>> obs = UnbinnedObs('UFO_filtered_gti.fits','/home/blinov/ssd/data/lat_spacecraft_merged.fits',expMap='UFO_expMap.fits',expCube='UFO_ltCube.fits',irfs='CALDB')
>>> like = UnbinnedAnalysis(obs,'/home/blinov/Fermi/data/models/RBPLJ0221+2518_model.xml',optimizer='NewMinuit')

>>> like.tol = 0.0001
>>> likeobj = pyLike.NewMinuit(like.logLike)
>>> like.fit(verbosity=0,covar=True,optObject=likeobj)

>>> print likeobj.getRetCode()

