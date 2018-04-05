import numpy as np
import matplotlib.pyplot as plt
import sys, os, shutil
import glob
import aplpy
from astropy import coordinates
from astropy import units
import matplotlib.cm as cm          
#from general_functions import stokes_parameters_csv


def coordinates_os_stars(path):
    
    ra, dec, Ps, angles = [], [], [], []
    fop = open(path)
    fop.readline()
    for line in fop.readlines():
        sl = line.split()
        c = coordinates.ICRSCoordinates(ra=sl[1],dec=sl[2],unit=(units.hour,units.degree))
        x = c.ra.degree
        y = c.dec.degree
        P = float(sl[5])
        theta = float(sl[7])
        ra.append(x), dec.append(y), angles.append(theta), Ps.append(P)
    fop.close()    
    # in order to show the scale
    ra.append(91.91), dec.append(22.05), angles.append(90), Ps.append(1.0)
    return ra, dec, Ps, angles

def segments_on_map(fitsfile, ra, dec, Ps, pas, scale = 1000):
        '''
        Plot polarization segments on fits file.
        Input: 
        fitsfile: string, name of fits file
        ra: list of ra
        dec: list of dec
        pas: list of angles in radians. Angles are with respect to y axis of image.
             You need to rotate the real angle to have angles with respect to north.
        savename: string, name of plot to be saved
        coord: coordinate system of fits image
        scale: float/int, determines how long the segments will be in pixels
        '''
        fig = aplpy.FITSFigure(fitsfile,figsize = (9,9))
        fig.show_grayscale(invert=True)
        fig.add_grid()
        fig.grid.show()
        
        linelist1 = []
        kukloi_ra, kukloi_dec, kukloi_polosi = [], [], []
        for iv in range(len(pas)):
            xpix,ypix=fig.world2pixel(ra[iv],dec[iv])
            linelength_half= scale*Ps[iv]

            y=[ypix-linelength_half*np.cos(pas[iv]),ypix+linelength_half*np.cos(pas[iv])]
            x=[xpix+linelength_half*np.sin(pas[iv]),xpix-linelength_half*np.sin(pas[iv])]
            x_world,y_world=fig.pixel2world([x[0],x[1]],[y[0],y[1]])
            line=np.array([x_world,y_world])
            linelist1.append(line)
        
        fig.show_lines(linelist1, layer='line', color='r')
        fig.add_label(91.858,22.08, 'fs1', color='r')
        fig.add_label(91.8707,22.0966, 'fs2', color='r')
        fig.add_label(91.889,22.11805555555, 'fs3', color='r')
        fig.add_label(91.836,22.10222222222, 'fs4', color='r')
        fig.add_label(91.91,22.052,'1\%', color='r')
        
        fig.save('segments_SDSS_field.eps')
        return



if __name__ == '__main__':
    # Load in the coordinates of the stars
    path = 'IGRJ06074+2205.txt'
    eikona_fits = 'dss.06.07.26.60+22.05.48.0.fits'
    ras_all, decs_all, Ps, angles = coordinates_os_stars (path)  
    #q, u, sq, su, polosis = stokes_parameters_csv(path)
    angles = np.array(angles)
    angles = np.radians(angles)
    ra = ras_all
    dec = decs_all
    # Plot polarization segments on DSS 
    coordang = 0.
    pas_forplot = angles + coordang
   
    segments_on_map( eikona_fits, ra, dec, Ps, pas_forplot, scale = 19)

    
