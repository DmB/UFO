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
    #fop.readline()
    for line in fop.readlines():
        if line.startswith('#'):
            continue
        sl = line.split()
        c = coordinates.ICRSCoordinates(ra=sl[0],dec=sl[1],unit=(units.degree,units.degree))
        x = c.ra.degree
        y = c.dec.degree
        P = 100*float(sl[8])
        theta = float(sl[10])
        ra.append(x), dec.append(y), angles.append(theta), Ps.append(P)
    fop.close() 
    
    """ Convert angles from range [-pi,pi] or [-pi/2,pi/2]  
        to [0,pi).
    """
    for aa in range(len(angles)):
        if angles[aa] < 0:
            angles[aa] = np.pi - abs(angles[aa])
        if angles[aa] == np.pi:
            angles[aa] = 0
        angles[aa] = angles[aa]*180.0/np.pi
    
   
    # in order to show the scale
    ra.append(35.18), dec.append(25.17), angles.append(90), Ps.append(1.0)
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
        fig.show_regions('ds9.reg')
        fig.add_label(35.38918,25.175,'UFS', color='black')
        fig.add_label(35.34918,25.21,'RG', color='black')
        fig.add_label(35.18,25.175,'1%', color='r')
        
        fig.save('segments_SDSS_field.eps')
        return



if __name__ == '__main__':
    # Load in the coordinates of the stars
    path = 'final.dat'
    eikona_fits = 'dss.02.21.16.+25.18.20.fits'
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

    
