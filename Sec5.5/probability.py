#import numpy  as np
import math

'''
#Use this https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html to find x^2
#http://exoplanetsdigest.com/2012/02/06/how-many-sigma/   to find Q
#Q = 1-probability you want. eg for probability 95% (2sigma) Q = 0.05
#5sigma --> s=30.0389   4sigma --> s=19.3352   3sigma --> s=11.6182   2sigma --> s=5.9914   1sigma --> s=2.2788
#d = dimension of data (for x,y data d=2)
'''


s = 11.6182  #the x^2 value
semi_maj_2sigma = 5.73 #value in arcmins
semi_min_2sigma = 4.8
s_x = semi_maj_2sigma/5.9914**0.5
s_y = semi_min_2sigma/5.9914**0.5 #5.991 = x^2 value for 2sigma
semi_maj_new = float(s_x) * s**0.5
semi_min_new = float(s_y) * s**0.5

    
#print 'semi_maj_new = ', semi_maj_new, 'arcmins'
#print 'semi_maj_new = ', semi_maj_new/60.0, 'degrees'
#print 'semi_min_new = ', semi_min_new, 'arcmins'
#print 'semi_min_new = ', semi_min_new/60.0, 'degrees'


el_area = math.pi*semi_maj_new*semi_min_new
el_area_deg = el_area/3600

SDSS_area_deg =9376.0
quasars = 526356.0
quasars_mag = 8767   #mag <= 17.59
quasars_z = 11 #z <= 
SDSS_density = quasars/SDSS_area_deg
SDSS_density_mag = quasars_mag/SDSS_area_deg
SDSS_density_z = quasars_z/SDSS_area_deg


probability = SDSS_density*el_area_deg
probability_mag = SDSS_density_mag*el_area_deg
probability_z = SDSS_density_z*el_area_deg
probability_comb = probability_mag*probability_z


print 'probability of quasars = ', probability
print 'probability with right mag = ', probability_mag
print 'probability with right z = ', probability_z
print 'probability of UFO = ', probability_comb

