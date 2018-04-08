import numpy as np



u=[]
uerr=[]
q=[]
qerr=[]
p=[]
perr=[]
n=[]
ra=[]
dec=[]
angle=[]
fop = open('final.dat')
#fop2 = open('centroid.dat','w+')
for line in fop.readlines():
    if line.startswith('#'):
        continue
    sl = line.split()
    #q.append(float(sl[4]))
    #qerr.append(float(sl[5]))
    #u.append(float(sl[6]))
    #uerr.append(float(sl[7]))
    #n.append(float(sl[2]))
    p.append(float(sl[8]))
    perr.append(float(sl[9]))
    #ra.append(float(sl[0]))
    #dec.append(float(sl[1]))
    #angle.append(float(sl[10]))


xos2=0
sm2=0
for i in range(len(p)):
    xos2 += (p[i]/perr[i]**2)
    sm2 += 1.0/perr[i]**2
wmean = xos2/sm2
wdev = np.sqrt(1.0/sm2)

print 'Weighted mean = ',wmean,' Error = ',wdev  



# CORRECT BEFORE USE!!!!!!!!

##print 'p', p
##print 'perr', p_err
#fop2.write('# RA   DEC   star   p  some bullshit    some bullshit   sigmas from average p   angle>50deg from average   degree seperation')
#fop2.write('\n')
#for i in range(len(q)):
    #fop2.write(str(ra[i]))
    #fop2.write('   ')
    #fop2.write(str(dec[i]))
    #fop2.write('   ')
    #fop2.write(str(int(n[i])))
    #fop2.write('   ')
    #fop2.write(str(p[i]))
    #fop2.write('   ')
    #fop2.write(str(float(abs(p[i]-p_av)/perr[i])))  #Using as sigma perr of source 
    #fop2.write('   ')
    #fop2.write(str(float(abs(p[i]-p_av)/p_err)))  #Using as sigma p_err 
    #fop2.write('   ')
    #fop2.write(str(float(abs(p[i]-p_av)/dev_p)))  #Using as sigma dev_p
    #if (angle[i] - a_av) > 0.872665:
        #fop2.write('   ')
        #fop2.write('YES')
        #fop2.write('   ')
        #if (angle[i] - a_av)*57.2958<=90:
            #fop2.write(str((angle[i] - a_av)*57.2958))
        #else:
            #fop2.write(str(180-(angle[i] - a_av)*57.2958))
    #fop2.write('\n')








fop.close()
#fop2.close()