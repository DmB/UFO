import math

name=[]
ra=[]
dec=[]
semi_major=[]
semi_minor=[]
flux=[]
var=[]
bl=[]
bl2=[]
name_new=[]
ra_new=[]
dec_new=[]
semi_major_new=[]
semi_minor_new=[]
flux_new=[]
var_new=[]
bl_new=[]

fop = open('U3FGL.csv')
fop2 = open('ranking.csv','w+')

for line in fop.readlines():
    if line.startswith('#'):
        continue
    sl=line.split(',')
    name.append(str(sl[0]))
    ra.append(str(sl[1]))
    dec.append(str(sl[2]))
    semi_major.append(str(sl[8]))
    semi_minor.append(str(sl[9]))
    flux.append(str(sl[18]))
    var.append(str(sl[55]))
for i in range(len(ra)):
    bl.append(float(10.0**22.0*float(var[i])*float(flux[i])/float(math.pi*float(semi_major[i])*float(semi_minor[i]))))    
    bl2.append(float(10.0**22.0*float(var[i])*float(flux[i])/float(math.pi*float(semi_major[i])*float(semi_minor[i]))))

bl2.sort(reverse=True)
for i in range(len(bl2)):
    for j in range(len(bl2)):
        if bl2[i]==bl[j]:
            name_new.append(name[j])
            ra_new.append(ra[j])
            dec_new.append(dec[j])
            semi_major_new.append(semi_major[j])
            semi_minor_new.append(semi_minor[j])
            flux_new.append(flux[j])
            var_new.append(var[j])
            break

fop2.write(str('#Name ranking RA DEC Semi_major Semi_minor Flux Var'))
fop2.write('\n')
order=0
for i in range(len(ra)):
    #bl.append(float(10.0**22.0*float(var[i])*float(flux[i])/float(math.pi*float(semi_major[i])*float(semi_minor[i]))))
    fop2.write(str(name_new[i]))
    fop2.write(',')
    fop2.write(str(bl2[i]))
    fop2.write(',')
    fop2.write(str(ra_new[i]))
    fop2.write(',')
    fop2.write(str(dec_new[i]))
    fop2.write(',')
    fop2.write(str(semi_major_new[i]))
    fop2.write(',')
    fop2.write(str(semi_minor_new[i]))
    fop2.write(',')
    fop2.write(str(flux_new[i]))
    fop2.write(',')
    fop2.write(str(var_new[i]))
    fop2.write('\n')
    order+=1
    if name_new[i] == '3FGL J0221.2+2518':
        print 'parameter=', str(bl2[i]), ' order=',order



fop.close()
fop2.close()