#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as plt
#star   ra   dec  nobs wq qerr  wu  uerr p perr angle angleerr p/sp pdeb
#f=np.genfromtxt('ovro_sample.txt',dtype=None, names=True)
f=np.genfromtxt('plus.dat',dtype=None, names=True)
q=f['wq']
eq=f['qerr']
eu=f['uerr']
u=f['wu']
p=f['p']

sq=[]
su=[]
esq=[]
esu=[]
for pp,ppp in enumerate(p):
    if ppp>0.04:
        print q[pp],u[pp]
        su.append(u[pp])
        sq.append(q[pp])
        esu.append(eu[pp])
        esq.append(eq[pp])


test=[100*mm for mm in p if mm<0.04]
p2=[np.sqrt(qq**2+uu**2) for qq,uu in zip(q,u)]
p=[100*mm for mm in p]
print np.median(p),np.mean(p),np.std(p)
print np.median(test),np.mean(test),np.std(test)
p3=[np.sqrt(qq**2+uu**2) for qq,uu in zip(sq,su)]
plt.rcParams["figure.figsize"] = [7,7]
plt.figure()
#plt.subplot(121)

plt.arrow( 5, 0.1, 0, -0.06, fc="k", ec="k", width=0.01,length_includes_head=True, head_width=0.2, head_length=0.02 ,overhang=0.3)
#plt.annotate('Candidate', xy=(5, 0.1), xytext=(4.68, 0.14),fontsize=14)
plt.annotate('UGSC', xy=(5, 0.1), xytext=(4.7, 0.12),fontsize=20)
plt.hist(p,histtype='step',normed=True)#,label='from file')
plt.ylabel(r"$\mathrm{Probability}$ $\mathrm{density}$",fontsize=18)
plt.xlabel(r"$\mathrm{Polarization}$ $\mathrm{degree}$ $\mathrm{(\%)}$",fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.hist(p2,histtype='step',label='from q,u')
#plt.hist(p3,histtype='step',label='special')
#plt.legend(['Field stars','Candidate source'])
#plt.show()
#plt.subplot(122)
#plt.errorbar(q,u,yerr=eu,xerr=eq,fmt='^',label='Field stars')
#plt.errorbar(sq,su,yerr=esu,xerr=esq,fmt='^',label='Candidate source')
#plt.xlabel(r"$\mathrm{Q/I}$",fontsize=18)
#plt.ylabel(r"$\mathrm{U/I}$",fontsize=18)
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)
#plt.legend( prop={'size': 15})
#plt.tight_layout()
#plt.show()
plt.savefig('p_hist.eps')


