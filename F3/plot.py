#!/usr/bin/python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class Star():
    def __init__(self):
        self.q    = None
        self.qerr = None
        self.u    = None
        self.uerr = None

def ReadData(fn):
    data = []
    fop = open(fn)
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        s = Star()
        s.q     = float(sl[4])
        s.qerr = float(sl[5])
        s.u     = float(sl[6])
        s.uerr = float(sl[7])
        data.append(s)
    fop.close()
    return data

def plot(data,ufo):
    Qs = map(lambda x: x.q , data)
    Qerrs = map(lambda x: x.qerr , data)
    Us = map(lambda x: x.u, data)
    Uerrs = map(lambda x: x.uerr , data)
    
    
    ufo_Qs = map(lambda x: x.q , ufo)
    ufo_Qerrs = map(lambda x: x.qerr , ufo)
    ufo_Us = map(lambda x: x.u, ufo)
    ufo_Uerrs = map(lambda x: x.uerr , ufo)
    
    figure = plt.figure(figsize=(5, 5), dpi=150)
    figure.subplots_adjust(hspace=0.1)
    fig = plt.subplot(111)
    
    plt.xlabel('Q/I',fontsize=16)
    plt.ylabel('U/I',fontsize=16)
    plt.errorbar(Qs, Us, xerr=Qerrs, yerr=Uerrs, capthick=0, fmt='o', linewidth=1)
    plt.errorbar(ufo_Qs, ufo_Us, xerr=ufo_Qerrs, yerr=ufo_Uerrs, capthick=0, fmt='o', linewidth=1, c='r')
    plt.axvline(0,linestyle='dashed',c='k',linewidth=1)
    plt.axhline(0,linestyle='dashed',c='k',linewidth=1)
    plt.savefig('Q_U.eps',bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

if __name__ == "__main__":
    data = ReadData('field.dat')
    ufo  = ReadData('ufo.dat')
    plot(data,ufo)


