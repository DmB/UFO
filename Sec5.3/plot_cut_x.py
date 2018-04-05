#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

class Obs():
    def __init__(self):
        self.JD = None
        self.W1    = None
        self.W1err = None
        self.W2    = None
        self.W2err = None

def ReadData():
    data = []
    fop = open('w1w2_multiepoch.dat')
    for line in fop.readlines():
        if line.startswith("#"):
            continue
        sl = line.split()
        o = Obs()
        o.JD = float(sl[0]) - 55000
        o.W1    = float(sl[1])
        o.W1err = float(sl[2])
        o.W2    = float(sl[3])
        o.W2err = float(sl[4])
        data.append(o)
    fop.close()
    return data

def AverGroups(data):
    return None

def Plot(data):
    JDs = map(lambda x: x.JD, data)
    W1s = map(lambda x: x.W1, data)
    W1errs = map(lambda x: x.W1err, data)
    W2s = map(lambda x: x.W2, data)
    W2errs = map(lambda x: x.W2err, data)
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
    ax1.errorbar(JDs,W1s,yerr=W1errs,linestyle='None')
    ax1.scatter(JDs,W1s)
    ax2.errorbar(JDs,W1s,yerr=W1errs,linestyle='None')
    ax2.scatter(JDs,W1s)
    ax3.errorbar(JDs,W1s,yerr=W1errs,linestyle='None')
    ax3.scatter(JDs,W1s)
    ax1.set_xlim([223, 224])
    ax2.set_xlim([412, 414])
    ax3.set_xlim([587, 588.5])
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax1.tick_params(labelright='off')
    ax2.yaxis.set_ticks_position('none') 
    ax3.yaxis.tick_right()
    
    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)
    ax2.plot((1-d,1+d), (-d,+d), **kwargs)
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs)
    
    
    kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
    ax3.plot((-d,+d), (1-d,1+d), **kwargs)
    ax3.plot((-d,+d), (-d,+d), **kwargs)

    ax2.set_xlabel("JD (-2455000)")
    ax1.set_ylabel("W1 (mag)")
    plt.savefig('wise.png',bbox_tight="True")

if __name__ == "__main__":
    data = ReadData()
    averages = AverGroups(data)
    Plot(data)


