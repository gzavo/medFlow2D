#!/usr/bin/python
# -*- coding: utf-8 -*-

#   Calculate several measures of a fractal structure based on the residence times
#   alive - single array with residence time values of the number-of-alive-particles
#   tmin - calculate from this time (no need to calculate until some particle have left the domain)
#   tmax - calculate up to this time value (after a while the few number of alive particles spoil the exponential behaviour)

import sys

from numba import jit

import Tkinter as Tk
from Tkinter import *
import tkFileDialog
import numpy as np
import matplotlib.pyplot as plt
import time

@jit
def freeEnergy(data,t,beta):
    search = np.array(data > t, dtype=int)
    maximum = len(search)
    sign = search[0]
    li = 0.0
    betaf = 0.0
    for i in range(1, maximum):
        li = li+1
        if (search[i] != sign) or (i == maximum-1 and search[i] == 1):
            sign = search[i]
            if (sign == 0) or (i == maximum-1 and search[i] ==1):
                betaf = betaf + li**beta
            li = 0
    betaf = np.log(betaf)
    return betaf

def Solver():
    global alive, ey, efit, escapeRate, bFb, lyapunov, metricEntropy, infDimension
    global tesc, logesc, xd, ey, eyfit_eval, rsq, betapl
    # Computing escape rate, at Beta=1
    start = time.time()
    ey = []
    for k in xd:
        ey= np.append(ey, freeEnergy(alive, k, 1.0))
    efit = np.polyfit(xd, ey, 1)
    escapeRate = -efit[0]
    strescaperate = str('%10.8f'%escapeRate)
    L18.config(text=strescaperate)

    #Compute Lyupanov exponent, metric entropy and information dimension
    bFb = np.zeros((1,nstep))
    Fb = np.zeros((1,nstep))

    for number in range(1, nstep+1):
        yd = []
        betal = number/(nstep/2.0)
        for k in xd: # yd=arrayfun in matlab
            yd = np.append(yd, -freeEnergy(alive, k, betal))
        fitd = np.polyfit(xd, yd, 1)
        bFb[0, number-1] = fitd[0]
        Fb[0, number-1] = bFb[0, number-1]/(number/(nstep/2.0))
    dbFb = np.diff(bFb)/wstep
    dFb = np.diff(Fb)/wstep
    lyapunov = dbFb[0, round(nstep/2.0)-1]
    metricEntropy = dFb[0, round(nstep/2.0)-1]
    infDimension = metricEntropy/lyapunov

    strlyapunov = str('%10.8f'%lyapunov)
    L20.config(text=strlyapunov)

    strmetricEntropy = str('%10.8f'%metricEntropy)
    L22.config(text=strmetricEntropy)

    strinfDimension = str('%10.8f'%infDimension)
    L26.config(text=strinfDimension)

    tTot = np.max(alive)
    slices = 200
    tesc = np.linspace(0, tTot, slices)
    esc = np.zeros((slices,1))
    for t in range(1,slices+1):
        esc[t-1] = np.sum(alive > (tTot/slices)*(t-1))
    logesc = np.log(esc)

    betapl = np.linspace(2.0/nstep, 2.0, nstep)
    betapl = betapl.reshape(1, len(betapl))

    eyfit_eval = np.polyval(efit, xd)
    yresid = ey-eyfit_eval
    SSresid = np.sum(np.square(yresid))
    SStotal = len(ey) * np.var(ey)
    rsq = 1-(SSresid/SStotal)

    strrsq = str('%10.5f'%rsq)
    L24.config(text=strrsq)

    nearD1 = (np.abs(betapl- infDimension)).argmin() #give the index of the closest value to zero 
    
    #Linear interpolation between two points
    if betapl[0, nearD1] == infDimension:
        dy = bFb[0, nearD1]
    elif betapl[0, nearD1] > infDimension:
        if nearD1 != 0:
            x1 = betapl[0, nearD1-1]
            x2 = betapl[0, nearD1]
            y1 = bFb[0, nearD1-1]
            y2 = bFb[0, nearD1]
            dy = y1 + (y2-y1)*(infDimension-x1)/(x2-x1) # tgalfa = (y2-y1) / (x2-x1) = (y-y1) / (x-x1)  solve x=infDimension   
        else:
            dy = '####'
            print 'There aro no enough points to interpolate at x=D1'
    else:
        x1 = betapl[0, nearD1]
        x2 = betapl[0, nearD1+1]
        y1 = bFb[0, nearD1]
        y2 = bFb[0, nearD1+1]
        dy = y1 + (y2-y1)*(infDimension-x1)/(x2-x1)

    dyD1 = 0-dy # distance from y=0

    strdyD1 = str('%10.5f'%dyD1)
    L28.config(text=strdyD1)
    L16.config(text="**Finished**")

    end = time.time()
    print end - start

def getEscapeRate(slices, tTot):
    esc = np.zeros((slices,1))
    for t in range(0,slices):
        esc[t] = np.sum(alive > (tTot/slices)*t)

    return esc

def Plotfile(gr):

    if gr == 0:
        plt.plot(alive0)
        plt.xlabel('Traced particle')
        plt.ylabel('Residence time')
        plt.title('All particle traced')
    elif gr == 1:
        plt.plot(alive)
        plt.xlabel('Traced particle')
        plt.ylabel('Residence time')
        plt.title('%d-%d particle traced' % (cfrom+1, cto))
    elif gr == 2:

        f, (ax1, ax2) = plt.subplots(2)
        ax1.plot(tesc,logesc[:,0], label='Kappa from raw data')
        ax1.plot(xd, ey, 'bo', label='Kappa from free energy')
        ax1.plot(xd, eyfit_eval, 'r', linewidth=2.0, label='Fitted line')
        text1 = '$R^2:$ %10.5f' %rsq
        ax1.text(0.70, 0.85, text1, transform=ax1.transAxes, bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
        legend = ax1.legend()
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        for label in legend.get_texts():
            label.set_fontsize('small')
        ax1.set(xlabel='time [s]', ylabel='log(number of particles)')
        ax1.set_title('Escape rate (Kappa)')
        ax2.plot(betapl[0,:], bFb[0,:])

        ax2.set(xlabel='Beta', ylabel='Beta F(Beta)')
        ax2.set_title('Beta F(Beta)')
        text2 = ' $n:$   %d \n $\kappa:$ %10.5f \n $h:$ %10.5f \n $K_1:$ %10.5f \n $D_1:$ %10.5f' % (nstep, escapeRate, lyapunov, metricEntropy, infDimension)
        ax2.text(0.05, 0.45, text2, transform=ax2.transAxes, bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
        ax2.plot(infDimension, 0, 'r+', ms=20)
        plt.tight_layout()

    plt.show()


def AskFile():
    global alive0, alive_at_end

    alive0=[]
    L1.config(text="**Wait**")
    ask = tkFileDialog.askopenfile()
    alive0 = np.loadtxt(ask, usecols=(0,))
    npoints = len(alive0)
    strstuckpoints = str((-1 == alive0).sum())
    strstillactive = str((0 == alive0).sum())
    strallpoints = str(npoints)
    L3.config(text=strallpoints)
    L7.config(text=strstillactive)
    L5.config(text=strstuckpoints)

    max_alive = np.max(alive0)
    alive_at_end = len(alive0 == 0)
    alive0[alive0==0] = max_alive
    L1.config(text="**Loaded**")

def changeInput():
    global cfrom, cto, tmin, tmax, nstep, wstep, xd, alive

    cfrom = int(cfromstr.get())-1
    cto = int(ctostr.get())
    tmin = float(tminstr.get())
    tmax = float(tmaxstr.get())
    nstep = int(nstepstr.get())
    wstep = (float(tmax)-float(tmin))/nstep
    xd = np.linspace(tmin, tmax, nstep)
    alive = []
    alive = alive0[cfrom:cto]

    L14.config(text="**Done**")

    print 'Changed:', cfrom+1, cto, tmin, tmax, nstep


if __name__ == "__main__":

    w = Tk()
    w.wm_title("Fractal dimensions calculator")
    w.geometry("350x400")

    alive = []
    alive_at_end = 0
    tmin = 6.77
    tmax = 12.01
    nstep = 50
    cfrom = 316000
    cto = 420000
    wstep = (float(tmax)-float(tmin))/nstep
    xd = np.linspace(tmin, tmax, nstep)
    search = []


    strescaperate = ''
    strlyapunov = ''
    strmetricEntropy = ''
    strinfDimension = ''
    strrsq = ''
    strdyD1 = ''
    strallpoints = ''
    strstuckpoints = ''
    strstillactive = ''

    L0 = Label(w, text='Open and plot the complete file:')
    L0.grid(row=0, column=0, sticky='W', columnspan=3)

    L1 = Label(w, text='**Open file**')
    L1.grid(row=1, column=3, sticky='W', columnspan=1)


    b1 = Button(w, text="Open", command = lambda: AskFile() )
    b1.grid(row=1, column=2)

    b2 = Button(w, text="Plot", command = lambda: Plotfile(0) )
    b2.grid(row=1, column=4)

    L2 = Label(w, text="All points:")
    L2.grid(row=2, column=2, sticky='E')
    L3 = Label(w, text=strallpoints)
    L3.grid(row=2, column=3, sticky='W')

    L4 = Label(w, text="Stuck points:")
    L4.grid(row=3, column=2, sticky='E')
    L5 = Label(w, text=strstuckpoints)
    L5.grid(row=3, column=3, sticky='W')

    L6 = Label(w, text="Still active:")
    L6.grid(row=4, column=2, sticky='E')
    L7 = Label(w, text=strstillactive)
    L7.grid(row=4, column=3, sticky='W')

    L8 = Label(w, text='Selecting chaotic part of the file:')
    L8.grid(row=5, column=0, sticky='W', columnspan=3)

    L9 = Label(w, text='From:')
    L9.grid(row=6, column=1, sticky='E')
    cfromstr = StringVar()
    e1 = Entry(w, text="75000", width=10, justify="center", textvariable=cfromstr)
    cfromstr.set(str(cfrom))
    e1.grid(row=6, column=2)

    L10 = Label(w, text='To:')
    L10.grid(row=6, column=3, sticky='E')
    ctostr = StringVar()
    e2 = Entry(w, text="210000", width=10, justify="center", textvariable=ctostr)
    ctostr.set(str(cto))
    e2.grid(row=6, column=4)


    L11 = Label(w, text='tmin:')
    L11.grid(row=7, column=1, sticky='E')
    tminstr = StringVar()
    e3 = Entry(w, text="1.04", width=10, justify="center", textvariable=tminstr)
    tminstr.set(str(tmin))
    e3.grid(row=7, column=2)

    L12 = Label(w, text='tmax:')
    L12.grid(row=7, column=3, sticky='E')
    tmaxstr = StringVar()
    e4 = Entry(w, text="2.09", width=10, justify="center", textvariable=tmaxstr)
    tmaxstr.set(str(tmax))
    e4.grid(row=7, column=4)

    L13 = Label(w, text='nstep:')
    L13.grid(row=8, column=1, sticky='E')
    nstepstr = StringVar()
    e5 = Entry(w, text="50", width=10, justify="center", textvariable=nstepstr)
    nstepstr.set(str(nstep))
    e5.grid(row=8, column=2)

    b3 = Button(w, text="Set", command = lambda: changeInput() )
    b3.grid(row=9, column=2)
    L14 = Label(w, text='**Click set**')
    L14.grid(row=9, column=3, sticky='W', columnspan=1)

    b4 = Button(w, text="Plot", command = lambda: Plotfile(1) )
    b4.grid(row=9, column=4)

    L15 = Label(w, text='Calculate dimensions:')
    L15.grid(row=10, column=0, sticky='W', columnspan=3)

    b5 = Button(w, text="Solve", command = lambda: Solver() )
    b5.grid(row=11, column=4)

    L16 = Label(w, text='**Click Solve**')
    L16.grid(row=11, column=3, sticky='W', columnspan=1)

    L17 = Label(w, text="Escape Rate:" )
    L17.grid(row=12, column=1, sticky='E')
    L18 = Label(w, text=strescaperate)
    L18.grid(row=12, column=2, sticky='W')

    L19 = Label(w, text="Lyapunov:")
    L19.grid(row=13, column=1, sticky='E')
    L20 = Label(w, text=strlyapunov)
    L20.grid(row=13, column=2, sticky='W')

    b6 = Button(w, text="Plot results", command = lambda: Plotfile(2) )
    b6.grid(row=13, column=4)

    L21 = Label(w, text="Metric Entropy:")
    L21.grid(row=14, column=1, sticky='E')
    L22 = Label(w, text=strmetricEntropy)
    L22.grid(row=14, column=2, sticky='W')

    L23 = Label(w, text=('R^2: '))
    L23.grid(row=14, column=3, sticky='E')
    L24 = Label(w, text=strrsq)
    L24.grid(row=14, column=4, sticky='W')

    L25 = Label(w, text="Inf Dimension:")
    L25.grid(row=15, column=1, sticky='E')
    L26 = Label(w, text=strinfDimension)
    L26.grid(row=15, column=2, sticky='W')

    L27 = Label(w, text=('dy(x=D1,y=0): ')) #Distance from the D1 at y=0
    L27.grid(row=15, column=3, sticky='E')
    L28 = Label(w, text=strdyD1)
    L28.grid(row=15, column=4, sticky='W')

    w.mainloop()
