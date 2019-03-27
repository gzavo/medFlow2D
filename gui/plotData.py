#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014, All Right Reserved, Gábor Závodszky, gabor@zavodszky.com
#
# This source is subject to the
# Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.
# Please see the License.txt file for more information.
# All other rights reserved.
# You should have received a copy of the license along with this
# work in the License.txt file.
# If not, see <http://creativecommons.org/licenses/by-nc-nd/4.0/>.

import os

import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

import tkFileDialog
import zipfile

u=None
v=None
mag=None
cX=None
cY=None
rho=None
pres=None
conc=None
fx=None
fy=None
u2=None
u3=None

cqX=None
cqY=None
qRho=None

a=None
f=None

t=None
s=None

ax=None
ax2=None

canvas=None
root = None

slices = 450
x=0
y=0    
axis=False  # 0 - x; 1 -y
sliceVal=None

qstep=10.0
xRes=2.89337e-5

def on_key_event(event):
    global canvas, toolbar
    print('you pressed %s'%event.key)
    key_press_handler(event, canvas, toolbar)

def chooseDirToSave():
    fileName = tkFileDialog.asksaveasfile(mode='w', defaultextension='txt')

    if(fileName):        
        for i in range(len(s)):
            fileName.write("%f %f\n" % (t[i],s[i]))
        fileName.close()

def changeSlice():
    global axis
    axis = not axis

def changeSliceVal():
    global sliceVal, slices
    slices = round(float( sliceVal.get() ))
    ax2.set_xlabel('Pixel')
    ax2.set_ylabel('Pixel', rotation=0)
    ax2.yaxis.set_label_coords(-0.04, 1.2)
    ax2.set_title('Slice')

def changexRes():
    global xResval, xRes, dudy, u3
    xRes = float(xResVal.get())
    dudy=(np.divide(u3,xRes))*(-1)

def _quit():
    global root
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate

def plotData(slices, data):
    if axis:
        plotDataY(slices, data)
    else:
        plotDataX(slices, data)


def plotDataX(slices, data):
    global f, root, canvas, ax, ax2, t, s
    t=cY[slices, :]    
    s=data[slices, :]
    ax.clear()
    ax.plot(t,s)
    
    if (data==u).all():
        ttl='x velocity'
        ylbl='v_x [m/s]'
    elif (data==v).all():
        ttl='y velocity'
        ylbl='v_y [m/s]'
    elif (data==mag).all():
        ttl='Velocity'
        ylbl='Velocity [m/s]'
    elif (data==pres).all():
        ttl='Pressure'
        ylbl='p [Pa]'
    elif (data==conc).all():
        ttl='Platelets concentration'
        ylbl='Concentration'    
    elif (data==fx).all():
        ttl='Margination force_x'
        ylbl='Force_x [N]' 
    elif (data==fy).all():
        ttl='Margination force_y'
        ylbl='Force_y [N]' 
    elif (data==dudy).all():
        ttl='Shear velocity'
        ylbl='du/dy [1/s]' 
    else:
        ttl='Sg wrong'
        ylbl='sg wrong'

    ax.set_xlabel('Pixel')
    ax.set_title('%s'%ttl)
    ax.set_ylabel('%s'%ylbl)

    ax2.clear()
    ax2.pcolor(cqX,cqY,qRho)
    ax2.plot([slices, slices], [0, y-1], 'y-', lw=2)
    ax2.set_xlabel('Pixel')
    ax2.set_ylabel('Pixel', rotation=0)
    ax2.yaxis.set_label_coords(-0.04, 1.2)
    ax2.set_title('Slice')
    canvas.draw()

def plotDataY(slices, data):
    global f, root, canvas, ax, ax2, t, s
    t=cX[:, y-slices]    
    s=data[:, y-slices]
    ax.clear()    
    ax.plot(t,s)

    if (data==u).all():
        ttl='x velocity'
        ylbl='v_x [m/s]'
    elif (data==v).all():
        ttl='y velocity'
        ylbl='v_y [m/s]'
    elif (data==mag).all():
        ttl='Velocity'
        ylbl='Velocity [m/s]'
    elif (data==pres).all():
        ttl='Pressure'
        ylbl='p [Pa]'
    elif (data==conc).all():
        ttl='Platelets concentration'
        ylbl='Concentration'    
    elif (data==fx).all():
        ttl='Margination force_x'
        ylbl='Force_x [N]' 
    elif (data==fy).all():
        ttl='Margination force_y'
        ylbl='Force_y [N]' 
    elif (data==dudy).all():
        ttl='Shear velocity'
        ylbl='du/dy [1/s]' 
    else:
        ttl='Sg wrong'
        ylbl='sg wrong'
        
    ax.set_xlabel('Pixel')
    ax.set_title('%s'%ttl)
    ax.set_ylabel('%s'%ylbl)

    ax2.clear()    
    ax2.pcolor(cqX,cqY,qRho)
    ax2.plot([0, x-1], [slices, slices], 'y-', lw=2)
    ax2.set_xlabel('Pixel')
    ax2.set_ylabel('Pixel', rotation=0)
    ax2.yaxis.set_label_coords(-0.04, 1.2)
    ax2.set_title('Slice')
    canvas.draw()

def plotContour(contourdata,caxislabel):
    global f, root, canvas, ax, cX, cY
    ax.clear() 
    if caxislabel == 0:
        cxlabel='Pixel'
        cylabel='Pixel'
        ctitle='Velocity contour 2D'
    elif caxislabel == 1:
        cxlabel='Pixel'
        cylabel='Pixel'
        ctitle='Pressure contour 2D'
    elif caxislabel == 2:
        cxlabel='Pixel'
        cylabel='Pixel'
        ctitle='Fluid age contour2D'
    elif caxislabel == 3:
        cxlabel='Pixel'
        cylabel='Pixel'
        ctitle='Fluid age'
    
    cdmin=contourdata.min() + 0.000001
    cdmax=contourdata.max()
    levels = np.linspace(cdmin, cdmax, 20)

    if caxislabel == 3:
        ax.contourf(np.swapaxes(contourdata, 0,1)[::-1,:], levels)
        #ax.set_ylim(ax.get_ylim()[::-1])
    else:
        ax.contour(np.swapaxes(contourdata, 0,1)[::-1,:], levels)
        #ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel('%s'%cxlabel)
        ax.set_title('%s'%ctitle)
        ax.set_ylabel('%s'%cylabel)
    canvas.draw()

    
def loadDataFile(fileName):
    global u,v,mag,cX,cY,rho,pres,conc,fx,fy,x,y,u2,u3,dudy, cqX, cqY, qRho

    #fin = open(fileName)
    zip_fin = zipfile.ZipFile(fileName)
    zip_name = os.path.basename(fileName)
    text_vel_name = zip_name[:-4]

    fin = zip_fin.open(text_vel_name)

    x,y = map(int, fin.readline().split())

    u=np.zeros((x,y))
    v=np.zeros((x,y))
    mag=np.zeros((x,y))
    cX=np.zeros((x,y))
    cY=np.zeros((x,y))
    rho=np.zeros((x,y))
    pres=np.zeros((x,y))
    conc=np.zeros((x,y))
    fx=np.zeros((x,y))
    fy=np.zeros((x,y))
    u2=np.zeros((x,y))
    u3=np.zeros((x,y))
    dudy=np.zeros((x,y))
    
    cqX=np.zeros((int(x/qstep),int(y/qstep)))
    cqY=np.zeros((int(x/qstep),int(y/qstep)))
    qRho=np.zeros((int(x/qstep),int(y/qstep)))
    
    for j in range(0, y):
        for i in range(0,x):
            tu, tv, trho, tpres, tconc, tfx, tfy = map(float, fin.readline().split())
            u[i,j]=tu
            v[i,j]=-tv
            mag[i,j]=np.sqrt(tu**2+tv**2)
            pres[i,j]=tpres
            rho[i,j]=trho
            conc[i,j]=tconc
            cX[i,j]=i
            cY[i,j]=y-j
            fx[i,j]=tfx
            fy[i,j]=-tfy 

    fin.close()
    zip_fin.close()
    
    u2=np.diff(u) 
    u3[:,:-1] = u2
    u3=np.absolute(u3)
    dudy=(np.divide(u3,xRes))*(-1)  
            
    for j in range(0, int(y/qstep)):
        for i in range(0,int(x/qstep)):
            cqX[i,j] = cX[i*qstep, j*qstep]
            cqY[i,j] = cY[i*qstep, j*qstep]
            if (rho[i*qstep, j*qstep])==0:
                qRho[i,j]=0.0
            else:    
                qRho[i,j] = rho[i*qstep, j*qstep] - 1.0

if __name__ == "__main__":

    if(len(sys.argv) != 2):
        print "Argument needs to be a valid data file."
        sys.exit(-1)
    
    root = Tk.Tk()
    root.wm_title("Data plotter")

    f = Figure(figsize=(5,4), dpi=100, tight_layout=True, frameon=True)
        
    loadDataFile(sys.argv[1])
    
    slices = round(x/2.0)
    
    t=cY[slices, :]    
    s=u[slices, :]

    ax = f.add_subplot(121) 
    ax.plot(t,s)
    ax.set_xlabel('Pixel')
    ax.set_ylabel('Velocity [m/s]')
    ax.set_title('Velocity')

    ax2 = f.add_subplot(122)
    ax2.pcolor(cqX,cqY,qRho)
    ax2.plot([slices, slices], [0, y-1], 'y-', lw=2)
    ax2.set_xlabel('Pixel')
    ax2.set_ylabel('Pixel', rotation=0)
    ax2.yaxis.set_label_coords(-0.04, 1.2)
    ax2.set_title('Slice')

    cax = f.gca()
    cax.set_aspect('equal')


    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg( canvas, root )
    toolbar.update()
    canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

    canvas.mpl_connect('key_press_event', on_key_event)

    button1 = Tk.Button(master=root, text='VelX', command= lambda: plotData(slices, u) )
    button1.pack(side=Tk.LEFT)
    button2 = Tk.Button(master=root, text='VelY', command= lambda: plotData(slices, v) )
    button2.pack(side=Tk.LEFT)
    button3 = Tk.Button(master=root, text='Vel', command= lambda: plotData(slices, mag) )
    button3.pack(side=Tk.LEFT)
    button4 = Tk.Button(master=root, text='Pres', command= lambda: plotData(slices, pres) )
    button4.pack(side=Tk.LEFT)
    button5 = Tk.Button(master=root, text='Conc', command= lambda: plotData(slices, conc))
    button5.pack(side=Tk.LEFT)
    button6 = Tk.Button(master=root, text='fX', command= lambda: plotData(slices, fx) )
    button6.pack(side=Tk.LEFT)
    button7 = Tk.Button(master=root, text='fY', command= lambda: plotData(slices, fy) )
    button7.pack(side=Tk.LEFT)
    button8 = Tk.Button(master=root, text='du/dy', command= lambda: plotData(slices, dudy) )
    button8.pack(side=Tk.LEFT)
    button9 = Tk.Button(master=root, text='2D vel. contour', command= lambda: plotContour(mag,0) )
    button9.pack(side=Tk.LEFT)
    button10 = Tk.Button(master=root, text='2D pres. contour', command= lambda: plotContour(pres,1) )
    button10.pack(side=Tk.LEFT)
    button11 = Tk.Button(master=root, text='fluid age', command= lambda: plotContour(conc,3) )
    button11.pack(side=Tk.LEFT)
    button12 = Tk.Button(master=root, text='f.age contour', command= lambda: plotContour(conc,2) )
    button12.pack(side=Tk.LEFT)
    button13 = Tk.Button(master=root, text="Save data",command= lambda: chooseDirToSave())
    button13.pack(side=Tk.LEFT)


    button = Tk.Button(master=root, text='Quit', command=_quit)
    button.pack(side=Tk.RIGHT)

    sliceVal = Tk.StringVar()
    entry1 = Tk.Entry(master=root, text="100", textvariable=sliceVal)
    sliceVal.set(str(slices))
    entry1.pack(side=Tk.RIGHT)

    button14 = Tk.Button(master=root, text='Set coord', command= lambda: changeSliceVal() )
    button14.pack(side=Tk.RIGHT)

    button15 = Tk.Button(master=root, text='x/y', command= lambda: changeSlice() )
    button15.pack(side=Tk.RIGHT)
    
    xResVal = Tk.StringVar()
    entry2 = Tk.Entry(master=root, text="2.89337e-5", textvariable=xResVal)
    xResVal.set(str(xRes))
    entry2.pack(side=Tk.RIGHT)

    button16 = Tk.Button(master=root, text="Set xRes",command= lambda: changexRes())
    button16.pack(side=Tk.RIGHT)

    root.wm_iconbitmap('@gui/main.xbm')
    Tk.mainloop()